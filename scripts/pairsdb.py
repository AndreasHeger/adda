#! /bin/env python
################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_kamilah.py 2869 2010-03-03 10:20:13Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""

:Author: Andreas Heger
:Release: $Id: pairsdb.py 2869 2010-03-03 10:20:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Blast all-vs-all computation

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----
"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil
import fileinput, collections, gzip, time

import hashlib, base64

from ruffus import *
import csv
import sqlite3

import FastaIterator
import IOTools
import Pipeline as P
import Experiment as E
import Logfile
import Intervals

##------------------------------------------------------------
def computeHID ( sequence ):
    """returns a hash identifier for a sequence.
    """
    
    # do the encryption
    h = hashlib.md5(sequence).digest()
    # map to printable letters: hid has length 22, so the padded '=' are
    # truncated. You have to add them, if you ever want to decode,
    # but who would do such a thing :=)

    r = base64.encodestring(h)[0:22]
    
    # finally substitute some characters:
    # '/' for '_', so we have legal file names
    # '[' for '+' and ']' for '=' for internet-applications
    
    hid = r.replace('/', '_')
    hid = hid.replace('+', '[')
    hid = hid.replace('=', ']')
                                                
    return hid

@files((( 'uniref50.fasta.gz', 'nrdb50.fasta' ),)) 
def buildNrdb50( infile, outfile ):
    '''build nrdb50
    
    Renumber seqences.'''
    
    outf_fasta = IOTools.openFile( outfile, "w" )
    outf_table = IOTools.openFile( outfile + ".tsv", "w" )
    outf_table.write("nid\tpid\thid\tdescription\tval\tn\ttaxon\trepid\n" )

    rx = re.compile( "(\S+) (.*) n=(\d+) Tax=(.*) RepID=(\S+)" )

    nid = 1
    for entry in FastaIterator.iterate( IOTools.openFile( infile )):
        outf_fasta.write(">%i\n%s\n" % (nid, entry.sequence ) )
        pid, description, cluster_size, taxon, repid = rx.match( entry.title ).groups()
        hid = computeHID( entry.sequence )
        outf_table.write( "\t".join( (str(nid), pid, hid, description, cluster_size, taxon, repid)) + "\n" )
        nid += 1

    outf_fasta.close()
    outf_table.close()

@transform( buildNrdb50, suffix(".fasta"), ".asnb" )
def buildMask( infile, outfile ):
    '''build seg mask for protein sequences.'''

    to_cluster = True

    statement = '''
    segmasker -in %(infile)s
              -infmt fasta 
              -parse_seqids 
              -outfmt maskinfo_asn1_bin 
              -out %(outfile)s
    >& %(outfile)s.log
    '''
    P.run()

@files( ( (None, "scop.fasta"), ) )
def downloadSCOP( infile, outfile ):
    '''download the latest scop sequence set (< 40% identical)'''
    
    statement = '''
    wget -O %(outfile)s "http://astral.berkeley.edu/seq.cgi?get=scopdom-seqres-gd-sel-gs-bib;ver=1.75;item=seqs;cut=40"
    '''
    
    P.run()

@transform( downloadSCOP, suffix('.fasta'), '.links.gz' )
def mapSCOP( infile, outfile ):
    '''map scop against sequence database.
    '''

    to_cluster = True
    max_evalue = 0.00001
    num_results = 100
    mask = 21
    # upper case is critical, otherwise traceback fails!?
    matrix = "BLOSUM50"
    gop = 12
    gep = 2
    dbname = "/tmp/blast/nrdb50"
    num_jobs = 8
    
    job_options = '-pe dedicated %i -R y' % num_jobs

    statement = '''
    /ifs/devel/andreas/cgat/run.py --log=%(outfile)s.log
    'blastp
       -query %(infile)s
       -db %(dbname)s
       -evalue %(max_evalue)f
       -num_alignments %(num_results)i
       -num_descriptions %(num_results)i
       -db_soft_mask %(mask)i
       -matrix %(matrix)s
       -gapopen %(gop)i
       -gapextend %(gep)i
       -num_threads %(num_jobs)i
       -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq"
    | python /ifs/devel/andreas/cgat/blast2table.py 
        --alignment-format=blocks
    | gzip
    > %(outfile)s';
    checkpoint;
    echo "#//" | gzip >> %(outfile)s
    '''

    P.run()

@transform( mapSCOP, suffix('.links.gz'), 
            add_inputs( downloadSCOP ),
            '.domains.gz' )
def buildSCOPDomains( infiles, outfile ):
    '''reconcile mapped domains into a single domain file.

    * fragments are removed - a domain must map at least 90%
      of its length.

    * domains overlapping on the same sequence with the same
      superfamily classification are merged.
    '''
    
    linksfile, fastafile = infiles

    # filtering criteria
    min_coverage = 0.9
    # only take first four fold classes
    classes = 'abcd'

    rx = re.compile('(\S+)\s(\S+)\s(.*)' )
    id2class = {}
    with IOTools.openFile( fastafile ) as inf:
        for x in FastaIterator.iterate( inf ):
            pid, cls, description = rx.match(x.title).groups()
            id2class[pid] = (cls, len(x.sequence) )
            
    E.info('read mappings for %i sequences' % len(id2class))
    counter = E.Counter()

    with IOTools.openFile( linksfile ) as inf:
        nid2domains = collections.defaultdict( list )
        ndomains = 0
        for line in inf:
            if line.startswith('query_nid'): continue
            if line.startswith('#'): continue
            counter.links += 1
            
            domain_id, nid, evalue, domain_start, domain_end, sbjct_start, sbjct_end, \
                block_sizes, domain_starts, sbjct_starts, \
                bitscore, pid = line[:-1].split()
            
            nid, domain_start, domain_end, sbjct_start, sbjct_end = map(int, \
                                                                       ( nid, domain_start, domain_end, sbjct_start, sbjct_end ))

            family, length = id2class[domain_id]

            cls, fold, superfamily, family = family.split('.')
            if cls not in classes: continue
            if float(domain_end - domain_start) / length < min_coverage: continue
            counter.unmerged_domains += 1
            superfamily = '00%c%03i%03i' % (cls, int(fold), int(superfamily))

            nid2domains[nid].append( (superfamily, sbjct_start, sbjct_end ) )

        counter.sequences = len(nid2domains)

    E.info( 'merging %i domains in %i sequences' % (counter.unmerged_domains, counter.sequences))

    outf = IOTools.openFile( outfile, 'w' )
    outf.write('nid\tstart\tend\tfamily\n')
    for nid, dd in sorted(nid2domains.iteritems()):
        for family, domains in itertools.groupby( dd, key = lambda x: x[0] ):
            unmerged_domains = [ (x[1],x[2]) for x in domains ]
            merged_domains = Intervals.combine( unmerged_domains )
            for start, end in merged_domains:
                counter.domains += 1
                outf.write( '%i\t%i\t%i\t%s\n' % (nid, start, end, family ) )
    outf.close()

    E.info( counter )

@merge( (buildNrdb50, buildMask), "nrdb50.log" )
def prepareDatabase( infiles, outfile ):
    '''prepare the blast database.'''

    fastafile, maskfile = infiles
    to_cluster = True
    statement = '''
    makeblastdb 
            -in %(fastafile)s
            -dbtype prot 
            -parse_seqids
            -mask_data %(maskfile)s
            -out nrdb50
            -title "Uniref Protein Database"
    >& %(outfile)s
    '''
    P.run()
    
@follows( mkdir( "blast.dir") )
@split( buildNrdb50, "blast.dir/*.fasta" )
def splitFasta( infiles, outfiles):
    '''split fasta file.'''
    
    infile = infiles[0]
    chunk_size = 500
    statement = '''
    cat %(infile)s
    | perl /ifs/devel/andreas/cgat/split_fasta.pl 
       -a blast.dir/chunk_%%s.fasta
       %(chunk_size)i 
    > split.log
    '''
    
    P.run()

@transform( splitFasta, suffix(".fasta"), ".blast.gz" )
def runBlast( infile, outfile ):
    '''run blast
    '''
    
    to_cluster = True
    max_evalue = 1.0
    num_results = 1000000
    mask = 21
    dbsize = 1500000000
    # upper case is critical, otherwise traceback fails!?
    matrix = "BLOSUM50"
    gop = 12
    gep = 2
    dbname = "/tmp/blast/nrdb50"

    statement = '''
    /ifs/devel/andreas/cgat/run.py --log=%(outfile)s.log
    'blastp
       -query %(infile)s
       -db %(dbname)s
       -evalue %(max_evalue)f
       -num_alignments %(num_results)i
       -num_descriptions %(num_results)i
       -db_soft_mask %(mask)i
       -matrix %(matrix)s
       -gapopen %(gop)i
       -gapextend %(gep)i
       -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq"
    | python /ifs/devel/andreas/cgat/blast2table.py 
        --alignment-format=blocks
    | gzip
    > %(outfile)s';
    checkpoint;
    echo "#//" | gzip >> %(outfile)s
    '''

    P.run()

@merge( "blast.dir/*.log", None )
def removeBlastUnfinished( infiles, outfile ):
    '''remove aborted blast runs.'''

    deleted = 0

    for infile in infiles:
        line = IOTools.getLastLine( infile )
        
        if not re.search( "job finished", line ):
            fn = infile[:-len(".log")]
            if os.path.exists( fn ):
                P.info("deleting %s" % fn )
                os.unlink( fn )
                deleted += 1

    P.info("deleted %i files" % deleted)

@merge( runBlast, "blast.check.tsv" )
def checkBlastRuns( infiles, outfile ):
    '''check if output files are complete.
    '''
    
    outf = IOTools.openFile( outfile, "w" )

    outf.write( "chunkid\tquery_first\tquery_last\tfound_first\tfound_last\tfound_total\tfound_results\thas_finished\tattempts\t%s\n" %\
                    "\t".join(Logfile.RuntimeInformation._fields))

    for infile in infiles:
        E.debug( "processing %s" % infile)
        chunkid = P.snip( os.path.basename( infile ), ".blast.gz" )
        logfile = infile + ".log"
        chunkfile = P.snip( infile, ".blast.gz" ) + ".fasta"

        with IOTools.openFile( infile ) as inf:
            l = inf.readline()
            ids = set()
            total_results = 0
            for l in inf:
                if l.startswith("#//"): continue
                ids.add( int(l.split("\t")[0] ) )
                total_results += 1
            found_first = min(ids)
            found_last = max(ids)
            found_total = len(ids)

        l = IOTools.getFirstLine( chunkfile )
        query_first = l[1:-1]
        l2 = IOTools.getLastLine( chunkfile, nlines = 2).split("\n")
        query_last = l2[0][1:]

        logresults = Logfile.parse( logfile )
        
        outf.write( "\t".join( map(str, (\
                        chunkid, query_first, query_last,
                        found_first, found_last,
                        found_total, total_results,
                        logresults[-1].has_finished,
                        len(logresults),
                        "\t".join( map(str, logresults[-1]) ) ) ) ) + "\n" )
        
    outf.close()

@merge( runBlast, "pairsdb_50x50.links.gz" )
def mergeBlast( infiles, outfile ):
    '''merge blast results into a single file.'''

    to_cluster = True

    files = [ (int(re.match( ".*chunk_(\d+).blast.gz", x).groups()[0]), x) for x in infiles ]
    files.sort()

    files = " ".join( [ x[1] for x in files ] )

    statement = '''zcat %(files)s | awk '$1 == "query_nid" { if(a){ next;} a=1; } {print}' | gzip > %(outfile)s'''
    P.run()

    files = [ (int(re.match( ".*chunk_(\d+).blast.gz.log", x).groups()[0]), x) for x in infiles ]
    files.sort()

    files = " ".join( [ x[1] for x in files ] )

    statement = '''cat %(files)s  >> %(outfile)s.log'''
    P.run()

@transform( mergeBlast, suffix(".links.gz"), add_inputs(  buildNrdb50 ),
            ".summary")
def checkBlastRun( infiles, outfile ):
    '''build summary stats on file.'''

    pairsdbfile, seqfile = infiles
    
    nids = set()
    with IOTools.openFile( seqfile ) as inf:
        for r in FastaIterator.iterate( inf ):
            nids.add( int(r.title) )

    with IOTools.openFile( pairsdbfile ) as inf:
        query_ids, sbjct_ids = set(), set()
        total_results, self_links = 0, 0
        for l in inf:
            l = inf.readline()
            if l.startswith("#//"): continue
            query_id, sbjct_id = l.split("\t")[:2]
            query_ids.add( int(query_id) )
            sbjct_ids.add( int(sbjct_id) )
            if query_id == sbjct_id: self_links += 1
            total_results += 1

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "category\tcounts\n")
    outf.write( "\t".join( map(str, ('nids', len(nids)))) + "\n" )
    outf.write( "\t".join( map(str, ('links', total_results))) + "\n" )
    outf.write( "\t".join( map(str, ('self', self_links))) + "\n" )
    outf.write( "\t".join( map(str, ('queries', len(query_ids)))) + "\n" )
    outf.write( "\t".join( map(str, ('sbjcts', len(sbjct_ids)))) + "\n" )
    outf.close()

    outf = IOTools.openFile( outfile + '.missing_queries.gz', 'w' )
    outf.write( 'nid\n' )
    outf.write( "\n".join( map(str, sorted( list( nids.difference( query_ids )) ) )) + "\n" )
    outf.close()

    outf = IOTools.openFile( outfile + '.missing_sbjcts.gz', 'w' )
    outf.write( 'nid\n' )
    outf.write( "\n".join( map(str, sorted( list( nids.difference( sbjct_ids )) ) )) + "\n" )
    outf.close()

if __name__ == "__main__":
    sys.exit(P.main())
