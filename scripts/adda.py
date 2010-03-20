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
:Release: $Id: pipeline_kamilah.py 2869 2010-03-03 10:20:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

The main controlling script for the ADDA domain clustering method.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----
"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil
import fileinput, collections, gzip

import Adda.Experiment as E
import Adda.Pipeline as P
from ruffus import *
import csv
import sqlite3
from Adda import IndexedFasta, FastaIterator, IOTools

PARAMS = P.getParameters( "adda.ini" )

ADDA_STATEMENT="adda.py %(cmd)s"

def sequences(infile, outfile ):
    cmd = "sequences"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["input_graph"], PARAMS["output_graph"])
def graph(infile, outfile):
    cmd = "graph"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_graph"], PARAMS["output_fit"])
def fit(infile, outfile ):
    cmd = "fit"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_graph"], PARAMS["output_segments"])
def segment(infile, outfile):
    cmd = "segment"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_graph"], PARAMS["output_stats"])
def stats(infile, outfile):
    cmd = "stats"
    statement = ADDA_STATEMENT
    P.run()

@files( (PARAMS["output_fit"], PARAMS["output_segments"]), 
        PARAMS["output_domains"])
def optimise(infile, outfile):
    cmd = "optimise"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_domains"], PARAMS["output_domaingraph"])
def convert(infile, outfile):
    cmd = "convert"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_domaingraph"], PARAMS["output_mst"])
def mst(infile, outfile):
    cmd = "mst"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_mst"], PARAMS["output_mst"] + ".components")
def mst_components(infile, outfile):
    cmd = "mst-components"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_mst"], PARAMS["output_align"])
def align(infile, outfile):
    cmd = "align"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_align"], PARAMS["output_cluster"])
def cluster(infile, outfile):
    cmd = "cluster"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_cluster"], PARAMS["output_families"])
def families(infile, outfile):
    cmd = "families"
    statement = ADDA_STATEMENT
    P.run()

@files( PARAMS["output_families"], PARAMS["output_summary"])
def summary(infile, outfile):
    cmd = "summary"
    statement = ADDA_STATEMENT
    P.run()


#########################################################################
#########################################################################
#########################################################################
@merge( (PARAMS["filename_target_sequences"], 
         PARAMS["filename_adda_sequences"]),
         "target.new.fasta" )
def collectTargetSequences( infiles, outfile ):
    '''extract new sequences from input.'''
        
    filename_target, filename_adda = infiles
    statement = '''
	python %(scriptsdir)s/map_fasta2fasta.py 
		--filename-reference=%(filename_adda)s
                --output-filename-pattern=target.%%s
		%(filename_target)s > %(outfile)s.log
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@files( PARAMS["filename_adda_sequences"], "query.fasta" )
def collectADDASequences( infile, outfile ):
    '''unpack adda sequences.'''

    if infile.endswith(".gz"):
        statement = '''gunzip < %(infile)s > %(outfile)s'''
    else:
        statement = '''ln -s %(infile)s %(outfile)s'''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@files( (collectADDASequences, collectTargetSequences), "5.ooc" )
def buildBlatIndex( infiles, outfile):
    '''build blat index.'''
    infiles = " ".join( infiles )

    statement = '''
    blat -dots=100 -prot 
                -makeOoc=%(outfile)s 
		-minIdentity=%(map_min_identity)i
		%(infiles)s %(outfile)s.log < /dev/null >> %(outfile)s.log
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@files( PARAMS["filename_target_sequences"], "target.lengths" )
def collectSequenceLengths( infile, outfile ):
    '''get sequence lengths from input file.'''

    inf = gzip.open(infile)
    outf = open( outfile, "w" )
    outf.write( "id\tlength\n" )
    
    for seq in FastaIterator.FastaIterator(inf):
        pid = re.sub("\s.*", "", seq.title )
        outf.write( "%s\t%i\n" %  (pid, len(seq.sequence) ))
        
    outf.close()
    inf.close()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("blat.dir") )
@split( collectTargetSequences, "blat.dir/chunk_*.fasta" )
def splitSequenceFile( infile, outfiles ):

    # patch ruffus bug
    if type(infile) == type(list()):
        infile = infile[0]

    statement = '''
       perl %(scriptsdir)s/split_fasta.pl 
            -a blat.dir/chunk_%%s.fasta %(map_chunksize)i
            < %(infile)s > split.log
       '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( splitSequenceFile, 
            suffix( ".fasta"), 
            inputs( [".fasta", collectADDASequences] ),
            ".blat.gz" )
def runBlat( infiles, outfile ):
    '''run a blat job.'''

    to_cluster = True
    infile, fasta = infiles
    statement = '''
    blat  
	  -prot
	  -ooc=5.ooc
	  -noHead
	  -minIdentity=%(map_min_identity)i 
	  %(fasta)s
	  %(infile)s
          stdout | gzip > %(outfile)s
    '''
    
    P.run()

@transform( runBlat, suffix( ".blat.gz"), ".domains" )
def mapDomains( infile, outfile ):
    '''collect blat matching stats.'''

    to_cluster= True
    job_options = "-l mem_free=4000M"
    statement = '''gunzip 
        < %(infile)s 
	| python %(scriptsdir)s/map_blat2adda.py 
		--filename-domains=<( gunzip < %(map_filename_domains)s)
		--output-filename-pattern="%(outfile)s.%%s" 
		--log=%(outfile)s.log 
		--verbose=2 
        > %(outfile)s
        '''
    P.run()

def mergeWithHeader( infiles, outfile ):
    first = True
    outf = open( outfile, "w" )
    for line in fileinput.input( infiles ):
        if line.startswith("id"): 
            if not first: continue
            first = False
        outf.write(line)
    outf.close()
    
@merge( mapDomains, "mapped.stats" )
def collectMappingStats( infiles, outfile ):
    '''collect mapping stats.'''
    
    for x in ("full", "good", "partial", "log"):
        outf = open( "%s.%s" % (outfile,x), "w")
        for line in fileinput.input( [ "%s.%s" % (x,y) for y in infiles ] ):
            outf.write( line )
        outf.close()

    for x in ("mapped", "aggregate"):
        outf = open( "%s.%s" % (outfile,x), "w")
        outf.write( "bin\tcounts\n" )
        for line in fileinput.input( [ "%s.%s" % (x,y) for y in infiles ] ):
            if line.startswith("#"): continue
            if line.startswith("bin"): continue
            d = line[:-1].split("\t")
            outf.write( "%s\t%i\n" % (d[0], sum( map(int, d[1:] ))))
        outf.close()

@merge( mapDomains, "mapped.domains" )
def buildMappedDomains( infiles, outfile ):
    '''collect domains mapped via BLAT.'''
    mergeWithHeader( infiles, outfile )

@merge( (collectTargetSequences, PARAMS["map_filename_domains"]), "direct.domains" )
def buildDirectDomains( infiles, outfile ):
    '''collect domains that could be transfered without mapping.'''
    
    x, filename_domains = infiles

    statement = '''gunzip
        < %(filename_domains)s 
	| python %(scriptsdir)s/substitute_tokens.py 
		--apply=target.new2old.map 
		--invert 
		--column=1 
		--filter 
	> %(outfile)s
    '''
    P.run()

@merge( (buildMappedDomains, buildDirectDomains), "indirect.domains" )
def buildIndirectDomains( infiles, outfile ):
    '''collect domains mapped from domains mapped via BLAT.'''
    
    infiles = " ".join(infiles)
    statement = '''
	cat %(infiles)s |
	python %(scriptsdir)s/substitute_tokens.py 
		--apply=target.new2new.map
		--column=1 
		--invert \
		--filter > %(outfile)s
    '''
    P.run()

@merge( (buildMappedDomains, buildDirectDomains, buildIndirectDomains), "overlap.table" )
def buildOverlapTable( infiles, outfile ):
    '''calculate overlap between the different sources of domains.'''
    infiles = " ".join(infiles)
    statement = '''
    python %(scriptsdir)s/set_diff.py --add-percent %(infiles)s > %(outfile)s
    '''
    P.run()

@merge( (buildMappedDomains, buildDirectDomains, buildIndirectDomains), "all.domains" )
def buildMappingAllDomains( infiles, outfile ):
    '''calculate overlap between the different sources of domains.'''
    mergeWithHeader( infiles, outfile )

@files( buildMappingAllDomains, "all.families" )
def buildMappingAllFamilies( infile, outfile ):
    '''get family counts from a domains file.'''
    

    if type(infile) == type(list()):
        infile = infile[0]

    counts = collections.defaultdict( int )

    for line in open(infile):
        if line.startswith("id"): continue
        if line.startswith("#"): continue
        data = line[:-1].split("\t")
        counts[data[3]] += 1

    outf = open(outfile,"w")
    outf.write("family\tcounts\n" )
    for family, c in counts.iteritems():
        outf.write( "\t".join( (family,str(c)))+"\n")
    outf.close()

@merge( (buildMappingAllDomains, collectSequenceLengths), "mapping.coverage" )
def buildMappingCoverage( infiles, outfile ):
    '''compute coverage of target sequences with ADDA domains.'''
    
    filename_domains, filename_lengths = infiles

    statement = '''
    python %(scriptsdir)s/adda2coverage.py 
		--log=%(outfile)s.log 
		--filename-lengths=%(filename_lengths)s 
                --output-filename-pattern="%(outfile)s_%%s"
    < %(filename_domains)s 
    > %(outfile)s
    '''
    P.run()

@merge( (PARAMS["map_filename_domains"], 
         PARAMS["map_filename_adda_sequences"],
         PARAMS["map_filename_target_sequences"],
         collectTargetSequences,
         buildMappedDomains,
         buildDirectDomains,
         buildIndirectDomains,
         buildMappingAllDomains,
         buildMappingAllFamilies ), "mapping.summary" )
def buildMappingSummary( infiles, outfile ):
    '''build summary with mapping stats
    '''
        
    (adda_domains, 
     adda_sequences, 
     target_sequences,
     new_sequences,
     mapped_domains,
     direct_domains,
     indirect_domains,
     all_domains,
     all_families) = infiles

    outf = open( outfile, "w")

    nadda_sequences = len( list(FastaIterator.FastaIterator( gzip.open( adda_sequences))))
    ntarget_sequences = len( list(FastaIterator.FastaIterator( gzip.open( target_sequences))))

    outf.write( "set\tdomains\tsequences\tfamilies\n")
    outf.write( "%s\tna\t%i\tna\n" % \
                    (adda_sequences, nadda_sequences ))
    outf.write( "%s\tna\t%i\tna\n" % \
                    (target_sequences, ntarget_sequences) )
    
    def _countDomains( infile ):
        Counts = collections.namedtuple('Counts', 'ndomains, nsequences, nfamilies')

        sequences, ndomains, families = set(), 0, set()
        if infile.endswith(".gz"):
            inf = gzip.open(infile)
        else: 
            inf = open(infile)

        for line in inf:
            if line.startswith("#"): continue
            if line.startswith("id"): continue
            if line.startswith("nid"): continue
            data = line[:-1].split("\t")
            sequences.add( data[0] )
            ndomains += 1
            families.add( data[3] )

        return Counts._make( (ndomains, len(sequences), len(families) ) )

    for x,f in (
            ("adda", adda_domains), 
            ("mapped", mapped_domains), 
            ("direct", direct_domains), 
            ("indirect", indirect_domains),
            ("all", all_domains)):
        c = _countDomains(f)
        outf.write( "%s\t%i\t%i\t%i\n" % (x, c.ndomains, c.nsequences, c.nfamilies))
        
    outf.close()

@follows(
    collectTargetSequences,
    collectADDASequences,
    buildBlatIndex,
    splitSequenceFile,
    runBlat,
    mapDomains,
    buildMappingAllDomains,
    buildMappingAllFamilies,
    buildMappingCoverage,
    buildMappingSummary,
    )
def map(): pass


@follows( summary )
def adda(): pass

if __name__ == "__main__":
    # P.checkExecutables( ("blat", "gunzip", ))
    sys.exit(P.main())




