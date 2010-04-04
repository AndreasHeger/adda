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

Evaluate domain partitions against a reference set.

The script collects all overlapping domains between the ``test set`` and 
``reference set``. For each overlap, it collects the number of overlapping 
bases, the coverage of the ``test`` domain and the coverage of the ``reference``
domain. 

It outputs the following files:

overlaps.histogram
   histogram of relative coverage between ``test`` and ``reference`` domains
   (intersection/union).

domain_coverage.histogram
   distribution of coverage of ``test`` domains by ``reference`` domains.

reference_coverage.histogram
   distribution of coverage of ``reference`` domains by ``test`` domains.

For a good covering, both the coverage of ``reference`` domains by ``test``
domains and the coverage of ``test`` with ``reference`` domains should be
as close to 100% as possible.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----
"""

import string, getopt, sys, os, re, time, math, glob, optparse

import Adda.Experiment as E
import Adda.Histogram as Histogram
import Adda.Stats as Stats

from Adda.Pairsdb import *
# from Table_nrdb90_masks import Table_nrdb90_masks
from Adda.Table_nrdb import Table_nrdb
from Adda.TableDomains import TableDomains

def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )

    parser.add_option("-D", "--database", dest="database", type="string",          
                      help="tablename to use [default=%default]."  )

    parser.add_option("-t", "--trees", dest="table_name_trees", type="string",          
                      help="tablename with trees [default=%default]."  )

    parser.add_option("-r", "--parts", dest="table_name_parts", type="string",          
                      help="tablename with trees [default=%default]."  )

    parser.add_option("-b", "--bench", dest="table_name_bench", type="string",          
                      help="domain table to be benchmarked (for example: nrdb40_domo_domains_nr) [default=%default]."  )

    parser.add_option("-f", "--reference", dest="table_name_reference", type="string",          
                      help="table of reference table (for example: nrdb40_scop_domains_nr) [default=%default]."  )

    parser.add_option( "--bin-size", dest="bin_size", type="int",          
                      help="bin size [default=%default]."  )

    parser.add_option( "-o", "--resolution", dest="resolution", type="float",          
                      help="resolution for scaling of domains [default=%default]."  )

    parser.add_option("-s", "--switch", dest="switch", action = "store_true",
                      help="switch between coverage of reference and size ratio if coverage is 1 [default=%default]."  )

    parser.add_option("-k", "--skip-repeats", dest="skip_repeats", action = "store_true",
                      help="[default=%default]."  )

    parser.add_option("-m", "--skip-tms", dest="skip_tms", action = "store_true",
                      help="discard domains which contain transmembrane regions [default=%default]."  )

    parser.add_option("-e", "--check-selection", dest="check_selection", action = "store_true",
                      help="[default=%default]."  )

    parser.add_option("-q", "--quality", dest="quality", action = "store_true",
                      help="take only sequences which are curated [default=%default]."  )

    parser.add_option( "--no-full-length", dest="no_full_length", action = "store_true",
                      help="[default=%default]."  )

    parser.add_option( "--only-full-length", dest="only_full_length", action = "store_true",
                      help="[default=%default]."  )

    parser.add_option( "--check-if-comparable", dest="check_if_comparable", action = "store_true",
                      help="perform comparable check according to Islam95 (default level 85%) [default=%default]."  )

    parser.add_option( "--subset", dest="subset", type = "string",
                       help = "use only a subset of nids [default=%default]" )

    parser.set_defaults( 
        database = "pairsdb",
        table_name_reference = None,
        table_name_trees = None,
        table_name_parts = None,
        table_name_bench = None,
        resolution = None,
        loglevel = 1,
        min_overlap = 1,
        switch = 0,
        combine_repeats = 1,
        skip_repeats = 0,
        skip_tms = 0,
        discard_full_length = 0,
        check_selection = 0,
        selection_threshold = 0.9,
        quality = None,
        no_full_length = None,
        only_full_length = None,
        ## a full length domain should cover at least 90% of a sequence
        min_length_ratio = 0.9,
        check_comparable = None,
        check_comparable_level = 0.85,
        bin_size = 1,
        subset = None )

    (options, args) = E.Start( parser, 
                               argv = argv, 
                               add_output_options = True )

    dbhandle = Pairsdb()
    dbhandle.Connect()
    dbhandle.UseDatabase( options.database )
    
    tbl_reference = TableDomains(dbhandle, "generic")
    tbl_reference.SetName(options.table_name_reference)
    
    # tbl_masks = Table_nrdb90_masks(dbhandle)
    tbl_nrdb = Table_nrdb(dbhandle)

    if options.table_name_trees:

        nids_statement = "SELECT DISTINCT s.nid FROM %s AS s, %s AS t %%s WHERE t.rep_nid = s.nid %%s" %\
                         (options.table_name_trees, options.table_name_reference)

        if options.quality:
            nids_statement = nids_statement % (", nrdb_quality AS q", "AND q.nid = s.rep_nid AND q.is_curated = 'T'")
        else:
            nids_statement = nids_statement % ("","")
            
        statement = """
        SELECT t.node, t.parent, t.level, t.xfrom, t.xto,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (GREATEST( t.xto, %%i) - LEAST( t.xfrom, %%i))) AS ovl,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (t.xto - t.xfrom)) AS cov_dom,
        ((LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.xto - t.xfrom) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.nid = %%i
        AND (LEAST(t.xto, %%i) - GREATEST(t.xfrom, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (options.table_name_trees)
        
    elif options.table_name_parts or options.table_name_bench:

        if options.table_name_parts:
            table_name = options.table_name_parts
        else:
            table_name = options.table_name_bench

        if options.subset:
            nids_statement = "SELECT DISTINCT s.nid FROM %s AS s, %s AS t WHERE t.rep_nid = s.nid" % (options.subset, table_name)
        else:
            nids_statement = "SELECT DISTINCT s.rep_nid FROM %s AS s, %s AS t %%s WHERE t.rep_nid = s.rep_nid %%s" %\
                             (table_name, options.table_name_reference)

            if options.quality:
                nids_statement = nids_statement % (", nrdb_quality AS q", "AND q.nid = s.rep_nid AND q.is_curated = 'T'")
            else:
                nids_statement = nids_statement % ("","")

        statement = """
        SELECT 1, 0, 0, t.rep_from, t.rep_to,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (GREATEST( t.rep_to, %%i) - LEAST( t.rep_from, %%i))) AS ovl,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (t.rep_to - t.rep_from)) AS cov_dom,
        ((LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.rep_to - t.rep_from) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.rep_nid = %%i
        AND (LEAST(t.rep_to, %%i) - GREATEST(t.rep_from, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (table_name)
    else:
        print "what shall I compare?"
        sys.exit(1)

    if options.check_selection:
        selection_statement = """
        SELECT t.domain_from, t.domain_to,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (GREATEST( t.domain_to, %%i) - LEAST( t.domain_from, %%i))) AS ovl,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (t.domain_to - t.domain_from)) AS cov_dom,
        ((LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i)) / (%%i - %%i)) AS cov_ref,
        ((t.domain_to - t.domain_from) / (%%i - %%i)) AS rat_ref
        FROM %s AS t
        WHERE t.domain_nid = %%i
        AND (LEAST(t.domain_to, %%i) - GREATEST(t.domain_from, %%i) > %%i)
        ORDER BY ovl DESC
        LIMIT 1
        """ % (options.table_name_parts)
        options.table_name_parts = None
        
        parts_same_as_trees, parts_larger_than_trees, parts_smaller_than_trees, parts_much_smaller_than_trees =  0,0,0,0

    
    nids = map(lambda x:x[0],dbhandle.Execute(nids_statement).fetchall())

    overlaps = []
    cov_doms = []
    cov_refs = []
    touched  = {}

    if options.check_selection:
        options.stdout.write( "NID\tDNODE\tDPARENT\tDLEVEL\tDFROM\tDTO\tRID\tRFROM\tRTO\tOVL\tDCOV\tRCOV\tRRCOV\tMRCOV\n" )
    else:
        options.stdout.write( "NID\tDNODE\tDPARENT\tDLEVEL\tDFROM\tDTO\tRID\tRFROM\tRTO\tOVL\tDCOV\tRCOV\tRRCOV\tMRCOV\n" )

    E.info( "--> processing %i nids" % len(nids) )

    nskipped_no_assignments = 0
    nskipped_no_overlap = 0
    nskipped_wrong_domaintype = 0
    nfound = 0
    
    it = 0
    for nid in nids:

        it += 1

        E.debug( "--> processing %i" % nid )

        domains = tbl_reference.GetDomainBoundariesForNid( nid )

        length = tbl_nrdb.GetLength( nid )
        
        if not domains:
            nskipped_no_assignments +=1
            continue

        if options.no_full_length and len(domains) == 1:
            ## check if domain is actually full length, otherwise keep
            id, domain_from, domain_to = domains[0]
            if float(domain_to-domain_from+1) / float(length) >= options.min_length_ratio:
                nskipped_wrong_domaintype += 1
                continue
            
        if options.only_full_length:
            if len(domains) == 1:
                id, domain_from, domain_to = domains[0]
                if float(domain_to-domain_from+1) / float(length) <= options.min_length_ratio:
                    nskipped_wrong_domaintype += 1
                    continue
            else:
                nskipped_wrong_domaintype += 1                
                continue

        nfound += 1
        
        last_id = None
        x = 0

        while x < len(domains):
            
            id, domain_from, domain_to = domains[x]
                
            is_repeat = -1
            
            while x < len(domains) and domains[x][0] == id:
                domain_to = domains[x][2]
                x += 1
                is_repeat += 1

            if options.skip_repeats and is_repeat:
                continue

            # if options.skip_tms and tbl_masks.HasMask( nid, 2, domain_from, domain_to):
            #    continue

            if options.resolution:
                xdomain_from = int(float(domain_from-1)/options.resolution)
                xdomain_to   = int(float(domain_to-1)/options.resolution) + 1
            else:
                xdomain_from = domain_from
                xdomain_to   = domain_to

            E.debug( "processing domain %s (%i-%i) (%i-%i)" % ( id, domain_from, domain_to, xdomain_from, xdomain_to))

            s = statement % (xdomain_to, xdomain_from, xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from, xdomain_to, xdomain_from,
                             xdomain_to, xdomain_from,
                             nid,
                             xdomain_to, xdomain_from, options.min_overlap)

            if options.loglevel >= 4: print s
                
            result = dbhandle.Execute(s).fetchone()

            if not result: continue

            node, parent, level, xfrom, xto, overlap, cov_dom, cov_ref, rat_ref = result

            key = "%i-%s-%i-%i" % (nid, id, xfrom, xto)
            if touched.has_key(key):
                continue
            else:
                touched[key] = 1

            # discard full length domains
            if options.discard_full_length:
                if options.table_name_trees:            
                    if node == 0: continue
                else:
                    if length == xto - xfrom + 1: continue
            
            if options.switch and cov_ref == 1.0:
                xcov_ref = rat_ref
            else:
                xcov_ref = cov_ref
                
            # check, if selection did take a domain lower or further up
            if options.check_selection:
                yfrom = (xfrom * 10) + 1
                yto   = min(xto * 10 + 1, length)
                s = selection_statement % (yto, yfrom, yto, yfrom,
                                           yto, yfrom,
                                           yto, yfrom, yto, yfrom,
                                           yto, yfrom,
                                           nid,
                                           yto, yfrom, options.min_overlap)
                
                result = dbhandle.Execute(s).fetchone()
                if result:
                    parts_from, parts_to, ovl_parts, cov_parts, cov_tree, rat_parts = result


                    if rat_parts > 1.0:
                        parts_larger_than_trees += 1
                        token = ">"
                    elif rat_parts == 1.0:
                        parts_same_as_trees += 1
                        token = "="
                    else:
                        parts_smaller_than_trees += 1
                        token = "<"
                        if rat_parts < options.selection_threshold:
                            parts_much_smaller_than_trees += 1

                    options.stdout.write(string.join(map(str, (nid,
                                                               id, domain_from, domain_to,
                                                               level,
                                                               yfrom, yto,
                                                               parts_from, parts_to,
                                                               overlap, cov_dom, cov_ref, rat_ref, xcov_ref,
                                                               ovl_parts, cov_parts, cov_tree, rat_parts,
                                                               token)), "\t") + "\n")
                    
            else:
                options.stdout.write(string.join(map(str, (nid, node, parent, level, xfrom, xto,
                                                           id,
                                                           xdomain_from, xdomain_to,
                                                           overlap, cov_dom, cov_ref, rat_ref, xcov_ref)), "\t") + "\n")
                
                overlaps.append( int(overlap * 100) )
                cov_doms.append( int(cov_dom * 100) )
                cov_refs.append( int(xcov_ref * 100) )            


    E.info( "skipped nids because of no overlap with reference: %i" % nskipped_no_overlap )
    E.info( "skipped nids because of no assignments: %i" % nskipped_no_assignments )
    E.info( "skipped nids because of wrong domain type: %i" % nskipped_wrong_domaintype)
    E.info( "nids in comparison: %i" % nfound)
        
    if options.check_selection:
        E.info( " parts larger than trees=", parts_larger_than_trees )
        E.info( " parts like trees=", parts_same_as_trees )
        E.info( " parts smaller than trees=", parts_smaller_than_trees )
        E.info( " parts much smaller than trees (<%f)=" % options.selection_threshold, parts_much_smaller_than_trees )
        
    else:
        outfile_stats = E.openOutputFile( "stats" )
        outfile_stats.write("section\t%s\n" % Stats.Summary().getHeader())
        outfile_stats.write("overlaps\t%s\n" % str( Stats.Summary( overlaps ) ) )
        outfile_stats.write("domain_coverage\t%s\n" % str( Stats.Summary( cov_doms ) ) )
        outfile_stats.write("reference_coverage\t%s\n" % str( Stats.Summary( cov_refs ) ) )
        outfile_stats.close()

        outfile = E.openOutputFile( "overlaps.histogram" )
        outfile.write( "bin\tcounts\n")
        Histogram.Write(outfile, 
                        Histogram.Calculate( overlaps, 
                                             min_value=0, 
                                             increment=1, 
                                             no_empty_bins = True))
        outfile.close()

        outfile = E.openOutputFile( "domain_coverage.histogram" )
        outfile.write( "bin\tcounts\tfreq\tcumul_counts\tcumul_freq\treverse_counts\treverse_freq\n" )
        Histogram.Write(outfile,
                        Histogram.AddRelativeAndCumulativeDistributions(
                        Histogram.Calculate( cov_doms, 
                                             min_value=0, 
                                             increment=options.bin_size, 
                                             no_empty_bins = True)))
        outfile.close()

        outfile = E.openOutputFile( "reference_coverage.histogram" )
        outfile.write( "bin\tcounts\tfreq\tcumul_counts\tcumul_freq\treverse_counts\treverse_freq\n" )
        Histogram.Write(outfile,
                        Histogram.AddRelativeAndCumulativeDistributions(
                    Histogram.Calculate( cov_refs, 
                                         min_value=0, 
                                         increment=options.bin_size, 
                                         no_empty_bins = True)))
                        
        outfile.close()
    
    E.Stop()

if __name__ == "__main__":                                                                                                                                                                                     
    sys.exit( main( sys.argv) )  



















