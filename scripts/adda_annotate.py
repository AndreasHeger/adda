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
from Adda import IndexedFasta, FastaIterator, IOTools, AddaIO

PARAMS = P.getParameters( "adda.ini" )

@files( PARAMS["eval_alignment_graph"], "alignment_graph.gz" )
def annotateAlignmentGraph( infile, outfile ):

    # collect benchmark domains 
    E.info( "reading benchmark domains" )
    benchmark_domains = AddaIO.readMapNid2Domains( 
        gzip.open( PARAMS["eval_benchmark_domains"] ) )

    totuple = AddaIO.toTuple
    # build map of id to nid
    E.info( "reading map between pid and nid" )
    map_nid2pid = AddaIO.readMapPid2Nid( open(PARAMS["eval_adda_nids"], "r") )

    def getOverlappingDomains( pid, start, end ):
        '''get domains overlapping pid:start..end'''
        if pid not in benchmark_domains: return ()
        # greedy overlap testing
        r = []
        for family, domains in benchmark_domains[pid].iteritems():
            for other_start, other_end in domains:
                if start >= other_end or end <= other_start: continue
                r.append( (family, other_start, other_end) )
        return r

    counts = E.Counter()
    
    inf = gzip.open( infile )
    outf = gzip.open( outfile, "w" )
    
    outf.write( "%s\n" % "\t".join( ( "passed",
                                      "qdomain",
                                      "sdomain",
                                      "weight",
                                      "qstart",
                                      "qend",
                                      "qali",
                                      "sstart",
                                      "send",
                                      "sali",
                                      "score",
                                      "naligned",
                                      "ngaps",
                                      "zscore",
                                      "rfamilies",
                                      "sfamilies",
                                      "rdomains",
                                      "sdomains")) )
                
    for link in AddaIO.iterate_tested_links( inf ):
        qnid, qstart, qend = totuple(link.qdomain)
        snid, sstart, send = totuple(link.sdomain)
        qpid = map_nid2pid[qnid]
        spid = map_nid2pid[snid]
        qfamily = sorted(getOverlappingDomains( qpid, qstart, qend ))
        sfamily = sorted(getOverlappingDomains( spid, sstart, send ))

        if link.passed == "+": counts.passed += 1
        else: counts.failed += 1

        outf.write( "%s\t%s\t%s\t%s\t%s\n" % \
                        ("\t".join( map(str,link) ), 
                         ",".join( sorted(set([x[0] for x in qfamily])) ),
                         ",".join( sorted(set([x[0] for x in sfamily])) ),
                         ",".join("%s_%i_%i" % x for x in qfamily ),
                         ",".join("%s_%i_%i" % x for x in sfamily )))
    inf.close()

    E.info( "%s" % str( counts ) )

@follows( annotateAlignmentGraph )
def benchmark(): pass

if __name__ == "__main__":
    # P.checkExecutables( ("blat", "gunzip", ))
    sys.exit(P.main())
