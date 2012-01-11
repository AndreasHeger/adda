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
adda_tranlate.py - translate ADDA nids to sequence identifiers
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Translate ADDA internal numerical identifiers into
sequence identifiers.

The default maps the ADDA`s internal identifier (:term:`nid`) to 
the sequence identifier (:term:`pid`). The command line option
``-i/--invert" inverts this mapping.

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
from Adda import AddaIO

def main():
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"])

    parser.add_option( "-n", "--nids", dest="filename_nids", type="string",
                       help="filename with nids[default=%default].")

    parser.add_option( "-c", "--column", dest="columns", type="int", action="append",
                       help="columns with nids to translate (1-based) [default=%default].")

    parser.add_option( "-d", "--is-domains", dest="is_domains", action="store_true",
                       help="translate domain ids [default=%default].")

    parser.add_option( "-i", "--invert", dest="invert", action="store_true",
                       help="invert mapping [default=%default].")

    parser.add_option( "-e", "--no-header", dest="no_header", action="store_true",
                       help="file has no header [default=%default].")

    parser.set_defaults( 
        filename_nids = "adda.nids",
        columns = [],
        is_domains = False,
        invert = False,
        noheader = False,
        )
    
    (options, args) = E.Start( parser )
    
    map_nid2pid = AddaIO.readMapPid2Nid( open(options.filename_nids, "r") )
    if options.invert:
        E.info( "inverting mapping" )
        map_nid2pid = dict( [ (int(x[1]),str(x[0])) for x in map_nid2pid.iteritems()] )

    if len(options.columns) == 0: options.columns = [1]
    columns = [x-1 for x in options.columns ]

    toTuple, toDomain = AddaIO.toTuple, AddaIO.toDomain
    first = not options.no_header
    is_domains = options.is_domains
    ninput, noutput, nskipped = 0, 0, 0
    for line in options.stdin:
        if line.startswith("#"):
            options.stdout.write(line)
            continue

        if first:
            options.stdout.write(line)
            first = False
            continue
        
        ninput += 1

        data = line[:-1].split("\t")
        for x in columns:
            if is_domains:
                try:
                    d = toTuple(data[x])
                except ValueError:
                    E.warn( "could not parse domain `%s`" % data[x])
                    nskipped += 1
                    break

                try:
                    data[x] = toDomain( (str(map_nid2pid[d[0]]),d[1],d[2]) )
                except (IndexError, KeyError):
                    E.warn( "could not map domain `%s`" % data[x])
                    nskipped += 1
                    break
            else:
                try:
                    data[x] = str(map_nid2pid[int(data[x])])
                except IndexError:
                    E.warn( "could not map nid `%s`" % data[x])
                    nskipped += 1
                    break
        else:
            options.stdout.write("%s\n" % "\t".join(data))
            noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))
    E.Stop()

if __name__ == "__main__":    
    sys.exit(main())



