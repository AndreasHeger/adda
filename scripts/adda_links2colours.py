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
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Build a colour map from an edge list graph.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----
"""

colors=( "red", "blue", "green", "yellow", "cyan", "magenta",
         "lightblue", "lightred", "lightgreen", "orange")

increment_colors= ("lightgrey", "darkgrey")        

import os, sys, string, re, getopt, optparse

import Adda.Experiment as E

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )

    parser.add_option("-m", "--multi-labels", dest="multi_labels", action="store_true",
                      help="if label is a comma-separated list, use grey for for colours [default=%default]." )

    parser.add_option("-l", "--legend", dest="legend", action="store_true",
                      help="add legend to color list [default=%default]." )

    parser.add_option( "-1", "--label1", dest="label1", type=int,
                       help = "column to use for first label [default=%default]" )

    parser.add_option( "-2", "--label2", dest="label2", type=int,
                       help = "column to use for first label [default=%default]" )

    parser.set_defaults( 
        multi_labels = None,
        legend = None,
        label1 = 3,
        label2 = 4,
        )

    (options, args) = E.Start( parser )
    
    current_color = 0

    labels = { "na" : "white" }

    take = ( 0, 1, options.label1 - 1, options.label2 -1 )

    options.stdout.write( "node\tcolour\tlabel\n" )

    for line in sys.stdin:

        if line.startswith("#"): continue
        
        data = line[:-1].split("\t")
        node1, node2, label1, label2 = [ data[x] for x in take ]

        for node, label in ( (node1,  label1), 
                             (node2, label2 ) ):

            color = None

            if node in labels: continue
            
            n = label.count(",")

            if not labels.has_key(label):            

                if n and options.multi_labels:
                    n = n % len(increment_colors)
                    labels[label] = increment_colors[n]
                else:
                    labels[label] = colors[current_color]
                    current_color = (current_color + 1) % len(colors)
                    
            options.stdout.write("%s\t%s\t%s\n" % (node, labels[label], label))

    if options.legend:
        options.stdout.write("# legend:")
        for x in labels.keys():
            if x == "18": continue
            options.stdout.write( "%s=%s;" % (x, labels[x]) )
        options.stdout.write( "\n" )

    E.Stop()
