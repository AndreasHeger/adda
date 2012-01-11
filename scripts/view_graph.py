#!/usr/bin/env python

USAGE="""view_graph.py [OPTIONS] cmds

view links in indexed graph
"""

import sys, os, re, time, math, copy, glob, optparse, logging, traceback, shelve

import networkx as nx
import Adda.Experiment as E
import Adda.AddaIO
import cadda

def main():
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )


    parser.set_defaults( 
        filename_graph = "adda.graph",
        filename_index = "adda.graph.idx",
        )
    
    (options, args) = E.Start( parser )
    
    nid = int(args[0])

    index = cadda.IndexedNeighbours( options.filename_graph,
                                     options.filename_index )
            
    neighbours = index.getNeighbours( nid )
    
    for n in neighbours:
        print str(n)

    E.Stop()

if __name__ == "__main__":    
    sys.exit(main())
