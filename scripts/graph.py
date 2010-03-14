#!/usr/bin/env python

USAGE="""adda.py [OPTIONS] cmds

interface to compute adda
"""

import sys, os, re, time, math, copy, glob, optparse, logging, traceback, shelve

import networkx as nx
import Adda.Experiment as E
import Adda.AddaIO

def readAlignmentGraph( infile ):

    G = nx.Graph()
    for line in infile:
        if line.startswith("#"): continue
        if line.startswith("passed"): continue
        data = line[:-1].split("\t")
        try:
            (passed, start, end, 
             evalue,
             qstart, qend, qali, 
             sstart, send, sali,
             score, naligned, ngaps,
             zscore) = data
        except ValueError:
            print "parsing error in line `%s`" % data
            raise
        G.add_edge( start, end, (passed, float(score), int(naligned)) )
        
    return G

def main():
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "-f", "--format", dest="graph-format", type="choice",
                       choices=("alignments",),
                       help="graph format [default=%default].")

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       choices=("shortest-path", "translate", "components" ),
                       help="methods to apply [default=%default].")

    parser.add_option( "-a", "--filename-map", dest="filename_map", type="string",
                       help="filename mapping ids to nids (used for translation) [default=%default].")

    parser.add_option( "-1", "--node1", dest="node1", type="string",
                       help="first node for path calculation [default=%default].")

    parser.add_option( "-2", "--node2", dest="node2", type="string",
                       help="second node for path calculation [default=%default].")


    parser.set_defaults( 
        method = None,
        graph_format = "alignments",
        filename_map = None,
        node1 = None,
        node2 = None,
        )
    
    (options, args) = E.Start( parser )
    
            
    if options.method == "translate":
        
        if options.filename_map:
            E.info("reading map from %s" % options.filename_map)
            map_id2nid = Adda.AddaIO.readMapId2Nid( open( options.filename_map, "r") )
            map_nid2id = dict([[v,k] for k,v in map_id2nid.iteritems()])

        def translate_alignments( line ):        
            if line.startswith("passed"): return line
            data = line.split( "\t" )
            
            x = data[1].split("_")
            y = data[2].split("_")
            try:
                data[1] = "%s_%s_%s" % (map_nid2id[int(x[0])],x[1],x[2])
            except KeyError:
                sys.stderr.write("could not map: %s\n" % str(x) )
                raise
            try:
                data[2] = "%s_%s_%s" % (map_nid2id[int(y[0])],y[1],y[2])
            except KeyError:
                sys.stderr.write("could not map: %s\n" % str(y) )
                raise

            return "\t".join(data)

        if options.graph_format == "alignments":
            translator = translate_alignments
            
        for line in options.stdin:
            if not line.startswith("#"): 
                line = translator( line )
            options.stdout.write(line)
            
        E.Stop()
        return

    t = time.time()
    if options.graph_format == "alignments":
        G = readAlignmentGraph( options.stdin )
        
    E.info( "graph read in %i seconds" % (time.time() - t ))
    t = time.time()

    if options.method == "shortest-path":
        path = nx.shortest_path(G, options.node1, options.node2)
        last_node = path[0]
        for node in path[1:]:
            print "%s\t%s\t%s" % (last_node, node, str(G[last_node][node]))
            last_node = node

    elif options.method == "components":
        print "component\tnode"
        for id, component in enumerate(nx.connected_components( G )):
            for c in component:
                print "%i\t%s" % (id,c)
                
    E.info( "%s: %i seconds" % (options.method, time.time() - t ))
    E.Stop()

if __name__ == "__main__":    
    sys.exit(main())