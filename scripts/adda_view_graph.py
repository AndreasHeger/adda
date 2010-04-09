#!/usr/bin/env python

USAGE="""view_graph.py [OPTIONS] cmds

view links in indexed graph.
"""

import sys, os, re, time, optparse, logging

import cadda
import Adda.Experiment as E
from Adda import IndexedFasta, AddaProfiles, AddaIO, AddaAlign
import alignlib

def main():
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "--method", dest="method", type="choice",
                       choices=("view", "align" ),
                       help="method to perform [default=%default].")

    parser.set_defaults( 
        filename_graph = "adda.graph",
        filename_index = "adda.graph.idx",
        method = "view",
        filename_fasta= "adda",
        filename_config = "adda.ini",
        append = False,
        force = False,
        )
    
    (options, args) = E.Start( parser )
    
    config = AddaIO.ConfigParser()
    config.read( os.path.expanduser( options.filename_config ) )

    index = cadda.IndexedNeighbours( options.filename_graph,
                                     options.filename_index )

    if options.method == "view":
        for nid in  args:
            nid = int(args[0])
    
            neighbours = index.getNeighbours( nid )
    
            for n in neighbours:
                print str(n)

    elif options.method == "align":
        
        fasta = IndexedFasta.IndexedFasta( options.filename_fasta )
        nid1,start1, end1 = AddaIO.toTuple( args[0] )
        nid2,start2, end2 = AddaIO.toTuple( args[1] )

        align = AddaProfiles.AddaProfiles( options, config, fasta = fasta )

        logoddor    = alignlib.makeLogOddorDirichlet( 0.3 )
        regularizor = alignlib.makeRegularizorDirichletPrecomputed()
        weightor    = alignlib.makeWeightor()

        alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, 
                                                  -10.0,
                                                  -1.0 )

        toolkit = alignlib.makeToolkit()
        alignlib.setDefaultToolkit( toolkit )
        toolkit.setWeightor( weightor )
        toolkit.setLogOddor( logoddor )
        toolkit.setRegularizor( regularizor )

        def _buildProfile( nid, start, end ):
            neighbours = index.getNeighbours( nid )
            mali = align.buildMali( nid, neighbours )
            prof = alignlib.makeProfile( mali )
            E.info( "nid: %i, neighours=%i" % (nid, len(neighbours)))
            prof.useSegment( start, end )
            prof.prepare()
            seq = fasta.getSequence( nid )
            return alignlib.makeSequence( seq ), prof

        seq1, prof1 = _buildProfile( nid1, start1, end1 )
        seq2, prof2 = _buildProfile( nid2, start2, end2 )

        result = alignlib.makeAlignmentVector()
        
        alignator.align( result, prof1, prof2 )
        
        options.stdout.write( "%s vs %s: score=%5.2f, length=%i, numgaps=%i, row_from=%i, row_to=%i, col_from=%i, col_to=%i\n" %\
                                  (nid1, nid2,
                                   result.getScore(),
                                   result.getLength(),
                                   result.getNumGaps(),
                                   result.getRowFrom(), result.getRowTo(),
                                   result.getColFrom(), result.getColTo()))

        f = alignlib.AlignmentFormatExplicit( result, seq1, seq2 )

        options.stdout.write( "%s\n" % str(f))

    E.Stop()

if __name__ == "__main__":    
    sys.exit(main())
