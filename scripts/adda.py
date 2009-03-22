#!/usr/bin/env python

USAGE="""adda.py [OPTIONS] cmds

interface to compute adda
"""

import sys, os, re, time, math, copy, glob, optparse
from multiprocessing import Process
import fileinput
# segfault on my machine without the print statement - strange
print 
import Adda.Experiment as E
from Adda import *

def run( options, order, map_module, config, fasta = None ):

    if "all" in options.steps:
        steps = order
    else:
        steps = options.steps

    modules = []

    for step in order:
        if step in steps:
            module = map_module[step]( options, config, fasta = fasta )
            
            if not module.isComplete():
                modules.append( module )
            else:
                E.info( "%s complete" % step )
    print steps, order, modules        
    if len(modules) == 0: return

    E.info( "running apply on modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.run()
        module.finish()

def mergeGraph( options, 
                order, 
                map_module, 
                config, 
                fasta ):


    if "all" in options.steps:
        steps = order
    else:
        steps = options.steps

    modules = []
    for step in order:
        if step in steps:
            modules.append( map_module[step]( options, 
                                              config, 
                                              fasta = fasta ) )

    E.info( "performing merge on modules: %s" % str(modules) )

    for module in modules:
        module.merge()

class Runner( Process ):
    """run adda jobs."""

    def __init__(self, options, order, map_module, config, chunk, nchunks ):

        Process.__init__( self )

        fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

        if "all" in options.steps: steps = order
        else: steps = options.steps

        modules = []
        for step in order:
            if step in steps:
                module = map_module[step]( options, 
                                           config = config, 
                                           fasta = fasta,
                                           num_chunks = nchunks,
                                           chunk = chunk )
            
                module.startUp()
                
                if not module.isComplete():
                    modules.append( module )
                else:
                    E.info( "%s complete" % step )
        
        self.mModules = modules
        
class RunnerOnFile(Runner):
    """run jobs on a file per line."""

    def __init__(self, filename, options, order, map_module, config, chunk, nchunks ):

        Runner.__init__( self, options, order, map_module, config, chunk, nchunks )

        E.info( "opening file %s at chunk %i" % (filename, chunk) )
        
        self.mIterator = FileSlice.Iterator( filename, 
                                             nchunks,
                                             chunk,
                                             FileSlice.iterator )

    def run( self ):

        if len(self.mModules) == 0: return

        E.info( "running apply on modules: %s" % (",".join(map(str, self.mModules))) )

        for line in self.mIterator:
            for module in self.mModules:
                module.apply( line )

        for module in self.mModules:
            module.finish()
    
class RunnerOnGraph(Runner):
    """run jobs in parallel on graph."""

    def __init__(self, filename, options, order, map_module, config, chunk, nchunks ):

        Runner.__init__( self, options, order, map_module, config, chunk, nchunks )

        if len(self.mModules) == 0: return

        self.mOldAlignmentFormat = options.old_alignment_format

        E.info( "opening graph %s at chunk %i" % (filename, chunk) )

        self.mIterator = FileSlice.IteratorMultiline( filename, 
                                                      nchunks,
                                                      chunk,
                                                      FileSlice.groupby,
                                                      key = lambda x: x[:x.index("\t")] )

        self.mMapId2Nid = AddaIO.readMapId2Nid( open(config.get( "files", "output_nids", "adda.nids" ), "r" ) )

    def run( self ):

        if len(self.mModules) == 0: return

        E.info( "running apply on modules: %s" % (",".join(map(str, self.mModules))) )

        for record in self.mIterator:
            neighbours = []
            for line in record:
                n = AddaIO.NeighbourRecord( line )
                if self.mMapId2Nid:
                    if (n.mQueryToken not in self.mMapId2Nid or \
                            n.mSbjctToken not in self.mMapId2Nid ):
                        continue 
                    q = n.mQueryToken = self.mMapId2Nid[n.mQueryToken]
                    n.mSbjctToken = self.mMapId2Nid[n.mSbjctToken]
  
                if self.mOldAlignmentFormat:
                    n.mQueryFrom -= 1
                    m.mSbjctFrom -= 1
  
                neighbours.append( n )

            if neighbours:
                for module in self.mModules:
                    module.apply( AddaIO.NeighboursRecord( q, neighbours ) )

        E.info( "running finish on modules: %s" % (",".join(map(str, self.mModules))) )
        
        for module in self.mModules:
            module.finish()

def mergeConvert( options ):
    pass

def runParallel( runner, filename, options, order, map_module, config ):

    nchunks = config.get( "adda", "num_jobs", 4 )

    E.info( "running %i jobs" % nchunks )

    processes = []
    for chunk in range(nchunks):
        E.info( "starting job %i" % chunk )
        process = runner( filename, options, order, map_module, config, chunk, nchunks )
        process.start()
        processes.append( process )

    for p in processes:
        p.join()

    E.info( "all jobs finished" )

def main():
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "--config", dest="filename_config", type="string",
                      help="configuration file [default=%default].")

    parser.add_option( "--force", dest="force", action="store_true",
                      help="overwrite existing files [default=%default].")

    parser.add_option( "--continue", dest="append", action="store_true",
                      help="continue from an aborted run and append to existing files [default=%default].")

    parser.add_option( "--test", dest="test", type="int",
                      help="run a test with first # sequences [default=%default]")
    
    parser.add_option( "--old-alignment-format", dest="old_alignment_format", action="store_true",
                       help = "input graph is in old 1-based coordinates." )

    parser.add_option( "--steps", dest="steps", type="choice", action="append",
                       choices=("all", 
                                "sequences",
                                "blast",
                                "fit", 
                                "graph",
                                "index",
                                "check-index",
                                "profiles",
                                "segment", 
                                "merge-graph",
                                "optimise",
                                "convert",
                                "merge-convert",
                                "mst", 
                                "align",
                                "merge-align",
                                "cluster", ),
                       help="perform this step [default=%default]" )

    parser.add_option( "--start-at", dest="start_at", type="string",
                      help="start at sequence [default=%default]")

    parser.add_option( "--stop-at", dest="stop_at", type="string",
                      help="stop at sequenec [default=%default]")


    parser.set_defaults( 
                        filename_config = "adda.ini",
                        steps = [],
                        start_at = None,
                        stop_at = None,
                        old_alignment_format = False,
                        force = False,
                        append = False,
                        test = None,
                        temporary_directory = ".",
                        )
    
    (options, args) = E.Start( parser )

    config = AddaIO.ConfigParser()
    config.read( os.path.expanduser( options.filename_config ) )

    if args: options.steps = args
        
    ## collect modules and initialise them         
    map_module = { 'fit' : AddaFit.AddaFit,
                   'segment' : AddaSegment.AddaSegment,
                   'blast' : AddaBlast.AddaBlast,
                   'graph' : AddaGraph.AddaGraph,
                   'profiles' : AddaProfiles.AddaProfiles, 
                   'index' : AddaIndex.AddaIndexBuild,
                   'check-index' : AddaIndex.AddaIndexCheck,
                   'optimise' : AddaOptimise.AddaOptimise,  
                   'sequences' : AddaSequences.AddaSequences,
                   'convert' : AddaConvert.AddaConvert,
                   'mst' : AddaMst.AddaMst, 
                   'align' : AddaAlign.AddaAlign, 
                   'cluster' : AddaCluster.AddaCluster,
                   }

    # modules and their hierarchy
#    run( options, 
#          order = ( "sequences", ), 
#          map_module = map_module,
#          config = config )

    runParallel( 
        RunnerOnGraph,
        filename = config.get( "files", "input_graph", "adda.graph" ),
        options = options, 
        order = ("fit", "segment", "graph", "profiles" ),
        map_module = map_module,
        config = config )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )
    
    mergeGraph( options,
                order = ("fit", "segment", "graph", "profiles" ),
                map_module = map_module,
                config = config,
                fasta = fasta )

    run( options, 
         order = ("index", "check-index", "optimise" ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    E.info( "building domain graph" )

    run( options, 
         order = ( "convert", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    mergeConvert( options )

    E.info( "computing minimum spanning tree" )
    run( options, 
         order = ( "mst", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    E.info( "alignment of domains" )

    runParallel( 
        RunnerOnFile,
        config.get( "files", "output_mst", "adda.mst" ),
        options = options, 
        order = ( "align", ),
        map_module = map_module,
        config = config )

    mergeAlign( options )

    run( options, 
         order = ( "cluster", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    E.Stop()
    
                
if __name__ == "__main__":    
    sys.exit(main())









