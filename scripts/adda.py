#!/usr/bin/env python

USAGE="""adda.py [OPTIONS] cmds

interface to compute adda
"""

import sys, os, re, time, math, copy, glob, optparse, logging, traceback, shelve
from multiprocessing import Process, cpu_count, Pool
import multiprocessing

import fileinput
import cadda

import Adda.Experiment as E
from Adda import *

from logging import warn, info, debug

L = {}

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
                L.info( "%s complete" % step )

    if len(modules) == 0: return

    L.info( "working with modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.startUp()
        module.run()
        module.finish()

def merge( options, 
           order, 
           map_module, 
           config, 
           fasta ):
    """return True if all merging operations succeeded."""

    if "all" in options.steps: steps = order
    else: steps = options.steps

    nchunks = config.get( "adda", "num_slices", 10 )

    modules = []
    for step in order:
        if step in steps:
            modules.append( map_module[step]( options, 
                                              config, 
                                              fasta = fasta,
                                              chunk = 0,
                                              num_chunks = nchunks ) )

    L.info( "performing merge on modules: %s" % str(modules) )

    for module in modules:
        if not module.merge():
            return False
    return True


class Run(object):
    pass

class RunOnGraph(Run):
    
    def __init__(self, config, steps ):

        #from guppy import hpy
        #h = hpy()
        # ignore memory usage of previous object
        #h.setrelheap()

        L.info( "loading fasta sequences from %s" % config.get( "files", "output_fasta", "adda" ) ) 

        self.mFasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

        #print h.heap()

        if "all" in steps or "fit" in steps:

            L.info( "loading map_id2nid from %s" % config.get( "files", "output_nids", "adda.nids" ))
            infile = open( config.get( "files", "output_nids", "adda.nids" ) )
            self.mMapId2Nid = AddaIO.readMapId2Nid( infile, 
                                                    storage = config.get( "files", "storage_nids", "memory" ) )
            infile.close()

            L.info( "loading domain boundaries from %s" % config.get( "files", "input_reference") )
            infile = AddaIO.openStream( config.get( "files", "input_reference") )
            rx_include = config.get( "fit", "family_include", "") 

            self.mMapNid2Domains = AddaIO.readMapNid2Domains( infile, 
                                                              self.mMapId2Nid, 
                                                              rx_include,
                                                              storage = config.get( "files", "storage_domains", "memory" ) )
            infile.close()
            self.mMapId2Nid = None
        else:
            self.mMapNid2Domains = None
            self.mMapId2Nid = None
        #print h.heap()

        self.mFilenameGraph = config.get( "files", "output_graph", "adda.graph")
        self.mFilenameIndex = config.get( "files", "output_index", "adda.graph.index")

    def __call__(self, argv ):
        """run job, catching all exceptions and returning a tuple."""
        
        try:
            self.apply( argv )
            return None
        except:
            exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
            exception_stack  = traceback.format_exc(exceptionTraceback)
            exception_name   = exceptionType.__module__ + '.' + exceptionType.__name__
            exception_value  = str(exceptionValue)
            return (exception_name, exception_value, exception_stack)

    def apply( self, argv ):

        (chunk, nchunks, options, order, map_module, config ) = argv

        L.info( "chunk %i: setting up" % (chunk ))

        if "all" in options.steps: steps = order
        else: steps = options.steps

        # load all maps that were not inherited from the parent process
        if "all" in steps or "fit" in steps and self.mMapNid2Domains == None:
            L.info( "opening map_nid2domains from cache" )
            self.mMapNid2Domains = shelve.open( config.get( "files", "storage_domains", "memory" ), "r")

        # build the modules
        modules = []
        for step in order:
            if step in steps:
                if map_module[step]( options, 
                                     config = config, 
                                     fasta = self.mFasta ).isComplete():
                    L.info( "chunk %i: step %s is complete" % (chunk, step ))
                    continue

                module = map_module[step]( options, 
                                           config = config, 
                                           fasta = self.mFasta,
                                           num_chunks = nchunks,
                                           chunk = chunk,
                                           map_id2nid = self.mMapId2Nid,
                                           map_nid2domains = self.mMapNid2Domains )

                if not module.isComplete():
                    modules.append( module )
                else:
                    L.info( "chunk %i: step %s is complete" % (chunk,step) )

        if len(modules) == 0: 
            L.info( "chunk %i: nothing to be done" % chunk )
            return

        for module in modules: module.startUp()

        L.info( "chunk %i: starting work on modules: %s" % (chunk, ",".join(map(str, modules))) )

        # find out nids to work with
        nids = map(int, self.mFasta.keys())
        nids.sort()
        increment = int( math.ceil( len(nids) / float(nchunks) ) )
        start = chunk * increment
        nids = nids[start:start+increment]
        
        L.info( "chunk %i: starting work on %i nids from %s to %s" % (chunk, len(nids), str(nids[0]), str(nids[-1]) ) )

        index = cadda.IndexedNeighbours( self.mFilenameGraph, self.mFilenameIndex )

        iteration = 0
        for nid in nids:
            iteration += 1
            neighbours = index.getNeighbours( nid )

            L.info( "chunk %i: started nid=%s, neighbours=%i, progress=%i/%i" % (chunk, str(nid), len(neighbours), iteration, len(nids) ) )

            if neighbours:
                for module in modules:
                    module.run( AddaIO.NeighboursRecord( nid, neighbours ) )

            L.info( "chunk %i: finished nid=%s, neighbours=%i, progress=%i/%i" % (chunk, str(nid), len(neighbours), iteration, len(nids) ) )

        L.info( "chunk %i: running finish on modules: %s" % (chunk, ",".join(map(str, modules))) )

        for module in modules:
            module.finish()

        L.info( "chunk %i: finished  %i nids" % (chunk, len(nids)) )

def run_on_file( argv ):

    (filename, options, order, map_module, config, chunk, nchunks ) = argv

    L.info( "starting chunk %i on %s" % (chunk, filename) )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

    if "all" in options.steps: steps = order
    else: steps = options.steps

    modules = []
    for step in order:
        if step in steps:
            if map_module[step]( options, 
                                 config = config, 
                                 fasta = fasta ).isComplete():
                L.info( "%s complete" % step )
                continue

            module = map_module[step]( options, 
                                       config = config, 
                                       fasta = fasta,
                                       num_chunks = nchunks,
                                       chunk = chunk )

            if not module.isComplete():
                modules.append( module )
            else:
                L.info( "%s complete" % step )

    if len(modules) == 0: 
        L.info( "nothing to be done" )
        return

    L.info( "opening file %s at chunk %i" % (filename, chunk) )

    iterator = FileSlice.Iterator( filename, 
                                   nchunks,
                                   chunk,
                                   FileSlice.iterator )

    for module in modules: module.startUp()

    L.info( "working with modules: %s" % (",".join(map(str, modules))) )

    for line in iterator:
        if line.startswith("#"): continue

        for module in modules:
            module.run( line )

    for module in modules:
        module.finish()

    L.info( "running finish on modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.finish()

    L.info( "finished chunk %i on %s" % (chunk, filename) )

def getChunks( options, config ):
    """find out which chunks to compute from command line options."""
    nchunks = config.get( "adda", "num_slices", 10 )

    # set up the arguments for chunks to run
    if options.chunks == "all":
        chunks = range(nchunks) 
    else:
        ranges = options.chunks.split(",")
        chunks = []
        for r in ranges:
            s = r.split("-")
            if len(s) == 1: chunks.append( int(s[0]) )
            elif len(s) == 2: chunks.extend( list( range(int(s[0]), int(s[1]) ) ) )
            else: raise ValueError("can not parse range `%s`" % r )

        chunks = sorted(list(set(chunks)))
        if chunks[-1] >= nchunks: raise ValueError( "chunk `%i` out of range, maximum is " % (chunks[-1], nchunks-1 ) )

    return nchunks, chunks
        
def runParallel( runner, options, order, map_module, config ):
    """process filename in paralell."""

    if options.num_jobs:
        njobs = options.num_jobs 
    else:
        njobs = cpu_count()
    
    nchunks, chunks = getChunks( options, config )

    L.info( "running %i chunks in %i parallel jobs" % (len(chunks), njobs ))
    
    args = [ (chunk, nchunks, options, order, map_module, config ) for chunk in chunks ]

    logging.info('starting parallel jobs')

    pool = Pool( njobs )

    errors = pool.map( runner, args )
    pool.close()
    pool.join()

    errors = [ e for e in errors if e ]

    if errors:
        print "adda caught %i exceptions" % (len(errors))
        print "## start of exceptions"
        for exception_name, exception_value, exception_stack in errors:
            print exception_stack,
        print "## end of exceptions"
        sys.exit(1)

    L.info( "all jobs finished" )

def runSequentially( runner, options, order, map_module, config ):
    """process filename sequentially."""

    nchunks, chunks = getChunks( options, config )
    
    L.info( "running %i chunks sequentially" % (len(chunks) ))

    args = [ (chunk, nchunks, options, order, map_module, config ) for chunk in chunks ]

    for (job, argv) in enumerate(args):
        L.info( "job %i started" % job )
        error = runner( argv )

        if error:
            print "adda caught an exceptions"
            exception_name, exception_value, exception_stack = error
            print exception_stack,
            print "## end of exceptions"
            sys.exit(1)

        L.info( "job %i finished" % job )


    L.info( "all jobs finished" )

def oldrunSequentially( runner, filename, options, order, map_module, config ):
    """process filename sequentially."""

    nchunks, chunks = getChunks( options, config )
    
    L.info( "running %i chunks sequentially" % (len(chunks) ))

    args = [ (filename, options, order, map_module, config, chunk, nchunks ) for chunk in chunks ]

    for (chunk, argv) in enumerate(args):
        L.info( "job %i started" % chunk )
        error = runner( argv )

        if error:
            print "adda caught an exceptions"
            exception_name, exception_value, exception_stack = error
            print exception_stack,
            print "## end of exceptions"
            sys.exit(1)

        L.info( "job %i finished" % chunk )


    L.info( "all jobs finished" )
    
def main():
    global L
    
    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE )

    parser.add_option( "--config", dest="filename_config", type="string",
                      help="configuration file [default=%default].")

    parser.add_option( "--force", dest="force", action="store_true",
                      help="overwrite existing files [default=%default].")

    parser.add_option( "--continue", dest="append", action="store_true",
                      help="continue from an aborted run and append to existing files [default=%default].")

    parser.add_option( "--test", dest="test", type="int",
                      help="run a test with first # sequences [default=%default]")

    parser.add_option( "--num-jobs", dest="num_jobs", type="int",
                      help="use # processes. If not set, the number of CPUs/cores is taken [default=%default]")
    
    parser.add_option( "--alignment-format", dest="alignment_format", type="choice",
                       choices=("pairsdb", "pairsdb-old", "simap", "pairsdb-realign"),
                       help = "input format of graph. pairsdb-old: input graph is in old 1-based coordinates [default=%default]." )

    parser.add_option( "--chunks", dest="chunks", type="string",
                       help = "work on one or more chunks only. Provide a comma-separated list. [default=%default]" )

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
                                "optimise",
                                "convert",
                                "mst", 
                                "align",
                                "cluster", 
                                "realign",
                                "families", 
                                "summary"),
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
                        alignment_format = "pairsdb",
                        force = False,
                        append = False,
                        test = None,
                        num_jobs = None,
                        temporary_directory = ".",
                        chunks = "all",
                        )
    
    (options, args) = E.Start( parser )

    # setup logging
    if options.loglevel == 0:
        lvl = logging.ERROR
    elif options.loglevel == 1:
        lvl = logging.INFO
    else:
        lvl = logging.DEBUG

    logQueue = multiprocessing.Queue(100)
    handler = Logger.MultiProcessingLogHandler(logging.FileHandler( "adda.log", "w"), logQueue)
    handler.setFormatter( 
        logging.Formatter( '%(asctime)s pid=%(process)-8d %(name)-12s %(levelname)-8s %(message)s',
                           datefmt='%m-%d %H:%M' ) )
    logging.getLogger('adda').addHandler(handler)
    logging.getLogger('adda').setLevel( lvl )

    E.setLogger( logging.getLogger( "adda" ) )
    L = logging.getLogger( "adda" ) 

    config = AddaIO.ConfigParser()
    config.read( os.path.expanduser( options.filename_config ) )

    if args: options.steps = args
        
    ## collect modules and initialise them         
    map_module = { 'fit' : AddaFit.AddaFit,
                   'segment' : AddaSegment.AddaSegment,
                   'blast' : AddaBlast.AddaBlast,
                   'graph' : AddaGraph.AddaGraph,
                   'profiles' : AddaProfiles.AddaProfiles, 
                   'realign' : AddaRealignment.AddaRealignment,
                   'index' : AddaIndex.AddaIndexBuild,
                   'check-index' : AddaIndex.AddaIndexCheck,
                   'optimise' : AddaOptimise.AddaOptimise,  
                   'sequences' : AddaSequences.AddaSequences,
                   'convert' : AddaConvert.AddaConvert,
                   'mst' : AddaMst.AddaMst, 
                   'align' : AddaAlign.AddaAlign, 
                   'cluster' : AddaCluster.AddaCluster,
                   'families' : AddaFamilies.AddaFamilies,
                   'summary' : AddaSummary.AddaSummary,
                   }



    # modules and their hierarchy
    run( options, 
         order = ( "sequences", ), 
         map_module = map_module,
         config = config )

    mFilenameGraph = config.get( "files", "output_graph", "adda.graph")
    mFilenameIndex = config.get( "files", "output_index", "adda.graph.index")

    #i = cadda.IndexedNeighbours( mFilenameGraph, mFilenameIndex )
    #n = i.getNeighbours(1)
    #for nn in n: print str(nn)
    #sys.exit(1)



    # build the adda graph
    run( options,
         order = ( "index", ),
         map_module = map_module,
         config = config )

    if options.num_jobs == 1: 
        run_parallel = runSequentially
    else:
        run_parallel = runParallel

    run_on_graph = RunOnGraph( config, options.steps )

    if "realign" in options.steps:
        run_parallel( 
            run_on_graph,
            filename = config.get( "files", "input_graph", "adda.graph" ),
            options = options, 
            order = ("realign", ),
            map_module = map_module,
            config = config )

    run_parallel( 
        run_on_graph,
        options = options, 
        order = ("fit", "segment" ),
        map_module = map_module,
        config = config )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

    if not merge( options,
                  order = ("fit", "segment" ),
                  map_module = map_module,
                  config = config,
                  fasta = fasta ):
        L.info( "graph pre-processing incomplete - will not continue." )
        E.Stop()
        return
        
    run( options, 
         order = ( "optimise", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    L.info( "building domain graph" )

    run( options, 
         order = ( "convert", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    L.info( "computing minimum spanning tree" )
    run( options, 
         order = ( "mst", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    L.info( "alignment of domains" )

    run_parallel( 
        run_on_file,
        options = options, 
        order = ( "align", ),
        map_module = map_module,
        config = config )

    merge( options,
           order = ("align", ),
           map_module = map_module,
           config = config,
           fasta = fasta )

    run( options, 
         order = ( "cluster", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    run( options, 
         order = ( "families", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    run( options, 
         order = ( "summary", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    E.Stop()
    
                
if __name__ == "__main__":    
    sys.exit(main())









