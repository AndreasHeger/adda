#!/usr/bin/env python

USAGE="""adda.py [OPTIONS] cmds

interface to compute adda
"""

import sys, os, re, time, math, copy, glob, optparse, logging
from multiprocessing import Process, cpu_count, Pool
import multiprocessing

import fileinput
# segfault on my machine without the print statement - strange
print 
import Adda.Experiment as E
from Adda import *

from logging import warn, info, debug

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

    if len(modules) == 0: return

    E.info( "working with modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.startUp()
        module.run()
        module.finish()

def merge( options, 
           order, 
           map_module, 
           config, 
           fasta ):

    if "all" in options.steps: steps = order
    else: steps = options.steps

    modules = []
    for step in order:
        if step in steps:
            modules.append( map_module[step]( options, 
                                              config, 
                                              fasta = fasta ) )

    E.info( "performing merge on modules: %s" % str(modules) )

    for module in modules:
        module.merge()

class Run(object):
    pass

class RunOnGraph(Run):
    
    def __init__(self, config, steps ):

        manager = multiprocessing.Manager()

        dict_mapper = manager.dict

        E.info( "loading fasta sequences from %s" % config.get( "files", "output_fasta", "adda" ) ) 
        self.mFasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

        E.info( "loading map_id2nid from %s" % config.get( "files", "output_nids", "adda.nids" ))
        self.mMapId2Nid = dict_mapper(AddaIO.readMapId2Nid( open(config.get( "files", "output_nids", "adda.nids" ), "r" ) ))

        if "all" in steps or "fit" in steps:
            E.info( "loading domain boundaries from %s" % config.get( "files", "input_reference") )
            infile = AddaIO.openStream( config.get( "files", "input_reference") )
            rx_include = config.get( "fit", "family_include", "") 
            self.mMapNid2Domains = dict_mapper(AddaIO.readMapNid2Domains( infile, self.mMapId2Nid, rx_include ))
            infile.close()
        else:
            self.mMapNid2Domains = None

    def __call__(self, argv ):
        (filename, options, order, map_module, config, chunk, nchunks ) = argv
        E.info( "starting chunk %i on %s" % (chunk, filename) )

        if "all" in options.steps: steps = order
        else: steps = options.steps

        modules = []
        for step in order:
            if step in steps:
                if map_module[step]( options, 
                                     config = config, 
                                     fasta = self.mFasta ).isComplete():
                    E.info( "%s complete" % step )
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
                    E.info( "%s complete" % step )

        if len(modules) == 0: 
            E.info( "nothing to be done" )
            return

        E.info( "opening graph %s at chunk %i" % (filename, chunk) )

        gzip_factor = None
        if filename.endswith(".gz"):
            uncompressed_size = config.get( "files", "input_graph_uncompressed_size", 0 )
            if uncompressed_size > 0:
                compressed_size = FileSlice.getFileSize( filename )
                assert compressed_size < uncompressed_size, "file size of gzipped graph is larger than given uncompressed size" 
                gzip_factor = float(compressed_size) / uncompressed_size
                E.info( "setting graph compression to %5.2f (compressed = %i / uncompressed = %i)" % (gzip_factor, compressed_size, uncompressed_size) )

                assert 0.1 < gzip_factor < 0.8, "gzip factor unrealistic - values between 0.1 and 0.8 are usual, but is %f" % gzip_factor

        iterator = FileSlice.IteratorMultiline( filename, 
                                                nchunks,
                                                chunk,
                                                FileSlice.groupby,
                                                key = lambda x: x[:x.index("\t")],
                                                gzip_factor = gzip_factor )




        if options.alignment_format == "pairsdb":
            record_type = AddaIO.NeighbourRecordPairsdb
        elif options.alignment_format == "pairsdb-old":
            record_type = AddaIO.NeighbourRecordPairsdbOld
        elif options.alignment_format == "simap":
            record_type = AddaIO.NeighbourRecordSimap
        elif options.alignment_format == "pairsdb-realign":
            record_type = AddaIO.NeighbourRecordPairsdbRealign
        else:
            raise ValueError ("unknown record type %s" % options.alignment_format)

        for module in modules: module.startUp()

        E.info( "starting work on modules: %s" % (",".join(map(str, modules))) )

        map_id2nid = self.mMapId2Nid

        for record in iterator:
            neighbours = []
            for line in record:
                n = record_type( line )
                if map_id2nid:
                    if (n.mQueryToken not in map_id2nid or \
                            n.mSbjctToken not in map_id2nid ):
                        continue 
                    q = n.mQueryToken = map_id2nid[n.mQueryToken]
                    n.mSbjctToken = map_id2nid[n.mSbjctToken]

                neighbours.append( n )

            E.info( "started: nid=%s, neighbours=%i" % (str(q), len(neighbours) ) )

            if neighbours:
                for module in modules:
                    module.run( AddaIO.NeighboursRecord( q, neighbours ) )

            E.info( "finished: nid=%s, neighbours=%i" % (str(q), len(neighbours) ) )

        E.info( "running finish on modules: %s" % (",".join(map(str, modules))) )

        for module in modules:
            module.finish()

        E.info( "finished chunk %i on %s" % (chunk, filename) )
        

def old_run_on_graph( argv ):
    """process graph."""

    (filename, options, order, map_module, config, chunk, nchunks ) = argv
    E.info( "starting chunk %i on %s" % (chunk, filename) )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

    if "all" in options.steps: steps = order
    else: steps = options.steps

    modules = []
    for step in order:
        if step in steps:
            if map_module[step]( options, 
                                 config = config, 
                                 fasta = fasta ).isComplete():
                E.info( "%s complete" % step )
                continue

            module = map_module[step]( options, 
                                       config = config, 
                                       fasta = fasta,
                                       num_chunks = nchunks,
                                       chunk = chunk )

            if not module.isComplete():
                modules.append( module )
            else:
                E.info( "%s complete" % step )

    if len(modules) == 0: 
        E.info( "nothing to be done" )
        return

    E.info( "opening graph %s at chunk %i" % (filename, chunk) )

    gzip_factor = None
    if filename.endswith(".gz"):
        uncompressed_size = config.get( "files", "input_graph_uncompressed_size", 0 )
        if uncompressed_size > 0:
            compressed_size = FileSlice.getFileSize( filename )
            assert compressed_size < uncompressed_size, "file size of gzipped graph is larger than given uncompressed size" 
            gzip_factor = float(compressed_size) / uncompressed_size
            E.info( "setting graph compression to %5.2f (compressed = %i / uncompressed = %i)" % (gzip_factor, compressed_size, uncompressed_size) )

            assert 0.1 < gzip_factor < 0.8, "gzip factor unrealistic - values between 0.1 and 0.8 are usual, but is %f" % gzip_factor

    iterator = FileSlice.IteratorMultiline( filename, 
                                            nchunks,
                                            chunk,
                                            FileSlice.groupby,
                                            key = lambda x: x[:x.index("\t")],
                                            gzip_factor = gzip_factor )


    map_id2nid = AddaIO.readMapId2Nid( open(config.get( "files", "output_nids", "adda.nids" ), "r" ) )

    if options.alignment_format == "pairsdb":
        record_type = AddaIO.NeighbourRecordPairsdb
    elif options.alignment_format == "pairsdb-old":
        record_type = AddaIO.NeighbourRecordPairsdbOld
    elif options.alignment_format == "simap":
        record_type = AddaIO.NeighbourRecordSimap
    elif options.alignment_format == "pairsdb-realign":
        record_type = AddaIO.NeighbourRecordPairsdbRealign
    else:
        raise ValueError ("unknown record type %s" % options.alignment_format)

    for module in modules: module.startUp()

    E.info( "starting work on modules: %s" % (",".join(map(str, modules))) )

    for record in iterator:
        neighbours = []
        for line in record:
            n = record_type( line )
            if map_id2nid:
                if (n.mQueryToken not in map_id2nid or \
                        n.mSbjctToken not in map_id2nid ):
                    continue 
                q = n.mQueryToken = map_id2nid[n.mQueryToken]
                n.mSbjctToken = map_id2nid[n.mSbjctToken]

            neighbours.append( n )

        E.debug( "working on: %s with %i neighbours" % (str(q), len(neighbours) ) )

        if neighbours:
            for module in modules:
                module.run( AddaIO.NeighboursRecord( q, neighbours ) )

    E.info( "running finish on modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.finish()

    E.info( "finished chunk %i on %s" % (chunk, filename) )

def run_on_file( argv ):

    (filename, options, order, map_module, config, chunk, nchunks ) = argv

    E.info( "starting chunk %i on %s" % (chunk, filename) )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

    if "all" in options.steps: steps = order
    else: steps = options.steps

    modules = []
    for step in order:
        if step in steps:
            if map_module[step]( options, 
                                 config = config, 
                                 fasta = fasta ).isComplete():
                E.info( "%s complete" % step )
                continue

            module = map_module[step]( options, 
                                       config = config, 
                                       fasta = fasta,
                                       num_chunks = nchunks,
                                       chunk = chunk )

            if not module.isComplete():
                modules.append( module )
            else:
                E.info( "%s complete" % step )

    if len(modules) == 0: 
        E.info( "nothing to be done" )
        return

    E.info( "opening file %s at chunk %i" % (filename, chunk) )

    iterator = FileSlice.Iterator( filename, 
                                   nchunks,
                                   chunk,
                                   FileSlice.iterator )

    for module in modules: module.startUp()

    E.info( "working with modules: %s" % (",".join(map(str, modules))) )

    for line in iterator:
        if line.startswith("#"): continue

        for module in modules:
            module.run( line )

    for module in modules:
        module.finish()

    E.info( "running finish on modules: %s" % (",".join(map(str, modules))) )

    for module in modules:
        module.finish()

    E.info( "finished chunk %i on %s" % (chunk, filename) )

def runParallel( runner, filename, options, order, map_module, config ):
    """process filename in paralell."""

    nchunks = config.get( "adda", "num_slices", 10 )
    if options.num_jobs:
        njobs = options.num_jobs 
    else:
        njobs = cpu_count()
    
    E.info( "running %i jobs on %i slices" % (njobs, nchunks ))

    pool = Pool( njobs )

    logging.info('starting parallel jobs')

    args = [ (filename, options, order, map_module, config, chunk, nchunks ) for chunk in range(nchunks) ]

    pool.map( runner, args )
    pool.close()
    pool.join()

    E.info( "all jobs finished" )

def runSequentially( runner, filename, options, order, map_module, config ):
    """process filename sequentially."""

    nchunks = config.get( "adda", "num_slices", 4 )

    args = [ (filename, options, order, map_module, config, chunk, nchunks ) for chunk in range(nchunks) ]

    for (chunck, argv) in enumerate(args):
        E.info( "job %i started" % chunk )
        runner( argv )
        E.info( "job %i finished" % chunk )

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

    parser.add_option( "--num-jobs", dest="num_jobs", type="int",
                      help="use # processes. If not set, the number of CPUs/cores is taken [default=%default]")
    
    parser.add_option( "--alignment-format", dest="alignment_format", type="choice",
                       choices=("pairsdb", "pairsdb-old", "simap", "pairsdb-realign"),
                       help = "input format of graph. pairsdb-old: input graph is in old 1-based coordinates [default=%default]." )

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
    handler= Logger.MultiProcessingLogHandler(logging.FileHandler( "adda.log", "w"), logQueue)
    logging.getLogger('adda').addHandler(handler)
    logging.getLogger('adda').setLevel( lvl )

    E.setLogger( logging.getLogger( "adda" ) )

    #logging.basicConfig(
    #    lvl = lvl,
    #    format='%(asctime)s %(name)s %(levelname)s %(message)s',
    #    filename = "adda.log" )

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
        filename = config.get( "files", "input_graph", "adda.graph" ),
        options = options, 
        order = ("fit", "segment", "graph", "profiles" ),
        map_module = map_module,
        config = config )

    fasta = IndexedFasta.IndexedFasta( config.get( "files", "output_fasta", "adda" ) )

    merge( options,
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

    E.info( "computing minimum spanning tree" )
    run( options, 
         order = ( "mst", ),
         map_module = map_module,
         config = config,
         fasta = fasta)

    E.info( "alignment of domains" )

    run_parallel( 
        run_on_file,
        config.get( "files", "output_mst", "adda.mst" ),
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









