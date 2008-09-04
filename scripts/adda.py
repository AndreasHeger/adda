#!/usr/bin/env python

USAGE="""adda.py [OPTIONS] cmds

interface to compute adda
"""

import sys, os, re, time, math, copy, glob, optparse
import fileinput

from Adda import *

if __name__ == "__main__":
    
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
                                "fit", 
                                "graph",
                                "index",
                                "check-index",
                                "profiles", 
                                "segment", 
                                "optimise",
                                "convert",
                                "mst", 
                                "align", ),
                       help="perform this step [default=%default]" )

    parser.add_option( "--start-at", dest="start_at", type="string",
                      help="start at sequence [default=%default]")

    parser.add_option( "--stop-at", dest="stop_at", type="string",
                      help="stop at sequenec [default=%default]")


    parser.set_defaults( 
                        filename_config = "adda.ini",
                        steps = [],
                        suffix = "",
                        start_at = None,
                        stop_at = None,
                        old_alignment_format = False,
                        force = False,
                        append = False,
                        test = None,
                        temporary_directory = ".",
                        )
    
    (options, args) = Experiment.Start( parser )

    config = AddaIO.ConfigParser()
    config.read( os.path.expanduser( options.filename_config ) )

    if args: options.steps = args
        
    ## collect modules and initialise them         
    map_module = { 'fit' : AddaFit.AddaFit,
                   'segment' : AddaSegment.AddaSegment,
                   'graph' : AddaGraph.AddaGraph,
                   'profiles' : AddaProfiles.AddaProfiles, 
                   'index' : AddaIndex.AddaIndexBuild,
                   'check-index' : AddaIndex.AddaIndexCheck,
                   'optimise' : AddaOptimise.AddaOptimise,                    
                   'convert' : AddaConvert.AddaConvert,
                   'mst' : AddaMst.AddaMst, 
                   'align' : AddaAlign.AddaAlign, 
                   }
    
    # modules and their hierarchy
    pre_modules = ("fit", "segment", "graph", "profiles" )
    post_modules = ("index", "check-index", "optimise", "convert", "mst", "align" )
     
    if "all" in options.steps:
        options.steps = pre_modules + post_modules
        
    if len(options.steps) == 0:
        raise "nothing asked to be done."
        
    fasta = IndexedFasta.IndexedFasta( config.get( "files", "input_fasta") )
    tokens = set(fasta.getContigSizes().keys())
    
    modules_on_graph = [] 
    modules_on_adda = []
    for step in pre_modules:
        if step in options.steps:
            modules_on_graph.append( map_module[step]( options, config, fasta ) )
    for step in post_modules:    
        if step in options.steps:                
            modules_on_adda.append( map_module[step]( options, config, fasta ) )
            
    options.stdlog.write( "modules on graph: %s\n" % (",".join(map(str, modules_on_graph))))
    options.stdlog.write( "modules on adda: %s\n" % (",".join(map(str, modules_on_adda))))
            
    if modules_on_graph:
        
        files = config.get( "files", "input_graph").split(",")
        
        infile = fileinput.FileInput(files = files,
                                     openhook=fileinput.hook_compressed)
    
        iterator = AddaIO.NeighboursIterator(infile,
                                             tokens = tokens )
    
        ninput, noutput = 0, 0
        t0 = time.time()

        keep = options.start_at == None
        
        while 1:
            
            t1 = time.time()
            
            neighbours = iterator.next()
            
            if neighbours == None: break
        
            if options.start_at and neighbours.mQueryToken == options.start_at: 
                keep = True
                
            if not keep: 
                if options.loglevel >= 3:
                    options.stdlog.write( "# skipping: %s, will start at %s\n" %\
                                              ( neighbours.mQueryToken, options.start_at) )
                continue
                
            if options.stop_at and neighbours.mQueryToken == options.stop_at: break
            
            ninput += 1
        
            if options.old_alignment_format:
                for match in neighbours.mMatches:
                    match.mQueryFrom -= 1
                    match.mSbjctFrom -= 1
            
            for module in modules_on_graph:
                module.apply( neighbours )
    
            noutput += 1
                
            t2 = time.time()
            
            options.stdlog.write( "# iteration=%i, token=%s, time: this=%i, total=%i, avg=%f\n" %\
                                  (ninput, neighbours.mQueryToken, t2-t1, t2-t0, float(t2-t0) / ninput) )

            if options.test and ninput >= options.test:
                break
    
        for module in modules_on_graph:
            module.finish()
    
        if options.loglevel >= 1:
            options.stdlog.write( "# ninput=%i, noutput=%i\n" % (ninput, noutput ) )
        
    if modules_on_adda:
        
        for module in modules_on_adda:
            module.run()
            module.finish()
        
    Experiment.Stop()
    
                
    










