import sys, os, re, time, math, copy, glob, optparse, math, subprocess, tempfile
import shutil

import cadda

from AddaModule import AddaModuleBlock

class AddaMst( AddaModuleBlock ):
    """convert the sequence graph to a domain graph and construct a
        minimum spanning tree."""
    
    mName = "Mst"
    
    def __init__(self, *args, **kwargs ):

        AddaModuleBlock.__init__( self, *args, **kwargs )
                
        self.mFilenameDomainGraph = self.mConfig.get( "files", "output_domain_graph", "adda.domain_graph" )
        self.mFilenameMst = self.mConfig.get( "files", "output_mst", "adda.mst" )
                
        cadda.setLogLevel( self.mLogLevel )
        # cadda.setReportStep( 1 )

        self.mFilenames = ( self.mFilenameMst, )

    def startUp( self ):
        if self.isComplete(): return

    def applyMethod(self ):
        """index the graph.        
        """
        
        if self.isComplete(): return

        self.info( "construction of minimum spanning tree started" )
                
        cadda.dump_parameters()
        
        tmpdir = tempfile.mkdtemp( dir = self.mTemporaryDirectory )
        tmpfile = os.path.join( tmpdir, "sorted" )

        if not os.path.exists( tmpfile ):
            statement = "cut -f 1-3 %s | sort -T%s -k3,3n > %s" % ( self.mFilenameDomainGraph, 
                                                                    tmpdir,                                                  
                                                                    tmpfile ) 
        
            self.info( "sorting started" )
                        
            try:
                retcode = subprocess.call( statement , shell=True)
                if retcode < 0:
                    self.warn( "sorting was terminated by signal %i" % (-retcode) )
                elif retcode > 0:
                    self.warn( "sorting returned %i" % (retcode) )                
            except OSError, e:
                self.warn( "sorting failed with message: %s" % (e) )
            
            self.info( "sorting finished" )                
        else:
            self.info( "skipping sorting, because sorted output already exists" )                
            
        noutput = cadda.build_mst( self.mFilenameMst, tmpfile )
        
        if noutput == 0:
            self.warn( "mst construction failed" )
        else:
            self.info( "mst construction success: %i links output" % noutput )

        shutil.rmtree( tmpdir )
