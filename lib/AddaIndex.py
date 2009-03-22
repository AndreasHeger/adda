import sys, os, re

import cadda

from AddaModule import AddaModule

class AddaIndex( AddaModule ):
    """index graph for adda processing."""
    
    mName = "Index"
    
    def __init__(self, *args, **kwargs ):
        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameGraph = self.mConfig.get( "files", "output_graph", "adda.graph")
        self.mFilenameIndex = self.mConfig.get( "files", "output_index", "adda.graph.index")
                        
        cadda.setFilenameGraph( self.mFilenameGraph )
        cadda.setFilenameIndex( self.mFilenameIndex )
        cadda.setLogLevel( self.mLogLevel )
        cadda.setReportStep( self.mConfig.get( "adda", "report_step", 1000 ) )
        cadda.dump_parameters()
        
        self.mFilenames = (self.mFilenameIndex, )

class AddaIndexBuild( AddaIndex ):
    """index a graph."""
    
    mName = "BuildIndex"
    
    def __init__(self, *args, **kwargs ):
        AddaIndex.__init__( self, *args, **kwargs )

        self.mRequirements.append( self.mFilenameIndex )
                        
    def applyMethod(self ):
        """index the graph.        
        """
                
        self.info( "indexing of %s started" % self.mFilenameGraph )
        retval = cadda.build_index()
        if retval != 0:
            self.warn( "indexing of %s failed with error %i" % (self.mFilenameGraph, retval))
        else:
            self.info( "indexing of %s success" % (self.mFilenameGraph,))        

        
class AddaIndexCheck( AddaIndex ):
    """check indexed graph."""
    
    mName = "CheckIndex"
    
    def __init__(self, *args, **kwargs ):
        AddaIndex.__init__( self, *args, **kwargs )
                        
    def applyMethod(self ):
        """index the graph.        
        """
                
        self.info( "checking index %s agains %s: started" % (self.mFilenameIndex, self.mFilenameGraph ) )
        retval = cadda.check_index()
        if retval == 1:
            self.warn( "index check of %s failed" % self.mFilenameGraph)
        else:
            self.info( "indexing of %s passed" % (self.mFilenameGraph)) 

        
