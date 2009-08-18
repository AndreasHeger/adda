import sys, os, re

import cadda

from AddaModule import AddaModuleBlock
import AddaIO

class AddaIndex( AddaModuleBlock ):
    """index graph for adda.

    input
       ``files:input_graph``: the input graph
    output
       ``files:output_graph``: the indexed and compressed pairwise alignment graph.
       ``files:output_index``: an index of the pairwise alignment graph.
         
    Both output files are compressed.
    """
    
    mName = "Index"
    
    def __init__(self, *args, **kwargs ):
        AddaModuleBlock.__init__( self, *args, **kwargs )
                
        self.mFilenameInputGraph = self.mConfig.get( "files", "input_graph", "adda.graph")

        self.mFilenameOutputGraph = self.mConfig.get( "files", "output_graph", "adda.graph")
        self.mFilenameOutputIndex = self.mConfig.get( "files", "output_index", "adda.graph.index")
                        
        cadda.setLogLevel( self.mLogLevel )
        cadda.setReportStep( self.mConfig.get( "adda", "report_step", 1000 ) )
        cadda.dump_parameters()
        
        self.mFilenames = (self.mFilenameOutputIndex, )

        self.mAlignmentFormat = self.mConfig.get( "files", "graph_format", "pairsdb")

        if self.mAlignmentFormat == "pairsdb":
            self.mIterator = AddaIO.NeighboursIteratorPairsdb
        elif self.mAlignmentFormat == "pairsdb-old":
            self.mIterator = AddaIO.NeighboursIteratorPairsdbOld
        elif self.mAlignmentFormat == "simap":
            self.mIterator = AddaIO.NeighboursIteratorSimap
        elif self.mAlignmentFormat == "pairsdb-realign":
            self.mIterator = AddaIO.NeighbourRecordPairsdbRealign
        else:
            raise ValueError ("unknown record type %s" % self.mAlignmentFormat)

    def startUp( self ):
        pass

class AddaIndexBuild( AddaIndex ):
    """index a graph."""
    
    mName = "BuildIndex"
    
    def __init__(self, *args, **kwargs ):
        AddaIndex.__init__( self, *args, **kwargs )
        
    def isComplete( self ):
        return os.path.exists( self.mFilenameOutputIndex ) and os.path.exists( self.mFilenameOutputGraph)

    def applyMethod(self ):
        """index the graph.        
        """
        self.info( "indexing of %s started" % self.mFilenameInputGraph )

        self.info( "loading map_id2nid from %s" % self.mConfig.get( "files", "output_nids", "adda.nids" ))
        infile = open( self.mConfig.get( "files", "output_nids", "adda.nids" ) )
        map_id2nid = AddaIO.readMapId2Nid( infile, 
                                           storage = self.mConfig.get( "files", "storage_nids", "memory" ) )
        infile.close()
    
        infile = AddaIO.openStream( self.mFilenameInputGraph )
        iterator = self.mIterator( infile, map_id2nid )
        cadda.indexGraph( iterator, len(map_id2nid), self.mFilenameOutputGraph, self.mFilenameOutputIndex, self.mLogger )

        del map_id2nid

class AddaIndexCheck( AddaIndex ):
    """check indexed graph."""
    
    mName = "CheckIndex"
    
    def __init__(self, *args, **kwargs ):
        AddaIndex.__init__( self, *args, **kwargs )
                        
    def applyMethod(self ):
        """index the graph.        
        """

        if os.path.exists( filename_graph ): pass


        self.info( "checking index %s agains %s: started" % (self.mFilenameIndex, self.mFilenameGraph ) )
        retval = cadda.check_index()
        if retval == 1:
            self.warn( "index check of %s failed" % self.mFilenameGraph)
        else:
            self.info( "indexing of %s passed" % (self.mFilenameGraph)) 

        
