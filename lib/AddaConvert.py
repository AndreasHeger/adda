import sys, os, re, time, math, copy
import cadda
from AddaModule import AddaModule

class AddaConvert( AddaModule ):
    """convert the sequence graph to a domain graph and construct a
        minimum spanning tree."""
    
    mName = "Convert"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameGraph = self.mConfig.get( "files", "output_graph", "adda.graph" )
        self.mFilenameDomains = self.mConfig.get( "files", "output_domains", "adda.domains" )
        self.mEvalueThresholdTrustedLinks = float(self.mConfig.get( "convert", "evalue_threshold_trusted_links", -12.0 ))
        self.mFilenameDomainGraph = self.mConfig.get( "files", "output_domain_graph", "adda.domain_graph" )

        cadda.setFilenameGraph( self.mFilenameGraph )
        cadda.setFilenameDomains( self.mFilenameDomains )
        cadda.setLogLevel( self.mLogLevel )
        cadda.setEvalueThresholdTrustedLinks( self.mEvalueThresholdTrustedLinks )                         

        self.mFilenames = (self.mFilenameDomainGraph, )

    def isComplete( self ):
        return os.path.exists( self.mFilenameDomains )

    def startUp():
        if self.isComplete(): return

    def applyMethod(self ):
        """index the graph.        
        """
        
        if self.isComplete(): return

        self.info( "conversion of sequence graph to domain graph started" )
                
        cadda.dump_parameters()

        retval = cadda.convert( self.mFilenameDomainGraph )
        
        if retval == 0:
            self.warn( "domain graph construction failed" )
        else:
            self.info( "domain graph construction success" )

