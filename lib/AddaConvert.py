import sys, os, re, time, math, copy
import cadda
from AddaModule import AddaModule

class AddaConvert( AddaModule ):
    """convert the sequence graph to a domain graph and construct a
        minimum spanning tree."""
    
    mName = "Convert"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameGraph = self.mConfig.get( "files", "output_graph")
        self.mFilenameDomains = self.mConfig.get( "files", "output_domains" )
        self.mFilenameDomainGraph = self.mConfig.get( "files", "output_domain_graph" )
        self.mEvalueThresholdTrustedLinks = float(self.mConfig.get( "convert", "evalue_threshold_trusted_links" ))

        self.mRequirements.append( self.mFilenameGraph )
        self.mRequirements.append( self.mFilenameDomains )
            
        cadda.setFilenameGraph( self.mFilenameGraph )
        cadda.setFilenameDomains( self.mFilenameDomains )
        cadda.setLogLevel( self.mOptions.loglevel )
        cadda.setEvalueThresholdTrustedLinks( self.mEvalueThresholdTrustedLinks )                         
        
    def applyMethod(self ):
        """index the graph.        
        """
        
        self.info( "conversion of sequence graph to domain graph started" )
                
        cadda.dump_parameters()

        retval = cadda.convert( self.mFilenameDomainGraph )
        
        if retval == 0:
            self.warn( "domain graph construction failed" )
        else:
            self.info( "domain graph construction success" )

