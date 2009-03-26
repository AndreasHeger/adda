import sys, os, re, time, math, copy, glob, optparse, math
import pylab

import cadda

from AddaModule import AddaModule
import AddaIO
import Components
import SegmentedFile

class AddaCluster( AddaModule ):
    """Do quality control. Compare adda families to references."""
    
    mName = "Cluster"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameFamilies = self.mConfig.get( "files", "output_families", "adda.families" )
        self.mFilenameAlignments = self.mConfig.get("files","output_align", "adda.align" )
        self.mFilenames = (self.mFilenameFamilies, )

        self.mMinAlignedResidues = self.mConfig.get("cluster", "min_aligned_residues", 30 )
        self.mPatternFamily = self.mConfig.get("cluster", "pattern_family", "AD%06i" )

    def startUp(self):
        if self.isComplete(): return
        self.mOutfile = self.openOutputStream( self.mFilenameFamilies )

    def applyMethod(self ):
        """index the graph.        
        """
        componentor = Components.SComponents()

        infile = SegmentedFile.openfile( self.mFilenameAlignments, "r" )

        naccepted, nrejected_score, nrejected_aligned = 0, 0, 0
        for line in infile:
            if line[0] == "#": continue
            if line.startswith( "passed"): continue
            
            (code, qdomain, sdomain, estimate, 
             qstart, qend, qali, sstart, send, sali, 
             score, naligned, ngaps, zscore) =\
                line[:-1].split("\t")
            
            if code == "+":
                if int(ngaps) >= self.mMinAlignedResidues:
                    componentor.add( qdomain, sdomain )
                    naccepted += 1
                else:
                    nrejected_aligned += 1
            nrejected_score += 1
        
        self.info( "computing components with %i accepted links (%i rejected score, %i rejected alignment length)" %\
                   (naccepted, nrejected_score, nrejected_aligned ) )
        
        components = componentor.getComponents()
        self.mOutfile.write( "nid\tstart\tend\tfamily\n" )

        noutput = 0
        family_id = 0 
        for domains in components:
            family_id += 1
            for domain in domains:
                nid, start, end = domain.split("_")
                self.mOutfile.write( "%s\t%s\t%s\t%s\n" % \
                                     ( nid, start, end, self.mPatternFamily % family_id ) )

                noutput += 1

        self.mOutfile.close()

        self.info( "output %i families for %i domains" % (family_id, noutput) )
               
