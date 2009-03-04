import sys, os, re, time

import alignlib

from AddaModule import AddaModule
import IndexedFasta

class AddaBlast( AddaModule ):
    """run blast."""
    
    mName = "Blast"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mFilenameOutputFasta = self.mConfig.get( "files", "output_fasta", "adda" )
        self.mBlastDatabase = self.mConfig.get( "files", "database_blast", "adda" )
        self.mBlastCPUs = self.mConfig.get( "blast", "num_cpus", 2 )
        self.mBlastEvalue = self.mConfig.get( "blast", "evalue", 1.0 )
        self.mBlastNumResults = self.mConfig.get( "blast", "num_results", 100000 )
        self.mBlastResults = self.mConfig.get( "files", "output_blast", "adda.blast" )

    def applyMethod(self ):
        
        cmd = "formatdb -i %s.fasta -p T -n %s" % (self.mFilenameOutputFasta, self.mBlastDatabase )

        self.execute( cmd )

        cmd = "blastall -p blastp -i %s.fasta -d %s -a %i -e %f -b %i -v %i -m 8 > %s" %\
            (self.mFilenameOutputFasta, 
             self.mBlastDatabase,
             self.mBlastCPUs,
             self.mBlastEvalue,
             self.mBlastNumResults,
             self.mBlastNumResults,
             self.mBlastResults )

        self.execute( cmd )
