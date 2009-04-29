## $Id$

## load domains from a domain file

import sys, re ,string, os, gzip

from Domains import Domains
from Pairsdb import *

from TableDomainsADDA import TableDomainsAdda, TableFamiliesAdda
from TableDomainsAddaMalis import TableDomainsAddaMalis, TableFamiliesAddaMalis

class DomainsAdda (Domains):
    """load domains from a domain file.
    """

    name		= "DomainsFile"                 ## name of sender 
    requirements        = ()                            ## modules that should have finished before

    #--------------------------------------------------------------------------------
    def __init__ (self, dbhandle):
        
        # log parameters
        self.mLogLevel = 3

        self.mShortOptions  += 'i:'
        self.mLongOptions   += ['input=']

	Domains.__init__( self, dbhandle )

        self.mFileNameFamilies = None
        self.mFileNameDomains = os.tempnam( PATH_LOAD, "doma" )

        self.mTableFamilies = TableFamiliesAdda( dbhandle, "adda" )
        self.mTableDomains = TableDomainsAdda( dbhandle, "adda" )
        
        self.mTableFamilies.SetName( self.mTableNameFamilies )        
        self.mTableDomains.SetName( self.mTableNameDomains  )

    ##---------------------------------------------------------------------------------
    def ProcessOptions( self, optlist ):

        Domains.ProcessOptions( self, optlist)

        for o,a in optlist:        
            if o in ("-i", "--input"):
                self.mFileNameInput = a

    ##---------------------------------------------------------------------------------
    def Open( self, filename, mode= "r" ):
        """return opened file."""
        if filename.endswith(".gz"):
            return gzip.open( filename, "r" )
        else:
            return open( filename, "r" )

    ##---------------------------------------------------------------------------------
    def Create( self ):
        """select all members of a connected component and
        add to families. Each connected component gives the
        component_id for a given family.
        
        """

        self.mTableFamilies.Clear()
        self.mTableDomains.Clear()

        if self.mLogLevel >= 1:
            print "--> loading data"            
            print "--> at start: %i families with %i members" % (self.mTableFamilies.RowCount(),
                                                                 self.mTableDomains.RowCount())
            sys.stdout.flush()

        self.OpenOutfiles()

        infile = self.Open(self.mFileNameInput, "r")

        families = {}

        for line in infile:
            if line.startswith("#"): continue
            if line.startswith("nid"): continue
            (domain_nid, domain_from, domain_to, family) = line[:-1].split("\t")
            domain_from, domain_to = map( int, (domain_from, domain_to) )
            if not families.has_key(family):
                self.mTableFamilies.AddFamily(family)
                families[family] = 1

            l = domain_to - domain_from + 1
            self.mFileDomains.write( string.join( map(str, (domain_nid, domain_from, domain_to, "+%i" % l,
                                                            domain_nid, domain_from, domain_to, "+%i" % l,
                                                            family
                                                            )), "\t") + "\n")

        self.CloseOutfiles()
        self.Load()
        
        if self.mLogLevel >= 1:
            
            print "--> at the end: %i families with %i members" % (self.mTableFamilies.RowCount(),
                                                                   self.mTableDomains.RowCount())
            sys.stdout.flush()


    ##---------------------------------------------------------------------------------
    def CreateMalis( self ):
        """Build familes from a links file.
        """

        self.mTableFamilies = TableFamiliesAddaMalis( dbhandle, "adda" )
        self.mTableDomains = TableDomainsAddaMalis( dbhandle, "adda" )
        
        self.mTableFamilies.SetName( self.mTableNameFamilies )        
        self.mTableDomains.SetName( self.mTableNameDomains  )

        self.mTableFamilies.Clear()
        self.mTableDomains.Clear()

        if self.mLogLevel >= 1:
            print "--> loading data"            
            print "--> at start: %i families with %i members" % (self.mTableFamilies.RowCount(),
                                                                 self.mTableDomains.RowCount())
            sys.stdout.flush()

        ## simply create a copy
        os.system("cp %s %s" % (self.mFileNameInput, self.mFileNameDomains ))

        self.Load()
        self.mTableFamilies.AddDomains( self.mTableDomains )
        self.UpdateDomains()
        
        if self.mLogLevel >= 1:
            
            print "--> at the end: %i families with %i members" % (self.mTableFamilies.RowCount(),
                                                                   self.mTableDomains.RowCount())
            sys.stdout.flush()

#--------------------------------------< end of class definition >-------------------------------

if __name__ == '__main__':
    dbhandle = Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed" 
	sys.exit(1)

    x = DomainsAdda( dbhandle )
    
    x.Process()

                

