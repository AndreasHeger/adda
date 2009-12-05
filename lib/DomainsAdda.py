## $Id$

## load domains from a domain file

import sys, re ,string, os, gzip

from Domains import Domains
from Pairsdb import *

from TableDomainsADDA import TableDomainsAdda, TableFamiliesAdda
import Adda.AddaIO

class DomainsAdda (Domains):
    """load domains from a domain file.
    """

    name		= "DomainsFile"                 ## name of sender 
    requirements        = ()                            ## modules that should have finished before

    #--------------------------------------------------------------------------------
    def __init__ (self, dbhandle):
        
        self.mLogLevel = 3

	Domains.__init__( self, dbhandle )

        self.mTableFamilies = TableFamiliesAdda( dbhandle, "adda" )
        self.mTableDomains = TableDomainsAdda( dbhandle, "adda" )
        
        self.mTableFamilies.SetName( self.mTableNameFamilies )        
        self.mTableDomains.SetName( self.mTableNameDomains  )

        self.mFileNameDomains = "adda"

#--------------------------------------< end of class definition >-------------------------------

if __name__ == '__main__':
    dbhandle = Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed" 
	sys.exit(1)

    x = DomainsAdda( dbhandle )
    
    x.Process()

                

