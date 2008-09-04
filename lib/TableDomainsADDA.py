####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id$
##
##
####
####

import string
import Pairsdb
from TableDomains import TableDomains, TableFamilies

class TableDomainsAdda( TableDomains ):

    mTypeDomainId       = 'VARCHAR(30) DEFAULT ""'
    mTypeDomainClass    = 'INT UNSIGNED NOT NULL DEFAULT 0'
    mExtraFields        = (
        )
    
    mExtraIndices       = ()
    
    def __init__ ( self, handle, root = "adda" ):
        
	TableDomains.__init__( self, handle, root )

        ## select statement for additional info
        self.mAdditionalInfo = ""

    #---------------------------------------------------------------------------------------------------------------
    def GetAnotherInstance( self ):
        """return a handle to the same table."""
        return TableDomainsAdda( self.dbhandle )

    ##-----------------------------------------------------------------------------------------------
    def GetMaxFamily( self ):
        """retrieve maximum family."""
        return self.Execute( "SELECT MAX(family) FROM " + self.name).fetchone()[0]

    ##-----------------------------------------------------------------------------------------------
    def AddDomainsFromDomainsTable( self, table_name_source, offset = 0, extra_fields = "", subset = None):
        """adds domains from another table. Uses only the minimal information.
        Adds an offset to the domain family
        """
        if subset:
            s = " INNER JOIN %s AS subset ON subset.nid = rep_nid " % subset
        else:
            s = ""
            
        statement = """
        INSERT INTO %s
        SELECT
        rep_nid, rep_from, rep_to, rep_ali,
        domain_id, domain_from, domain_to, domain_ali,
        family + %i
        %s
        FROM %s
        %s
        """ % (self.name, offset, extra_fields, table_name_source, s)

        return self.Execute(statement)

##------------------------------------------------------------------------------------------------------
class TableFamiliesAdda( TableFamilies ):

    mTypeDomainClass    = 'INTEGER UNSIGNED NOT NULL DEFAULT 0'

    mExtraFields  = ()
    mExtraIndices = ()
    
    def __init__ ( self, handle, root = "adda" ):

	TableFamilies.__init__( self, handle, root )

    ##----------------------->start: common methods<----------------------
	
    def GetAnotherInstance( self ):
        """return a handle to the same table."""
        return TableFamiliesAdda( self.dbhandle )

    ##-----------------------------------------------------------------------------------------------
    def GetMaxFamily( self ):
        """retrieve maximum family."""
        return self.Execute( "SELECT MAX(family) FROM " + self.name).fetchone()[0]
