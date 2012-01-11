####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Table_nrdb40.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####


#----------------------------------------------------------------
# Name:		Table_nrdb40
#--------------General Information-------------------------------
# File:		Table_nrdb40.py
# Version:		1.0
# Description:	Table-Object, created automatically by Pairsdb
# Author:		Andreas Heger (heger@ebi.ac.uk)
#--------------Documentation-------------------------------------
# 
#
#--------------Change History------------------------------------
# 20.01.2000        	Created
#
#----------------------------------------------------------------

import Pairsdb
from Table import Table

class TableNids( Table ):
    def __init__ ( self, handle ):
	self.name   =  'nrdb'
        
	self.fields = (
			('nid', 'INT UNSIGNED NOT NULL'),
		)
	self.indices = (
			'UNIQUE (nid)',
		)
	Table.__init__( self, handle )

    #---------------------------------------------------------------------------------
    def GetLongSequences( self, size ):
        """retrieve all sequence above a maximum size.
        """
        statement = "SELECT n.nid FROM %s AS n, nrdb AS r" % self.name +\
                    " WHERE r.nid = n.nid AND r.length > %i" % size
        return map( lambda x:x[0], self.Execute(statement).fetchall())

    #--------------------------------------------------------------------------------
    def GetAllNids( self ):
        """retrieve all nids"""
        statement = "SELECT DISTINCTROW nid FROM %s " % self.name
        return map(lambda x: x[0], self.Execute( statement).fetchall())

    ##-------------------------------------------------------------
    def GetRandomSample( self, sample_size = 1):
        """retrieve a random nid"""
        statement = "SELECT nid FROM %s ORDER BY RAND() LIMIT %i" % (self.name, sample_size)
        return map( lambda x: x[0], self.Execute(statement).fetchall())
    
    #--------------------------------------------------------------------------------
    def SelectSequences( self ):
        result = self.Execute( "SELECT r.nid, r.sequence FROM %s.nrdb AS r, %s AS n" % (self.mPrefix, self.name) +\
                              " WHERE r.nid = n.nid")

        return result

    ##-----------------------------------------------------------------------------------------	
    def GetNids(self ):
        """return a tupe of all nids in database."""
        statement = "SELECT DISTINCTROW nid FROM " + self.name 
        query = self.Execute( statement)

        return tuple(map( lambda x: x[0], query.fetchall()))

    ##-----------------------------------------------------------------------------------------
    def GetNidsNotInSet( self, reference_table_name, field_name_nid = "nid", condition = None ):
        """retrieve nids in nrdb40, which have not been analysed using radar."""
        statement = """
        SELECT DISTINCTROW a.nid FROM %s AS a
        LEFT JOIN %s AS b
        ON a.nid = b.%s
        WHERE b.%s IS NULL""" % ( self.name, reference_table_name, field_name_nid, field_name_nid)

        if condition:
            statement +=" AND " + condition

        return tuple(map( lambda x: x[0], self.Execute(statement).fetchall()))

    ##-----------------------------------------------------------------------------------------
    def IsContained( self, nid ):
        """returns true, if nid is part of nrdb40."""
        statement = "SELECT nid FROM " + self.name + "WHERE nid = %i" % nid
        result = self.Execute(statement).fetchall()
        if len(result) > 0:
            return 1
        else:
            return 0
        

    ##-----------------------------------------------------------------------------------------------
    def GetSumLength( self ):
        """Get the sum of the length of all sequences."""

        statement = "SELECT SUM(r.length) FROM %s AS n, " % self.name +\
                    " pairsdb.nrdb AS r WHERE r.nid = n.nid " 
        
        return self.Execute(statement).fetchone()[0]

#---------------------------------------------------------------------------------------------------------------    
if __name__ == '__main__':
    
    dbhandle = Pairsdb.Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed"
	sys.exit(1)

    x = TableNids( dbhandle )














