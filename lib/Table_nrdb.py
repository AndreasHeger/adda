####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Table_nrdb.py,v 1.2 2002/11/18 13:03:28 heger Exp $
##
##
####
####


#----------------------------------------------------------------
# Name:		Table_nrdb
#--------------General Information-------------------------------
# File:		Table_nrdb.py
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

class Table_nrdb( Table ):
    def __init__ ( self, handle ):
	self.name   =  'nrdb'
	self.fields = (
			('updated', 'TIMESTAMP'),
			('nid', 'INTEGER UNSIGNED NOT NULL AUTO_INCREMENT'),
			('hid', 'CHAR(22) NOT NULL'),
			('sequence', 'TEXT NOT NULL'),
			('sequencedbs_id', 'TINYINT UNSIGNED NOT NULL'),
			('accessionnumber', 'VARCHAR(20) BINARY NOT NULL'),
			('identifier', "VARCHAR(20) BINARY  NOT NULL DEFAULT ''"),
			('description', "VARCHAR(255) NOT NULL DEFAULT ''"),
			('created', 'DATE NOT NULL'),
			('length', 'MEDIUMINT UNSIGNED NOT NULL DEFAULT 0'),
			('filter', 'TINYINT   UNSIGNED NOT NULL DEFAULT 100'),
		)
	self.indices = (
			'PRIMARY KEY (nid)',
			'INDEX  hid (hid(4))',
			'INDEX  updated (updated)',
			'INDEX  filter (filter)',
		)
        Table.__init__( self, handle )

    ##-------------------------------------------------------------
    def Check( self ):
        """Check if data in table is consistent"""
        # data in other tables will be used for the check
        # select sequences which are not a nrdb90 representative nor a nrdb90-groupie
        statement = "SELECT nrdb.nid, nrdb.identifier, nrdb.description, nrdb.length FROM nrdb " +\
                    "LEFT JOIN nrdb90         ON nrdb.nid = nrdb90.nid " + \
                    "LEFT JOIN pairsdb_100x90 ON nrdb.nid = pairsdb_100x90.mem_nid " +\
                    "WHERE nrdb90.nid is NULL AND pairsdb_100x90.mem_nid is NULL " +\
                    "AND nrdb.length >= 10"
        Table.Check( self, statement, "sequences neither nrdb90-representative nor groupie" )
        
              
    ##-------------------------------------------------------------
    def CheckSequence ( self, hid ):
	statement = "SELECT hid FROM " + self.name + " WHERE hid = '" + hid + "'"
	c = self.Execute( statement )
	if c.rowcount == 1:
	    return 1
	else:
	    return 0
	

    ##-------------------------------------------------------------
    def Get_NID_From_HID( self, hid ):
        """Return the nid correspondig to an hid"""
        query = self.Execute( "SELECT nid FROM " + self.name + " WHERE hid = '" + hid + "'")
        try:
            x = query.fetchone()[0]
        except:
            x = None
        return x

    ##-------------------------------------------------------------
    def SetIdentification( self, nid, dbs_id, identifier, accession, description ):
        """set preferred id."""
        statement = """
        UPDATE %s SET sequencedbs_id = %i, identifier='%s', accessionnumber='%s', description='%s'
        WHERE nid = %i""" %(self.name, dbs_id, identifier, accession,
                            self.dbhandle.QuoteString(description),
                            nid)
        return self.Execute( statement )

    ##-------------------------------------------------------------
    def GetIdentifier( self, nid ):
        """return identifier for nid"""
        statement = "SELECT identifier FROM %s  WHERE nid = %i " % (self.name, nid)
        return self.Execute( statement ).fetchone()[0]

    ##-------------------------------------------------------------
    def GetAccessionNumber( self, nid ):
        """return accessionnumber for nid"""
        statement = "SELECT accessionnumber FROM %s  WHERE nid = %i " % (self.name, nid)
        return self.Execute( statement ).fetchone()[0]

    ##-------------------------------------------------------------
    def GetDescription( self, nid ):
        """return description for nid"""
        statement = "SELECT description FROM %s WHERE nid = %i " % (self.name, nid)
        return self.Execute( statement ).fetchone()[0]

    ##-------------------------------------------------------------
    def Get_Last_NID ( self ):
        """Return the last nid used in database"""
        query = self.Execute( "SELECT LAST_INSERT_ID()")
        x = query.fetchone()[0]
        if x:
            return int(x)
        else:
            return 0
    ##-------------------------------------------------------------        
    def GetMaxNid( self ):
        """return maximum nid."""
        return self.Execute( "SELECT MAX(nid) FROM %s" % self.name).fetchone()[0]
    
    ##-------------------------------------------------------------
    def Get_Sequence_From_NID( self, nid ):
        """Return sequence for nid"""
        query = self.Execute( "SELECT sequence FROM " + self.name + " WHERE nid = '" + str(nid) + "'")
        try:
            x = query.fetchone()[0]
        except:
            x = None
        return x

    ##-------------------------------------------------------------
    def GetSequence( self, nid ):
        """return sequence.
        """
        return self.Execute("SELECT sequence FROM " + self.name + " WHERE nid = %i " % nid ).fetchone()[0]

    ##-------------------------------------------------------------
    def GetObsoleteSequences( self ):
        """return all sequences.
        """
        return self.Execute("SELECT nid, sequence FROM " + self.name + " WHERE filter = 0").fetchall()

    ##-------------------------------------------------------------
    def GetAllSequences( self ):
        """return all sequences.
        """
        return self.Execute("SELECT sequence FROM " + self.name + " WHERE filter > 0").fetchall()
    
    ##-------------------------------------------------------------
    def GetLength( self, nid ):
        """return sequence length."""
        result = self.Execute( "SELECT length FROM %s WHERE nid = '%i'" % (self.name, nid)).fetchone()

        if result:
            return result[0]
        else:
            return 0

    ##-------------------------------------------------------------
    def GetNumSequences( self ):
        """return numbers of sequences."""
        result = self.Execute( "SELECT COUNT(*) FROM %s WHERE filter != 0" % self.name).fetchone()

        if result:
            return result[0]
        else:
            return 0
    ##-------------------------------------------------------------
    def ResetFilter( self ):
        """Reset filter, i.e. set all values to 100"""
        statement = "UPDATE %s SET filter = 100 WHERE filter != 0" % self.name
        return self.Execute( statement )

    ##-------------------------------------------------------------
    def AddFilter( self, level ):
        """set filter according to nrdb% level database."""

        if self.dbhandle.Exists( "nrdb%i" % level):
            nids = map( lambda x: x[0], self.Execute("SELECT nid FROM nrdb%i" % level).fetchall())
            for nid in nids:
                self.Execute("UPDATE %s SET filter=%i WHERE nid = %i" % (self.name, level, nid))
                
    ##-------------------------------------------------------------
    def Flag_Old_Entries( self ):
        """Flags old entries as such by setting their filter to 0. Old entries are identified
        by looking at cross_references and looking up those entries from nrdb, which have no entry
        in cross_refernces."""
        statement = "SELECT DISTINCT n.nid FROM nrdb AS n LEFT JOIN cross_references AS c ON n.nid = c.nid "+\
                    " WHERE c.nid IS NULL"
        nids = map(lambda x: x[0], self.Execute( statement ).fetchall())

        for nid in nids:
            self.Execute( "UPDATE nrdb SET filter = 0 WHERE nrdb.nid = " + str(nid) )

        return len(nids)
    
    ##-------------------------------------------------------------
    def TSelect_Filtered_Sorted_By_Length( self, cutoff ):
        """Return list of nids and lengths which are below or equal to cutoff, sort descending by length"""
        query = self.Execute( "SELECT nid,length FROM nrdb WHERE filter > 0 AND filter <= " + str(cutoff) + " ORDER BY length DESC, nid ASC ")  
        return query.fetchall()

    ##-------------------------------------------------------------
    def TSelect_Filtered( self, cutoff ):
        """Return list of nids which are below or equal to cutoff"""
        query = self.Execute( "SELECT nid FROM nrdb WHERE filter > 0 AND filter <= " + str(cutoff) )  
        return query.fetchall()

    ##-------------------------------------------------------------
    def Filter_Update_Single( self, nid, cutoff):
        """Update a single sequence with filter-cutoff"""
        statement = "UPDATE nrdb SET filter = " + str(cutoff) + " WHERE nid = " + str( nid )
        return self.Execute( statement )

    ##-------------------------------------------------------------
    def GetAllNids( self ):
        statement = "SELECT nid FROM " +self.name + " WHERE filter > 0 "
        return map( lambda x: x[0], self.Execute( statement ).fetchall())

    ##-------------------------------------------------------------
    def GetAllActiveEntries( self ):
        statement = "SELECT nid, sequence FROM " + self.name + " WHERE filter > 0 "        
        return self.Execute(statement)
    
    ##--------------------------------------------------------------
    def GetAnnotationFromNid( self, nid ):
        return self.Execute("SELECT identifier, description FROM " + self.name +\
                            " WHERE nid = %i " % nid ).fetchall()[0]

    ##--------------------------------------------------------------
    def GetFragments( self ):
        statement = "SELECT nid FROM " + self.name + " WHERE description LIKE '%fragment%'"
        return map( lambda x: x[0], self.Execute(statement).fetchall())

    ##--------------------------------------------------------------
    def GetHypotheticalProteins( self ):
        statement = "SELECT nid FROM " + self.name + " WHERE description LIKE '%hypothetical%'"
        return map( lambda x: x[0], self.Execute(statement).fetchall())

    ##--------------------------------------------------------------
    def GetArtificialProteins( self ):
        statement = "SELECT nid FROM " + self.name + " WHERE description LIKE '%artificial%'"
        return map( lambda x: x[0], self.Execute(statement).fetchall())


    ##--------------------------------------------------------------
    # note: just synthetic also picks up biosynthetic and photosynthetic
    def GetSyntheticProteins( self ):    
        statement = "SELECT nid FROM " + self.name +\
                    " WHERE description LIKE '%synthetic%' " +\
                    " AND description NOT LIKE '%biosynthetic%' " +\
                    " AND description NOT LIKE '%photosynthetic%' "
        return map( lambda x: x[0], self.Execute(statement).fetchall())        

    ##--------------------------------------------------------------
    # note: test for swall and no _ in identifier
    def GetNonCuratedProteins( self ):    
        statement = "SELECT nid FROM " + self.name +\
                    " WHERE sequencedbs_id <> 1 " +\
                    " OR identifier NOT LIKE '%\_%'"
        return map( lambda x: x[0], self.Execute(statement).fetchall())        
        
    ##--------------------------------------------------------------
    def GetDatabase( self, dbs_id ):
        statement = "SELECT nid FROM " + self.name + " WHERE sequencedbs_id = %i " % dbs_id
        return map( lambda x: x[0], self.Execute(statement).fetchall())

    ##--------------------------------------------------------------
    def GetFilteredSequences(self, max_filter ):
        statement = "SELECT nid, sequence FROM " + self.name +\
                    " WHERE filter > 0 AND filter <= %i LIMIT 100" % max_filter
        
        
        return self.Execute(statement).fetchall()

    ##--------------------------------------------------------------
    def GetTotalLength(self, max_filter, subset = None ):
        """retrieve the total length of all sequences <= filter.
        """
        if subset:
            s = "INNER JOIN %s AS s ON s.nid = n.nid" % (subset)
        else:
            s = ""
            
        statement = "SELECT SUM(length) FROM %s AS n" % self.name +\
                    " %s WHERE n.filter > 0 AND n.filter <= %i " % (s,max_filter)
        
        return self.Execute(statement).fetchone()[0]

    ##--------------------------------------------------------------
    def GetTotalLengthNids(self, max_filter, subset = None ):
        """retrieve the total length of all sequences <= filter.
        """
        if subset:
            s = "INNER JOIN %s AS s ON s.nid = n.nid" % (subset)
        else:
            s = ""
        
        statement = "SELECT COUNT(DISTINCT n.nid), SUM(n.length) FROM %s AS n " % self.name +\
                    " %s WHERE n.filter > 0 AND n.filter <= %i " % (s, max_filter)
        return self.Execute(statement).fetchone()

#---------------------------------------------------------------------------------------------------------------    
#---------------------------------------------------------------------------------------------------------------    
if __name__ == '__main__':
    
    dbhandle = Pairsdb.Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed"
	sys.exit(1)

    x = Table_nrdb( dbhandle )
    ## print x.Flag_Old_Entries()
    x.Print_Statistics()

