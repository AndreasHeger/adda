## $Id$
##
##
import sys, re, string, os

import Intervalls
from Pairsdb import *
from Experiment import Experiment
from Table_nrdb import Table_nrdb
from TableDomainsCore import TableDomainsCore, TableFamiliesCore
from TableDomains import TableDomains, TableFamilies
from TableNids import TableNids

#-------------------------------------------
# Class:	       Brackets
# Superclasses:  Message
# Subclasses:    
# Function:      update nrdb-database
#
# Author:		Andreas Heger
#-------------------------------------------
class Domains (Experiment):

    def __init__ (self, dbhandle):

        ## source table, needed for adding singletons
        self.mTableNameSource = "nrdb40"

        ## mapping table, needed for mapping domains onto nrdb100
        self.mTableNameMapping = "pairsdb_100x40"
        
        ## destination tables
        self.mTableNameFamilies = None
        self.mTableNameDomains = None

        self.mTableNameMappedDomains = None
        self.mTableNameMappedFamilies = None
        self.mTableNameRepresentatives = None
        
        self.mWorkSpace = "temp"

        self.mDatabase = "pairsdb"
        
        ## whether to remove domains overlapping with repeats
        self.mFilterRepeats = None

        self.mShortOptions  += 'D:f:d:e:g:h:s:r:m:'
        self.mLongOptions   += ['Database=',
                                'families=','domains=',
                                'mapped_families=','mapped_domains=',                                
                                'workspace=', "representatives=",
                                'source=', 'repeats=', "mapping=",
                                'filter_repeats','combine_overlaps',
                                'min_domain_length=']

        self.mMinSingletonLength = 30
        self.mMinDomainLength = 20

        self.mCombineOverlaps = 0

        self.mDbHandle = dbhandle
        
	Experiment.__init__( self )

        self.mDbHandle.UseDatabase( self.mDatabase )
        
        self.mFileNameFamilies = None
        self.mFileNameDomains = None

        # get suffix
        x = string.rfind( self.mTableNameDomains, "_") 
        if x >= 1:
            self.mGeneralSuffix = self.mTableNameDomains[x:] 
        else:
            self.mGeneralSuffix = ""

        if not self.mTableNameDomains:
            raise "no table domains specified"
        if not self.mTableNameFamilies:
            raise "no table families specified"

        if string.find( self.mTableNameFamilies, ".") == -1:
            self.mTableNameFamilies = self.mWorkSpace + "." + self.mTableNameFamilies

        if string.find( self.mTableNameDomains, ".") == -1:
            self.mTableNameDomains = self.mWorkSpace + "." + self.mTableNameDomains

        self.mTableNrdb = Table_nrdb( self.mDbHandle )

        if self.mTableNameSource:
            self.mTableNids = TableNids( self.mDbHandle )
            self.mTableNids.SetName( self.mTableNameSource)
            
        self.mTableFamilies = TableFamilies( dbhandle, "generic" )
        self.mTableDomains = TableDomains( dbhandle, "generic" )

        self.mTableFamilies.SetName( self.mTableNameFamilies )        
        self.mTableDomains.SetName( self.mTableNameDomains  )
        
    ##---------------------------------------------------------------------------------
    def Finish(self):
        return os.path.exists("stop")

    ##---------------------------------------------------------------------------------        
    def ProcessOptions( self, optlist ):

        Experiment.ProcessOptions( self, optlist)

        for o,a in optlist:
            if o in ("-d", "--domains"):
                self.mTableNameDomains = a
            elif o in ("-f", "--families"):
                self.mTableNameFamilies = a
            elif o in ("-e", "--mapped_domains"):
                self.mTableNameMappedDomains = a
            elif o in ("-g", "--mapped_families"):
                self.mTableNameMappedFamilies = a
            elif o in ("-h", "--representatives"):
                self.mTableNameRepresentatives = a
            elif o in ("-w", "--workspace"):
                self.mWorkSpace = a
            elif o in ("-D", "--Database"):
                self.mDatabase = a
            elif o == "--min_domain_length":
                self.mMinDomainLength = string.atoi(a)
            elif o in ("-s", "--source"):
                self.mTableNameSource = a
            elif o in ("-m", "--mapping"):
                self.mTableNameMapping = a
            elif o in ("-r", "--repeats"):
                self.mTableNameDomainsCore = a
                self.mTableDomainsCore = TableDomainsCore( self.mDbHandle, "core" )
                self.mTableDomainsCore.SetName( self.mTableNameDomainsCore )

                self.mTableNameFamiliesCore = re.sub("domains", "families",a)
                self.mTableFamiliesCore = TableFamiliesCore( self.mDbHandle, "core" )
                self.mTableFamiliesCore.SetName( self.mTableNameFamiliesCore )
            elif o == "--filter_repeats":
                self.mFilterRepeats = 1
            elif o == "--combine_overlaps":
                self.mCombineOverlaps = 1
                
    #----------------------------------------------------------------------------------------------------------
    def Load( self ):
        """Load data into table."""

        if self.mFileNameFamilies:
            self.mTableFamilies.Drop()
            self.mTableFamilies.Create()
            self.mTableFamilies.Load( self.mFileNameFamilies )

        if self.mFileNameDomains:
            self.mTableDomains.Drop()
            self.mTableDomains.Create()
            self.mTableDomains.Load( self.mFileNameDomains )

    #----------------------------------------------------------------------------------------------------------
    def OpenOutfiles( self ):
        """Open output files.
        """

        if self.mFileNameFamilies:
            self.mFileFamilies = open(self.mFileNameFamilies, "w")
            
        if self.mFileNameDomains:
            self.mFileDomains    = open(self.mFileNameDomains, "w")

    #----------------------------------------------------------------------------------------------------------
    def CloseOutfiles( self ):
        """Close output files.
        """

        if self.mFileNameFamilies:
            self.mFileFamilies.close()
        if self.mFileNameDomains:
            self.mFileDomains.close()
        
    #----------------------------------------------------------------------------------------------------------
    def UpdateDomains( self ):
        """update domains table."""
        self.mTableFamilies.Update(self.mTableDomains )

    #----------------------------------------------------------------------------------------------------------
    def Renumber( self ):
        """renumber domains starting from 1.
        This procedure tries to keep the low numbered domains conserved
        and to minimize the number of changes.
        Starting from 1 look for empty slots and fill those slots with
        the last domain until no more domains are to be processed.
        """

        families = list(self.mTableFamilies.GetAllFamilies())
        
        new_family = 1
        
        while len(families) > 0:
            
            family = families[0]
            del family[0]

            # enters loop only, if family is empty
            while new_family < family:
                old_family = families[-1]
                del families[-1]
                self.mTableFamilies.RenumberDomain( old_family, new_family)
                self.mTableFamilies.RenumberDomain( old_family, new_family)                
                new_family += 1

            new_family = family + 1
            
    #----------------------------------------------------------------------------------------------------------    
    def Add( self ):
        """Add data into table."""

        if self.mFileNameFamilies:
            self.mTableFamilies.Load( self.mFileNameFamilies )
        if self.mFileNameDomains:
            self.mTableDomains.Load( self.mFileNameDomains )

    #--------------------------------------------------------------------------------
    def Dump( self, path = None, suffix = None):
        """Write a full dump of relevant tables to disk.
        """

        if not path or path == ".":
            path = os.getcwd()

        if not suffix: suffix = ""
        ##---------------------------------------------------------------
        filename = "%s/%s%s.txt" % (path, self.mTableNameDomains, suffix)

        if self.mLogLevel >= 1:
            print "--> dumping table %s into %s" % (self.mTableNameDomains, filename)
            sys.stdout.flush()
            
        self.mTableDomains.Dump( filename )
        
        ##---------------------------------------------------------------
        filename = "%s/%s%s.txt" % (path, self.mTableNameFamilies, suffix)
        if self.mLogLevel >= 1:
            print "--> dumping table %s into %s" % (self.mTableNameFamilies, filename)
            sys.stdout.flush()

        self.mTableFamilies.Dump( filename)

    #--------------------------------------------------------------------------------
    def LoadDump( self, path = None):
        """Write a full dump of relevant tables to disk.
        """

        if not path or path == ".":
            path = os.getcwd()
            
        ##---------------------------------------------------------------
        filename = "%s/%s.txt" % (path, self.mTableNameDomains)
        if self.mLogLevel >= 1:
            print "--> loading %s into table %s" % (filename, self.mTableNameDomains)
            sys.stdout.flush()
            
        self.mTableDomains.LoadDump( filename )
        
        ##---------------------------------------------------------------
        filename = "%s/%s.txt" % (path, self.mTableNameFamilies)
        if self.mLogLevel >= 1:
            print "--> loading %s into table %s" % (filename, self.mTableNameFamilies)
            sys.stdout.flush()

        self.mTableFamilies.LoadDump( filename)

    #--------------------------------------------------------------------------------
    def WriteStatistics( self ):

        print "--> domains: %i, alis: %i " % (self.mTableFamilies.RowCount(), self.mTableDomains.RowCount())
            
        sys.stdout.flush()

    #-------------------------------------------------------------------------------------------------------
    def WriteIntervalls( self, family, nid, intervalls, repeats = None):

        i = intervalls
        if repeats:
            r = []
            for rfamily, rfrom, rto in repeats:
                r.append( (rfrom, rto) )
            i = Intervalls.ShortenIntervallsOverlap( i, r )

        for first_res, last_res in i:
            l = last_res - first_res + 1
            if l >= self.mMinDomainLength:
                self.mFileDomains.write( string.join( map(str, (
                    nid, first_res, last_res, "+%i" % l,
                    nid, first_res, last_res, "+%i" % l, 
                    family)), "\t") + "\n")
                
    #-------------------------------------------------------------------------------------------------------        
    def Finalize( self ):
        """process each sequence and fine tune domain boundaries.
        adds singletons as well.
        """

        nids = self.mTableNids.GetAllNids()        

        if self.mLogLevel >= 1:
            print "--> at the beginning: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(self.mTableDomains.GetAllNids())
            print "--> in %s: %i" % (self.mTableNameSource, len(nids))
            sys.stdout.flush()

        new_family = self.mTableFamilies.GetMaxFamily() + 1

        self.OpenOutfiles()
        nsingletons = 0
        
        for nid in nids:

            if self.mFilterRepeats:
                repeats = self.mTableDomainsCore.GetDomainBoundaries( nid )
            else:
                repeats = None
            
            domains = list(self.mTableDomains.GetDomainBoundaries(nid))
            length = self.mTableNrdb.GetLength( nid )
            
            domains.sort()
            all_intervalls = []
            last_family = None
            for family, domain_from, domain_to in domains:

                if last_family != family:
                    if last_family:
                        if self.mCombineOverlaps:
                            i = Intervalls.CombineIntervallsLarge( family_intervalls )
                        else:
                            i = family_intervalls
                            
                        all_intervalls += i
                        self.WriteIntervalls( last_family, nid, i, repeats)

                    family_intervalls = []
                    
                last_family = family
                family_intervalls.append( (domain_from, domain_to) )

            if last_family:
                if self.mCombineOverlaps:
                    i = Intervalls.CombineIntervallsLarge( family_intervalls )
                else:
                    i = family_intervalls
                    
                all_intervalls += i                
                self.WriteIntervalls( last_family, nid, i, repeats)

            # remove all domains that overlap with repeats by adding the repeats
            if self.mFilterRepeats:

                for rfamily, rfrom, rto in repeats:
                    all_intervalls.append( (rfrom, rto) )
                    
            # add singletons
            i = Intervalls.ComplementIntervalls( all_intervalls, 1, length)

            if self.mLogLevel > 3:
                print "nid=%i" % nid, all_intervalls, repeats, domains, i

            for first_res, last_res in i:
                if last_res-first_res > self.mMinSingletonLength:
                    self.WriteNewSingleton( new_family, nid, first_res, last_res )
                    new_family += 1
                    nsingletons += 1
            
        self.CloseOutfiles()

        self.Load()

        if self.mLogLevel >= 1:
            print "--> at the end: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(self.mTableDomains.GetAllNids())
            print "--> singletons added: %i" % nsingletons
            sys.stdout.flush()
        
    ##-------------------------------------------------------------------------------------        
    def WriteNewSingleton( self, family, nid, first_res, last_res ):
            
        """add temporary files to tables.
        """

        l = last_res - first_res + 1
        
        self.mFileDomains.write( string.join( map(str, (
            nid, first_res, last_res, "+%i" % l,
            nid, first_res, last_res, "+%i" % l,
            family)), "\t") + "\n")
        
        self.mTableFamilies.AddFamily( family )
        
    ##-------------------------------------------------------------------        
    def AddSingletons( self ):
        """for each representative add singleton domains
        for each uncovered segment.
        """

        nids = self.mTableNids.GetAllNids()        

        if self.mLogLevel >= 1:
            print "--> at the beginning: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(self.mTableDomains.GetAllNids())
            print "--> in %s: %i" % (self.mTableNameSource, len(nids))
            sys.stdout.flush()

        family = self.mTableFamilies.GetMaxFamily() + 1

        self.OpenOutfiles()

        for nid in nids:
            domains = self.mTableDomains.GetDomainBoundaries( nid )
            length  = self.mTableNrdb.GetLength( nid )

            first_from = 1
            id = 1
            for (xfamily, domain_from, domain_to) in domains:
                if (domain_from - first_from > self.mMinSingletonLength):

                    self.WriteNewSingleton( family, nid,
                                            first_from,
                                            domain_from - 1)
                    
                    family += 1
                first_from = domain_to + 1


            if (length + 1 - first_from > self.mMinSingletonLength):
                self.WriteNewSingleton( family, nid,
                                        first_from,
                                        length)
                family += 1
                
        self.CloseOutfiles()
        self.Add()

        nids = self.mTableDomains.GetAllNids()        
        if self.mLogLevel >= 1:
            print "--> at the end: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(nids)
            sys.stdout.flush()

    ##---------------------------------------------------------------------------------
    def AddRepeats( self ):
        """add repeats from domain definition to table.
        """

        offset = self.mTableFamilies.GetMaxFamily() + 1
        
        nids = self.mTableNids.GetAllNids()        

        if self.mLogLevel >= 1:
            print "--> at the end: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(self.mTableDomains.GetAllNids())
            print "----> adding at offset %i" % offset            
            sys.stdout.flush()


        self.mTableDomains.AddDomainsFromDomainsTable( self.mTableNameDomainsCore,
                                                       offset=offset,
                                                       subset=self.mTableNameSource)

        families = self.mTableFamilies.GetMissingFamilies( self.mTableDomains )

        if self.mLogLevel >= 1:
            print "--> adding %i domains" % len(families)
            sys.stdout.flush()

        for family in families:
            self.mTableFamilies.AddFamilyFromFamiliesTable( self.mTableFamiliesCore,
                                                            family-offset,
                                                            family)
            
        if self.mLogLevel >= 1:
            print "--> at the end: "
            print "--> domains: %i" % (self.mTableFamilies.RowCount())
            print "--> alignments: %i" % (self.mTableDomains.RowCount())
            print "--> nids: %i" % len(self.mTableDomains.GetAllNids())
            sys.stdout.flush()

    ##---------------------------------------------------------------------------------
    def PropagateDown( self ):
        """propagate domains to level given by TableMapping."""

        if self.mLogLevel >= 1:
            print "--> mapping information."
            sys.stdout.flush()

        if not self.mTableNameMappedDomains:
            raise "please specify table of mapped domains."

        if not self.mTableNameMappedFamilies:
            raise "please specify table of mapped families."

        if not self.mTableNameRepresentatives:
            raise "please specify table of representatives."

        if not self.mTableNameMapping:
            raise "please specify table of alignments for mapping."
        
        result_table_domains  = self.mTableNameMappedDomains
        result_table_families = self.mTableNameMappedFamilies        
        
        if self.mLogLevel >= 1:
            print "--> mapping information from %s to %s" % (self.mTableNameDomains, result_table_domains)
            print "--> using alignments in %s" % (self.mTableNameMapping)
            print "--> using representatives in %s" % (self.mTableNameRepresentatives)
            sys.stdout.flush()

        new_tbl_domains = self.mTableDomains.GetAnotherInstance()
        new_tbl_domains.SetName( result_table_domains )
        new_tbl_domains.Drop()
        new_tbl_domains.Create()
        new_tbl_domains.mLogLevel = self.mLogLevel
        new_tbl_domains.Fill( self.mTableNameRepresentatives, self.mTableNameMapping, self.mTableNameDomains )

        if self.mLogLevel >= 1:
            print "--> entries in %s: %i" % (new_tbl_domains.GetName(),
                                             new_tbl_domains.RowCount())
            sys.stdout.flush()

        ##--------------------------------------------------------------------
        ## map annotation

        if self.mLogLevel >= 1:
            print "--> mapping information from %s to %s." % (self.mTableNameFamilies, result_table_families)
            sys.stdout.flush()

        new_tbl_families = self.mTableFamilies.GetAnotherInstance()
        new_tbl_families.SetName( result_table_families )
        new_tbl_families.Drop()
        new_tbl_families.Create()
        new_tbl_families.Fill( self.mTableFamilies.GetName(), new_tbl_domains )
        new_tbl_families.Update( new_tbl_domains )

        if self.mLogLevel >= 1:
            print "--> entries in %s: %i" % (new_tbl_families.GetName(),
                                             new_tbl_families.RowCount())
            sys.stdout.flush()
        
        return

    ##---------------------------------------------------------------------------------
    def PropagateUp( self ):
        """propagate domains to level given by TableMapping."""
        
        if self.mLogLevel >= 1:
            print "--> mapping information."
            sys.stdout.flush()

        src_table_domains = self.mTableDomains.GetName()

        if not self.mTableNameMappedDomains:
            prefix, suffix = self.mTableDomains.GetName().split(".")
            result_table_domains = prefix + "." + re.sub("\d+","",suffix,1)
        else:
            result_table_domains = self.mTableNameMappedDomains
        
        if self.mLogLevel >= 1:
            print "--> mapping information from %s to %s" % (self.mTableNameDomains, result_table_domains)
            print "--> using alignments in %s" % (self.mTableNameMapping)
            sys.stdout.flush()

            
        new_tbl_domains = self.mTableDomains.GetAnotherInstance()
        new_tbl_domains.SetName( result_table_domains )
        new_tbl_domains.Drop()
        new_tbl_domains.Create()
        
        new_tbl_domains.mLogLevel = self.mLogLevel
        
        new_tbl_domains.ReverseFill( self.mTableNameMapping, self.mTableNameDomains )

        if self.mLogLevel >= 1:
            print "--> entries in %s: %i" % (new_tbl_domains.GetName(),
                                             new_tbl_domains.RowCount())
            sys.stdout.flush()

        ##--------------------------------------------------------------------
        ## map annotation
        if not self.mTableNameMappedFamilies:
            prefix, suffix = self.mTableFamilies.GetName().split(".")
            result_table_families = prefix + "." + re.sub("\d+","",suffix,1)
        else:
            result_table_families = self.mTableNameMappedFamilies

        if self.mLogLevel >= 1:
            print "--> mapping information from %s to %s." % (self.mTableNameFamilies, result_table_families)
            sys.stdout.flush()

        new_tbl_families = self.mTableFamilies.GetAnotherInstance()
        new_tbl_families.SetName( result_table_families )
        new_tbl_families.Drop()
        new_tbl_families.Create()
        new_tbl_families.Fill( self.mTableFamilies.GetName(), new_tbl_domains )
        new_tbl_families.Update( new_tbl_domains )
        
        return

    ##---------------------------------------------------------------------------------
    def CreateTables( self ):
        """select all members of a connected component and
        add to families. Each connected component gives the
        component_id for a given family.
        
        """
        self.mTableFamilies.Drop()
        self.mTableDomains.Drop()
        self.mTableFamilies.Create()
        self.mTableDomains.Create()

#--------------------------------------< end of class definition >----------------------------------

if __name__ == '__main__':

    dbhandle = Pairsdb()
    if not dbhandle.Connect():
	print "Connection failed"
	sys.exit(1)

    x = Domains( dbhandle )

    x.Process()
    

    

