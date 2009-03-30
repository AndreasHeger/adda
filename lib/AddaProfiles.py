import sys, os, re, time, glob

import alignlib
import ProfileLibrary

from AddaModule import AddaModule
import AddaIO
import SegmentedFile

import Experiment as E

class AddaProfiles( AddaModule ):
    """write a links table for adda processing."""
    
    mName = "Profiles"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )

        self.mFilenameProfile = self.mConfig.get( "files", "output_profiles", "adda.profiles" )
        self.mScaleFactor  = self.mConfig.get( "profiles", "scale_factor", 0.3 )
        self.mMaxNumNeighbours = self.mConfig.get( "profiles", "max_neighbours", 1000)
        self.mPrepareProfile = self.mConfig.get( "profiles", "prepare_profile", False ) 

    def isComplete( self ):

        fn, fi = ProfileLibrary.getFileNames( self.mFilenameProfile + self.getSlice() )
        return SegmentedFile.isComplete( fi )

    def startUp( self ):

        if self.isComplete(): return

        if self.mAppend:
            self.mProfileLibrary = ProfileLibrary.ProfileLibrary( self.mFilenameProfile + self.getSlice(),
                                                                  "a" )
            self.mContinueAt = self.mProfileLibrary.getLastInsertedKey()
            self.info("processing will continue after %s" % (str( self.mContinueAt ) ) )
        else:
            self.mProfileLibrary = ProfileLibrary.ProfileLibrary( self.mFilenameProfile + self.getSlice(),
                                                                  "w",
                                                                  force=self.mForce )


        # set default values
        self.mProfileLibrary.setLogOddor( alignlib.makeLogOddorDirichlet( self.mScaleFactor ) )
        self.mProfileLibrary.setRegularizor( alignlib.makeRegularizorDirichletPrecomputed() )
        self.mProfileLibrary.setWeightor( alignlib.makeWeightor() )
        alignlib.setDefaultEncoder( alignlib.getEncoder( alignlib.Protein20 ) )

    def buildMali(self, neighbours):
        """build a multiple alignment from a set of neighbours.
        """
        # build multiple alignment
        mali = alignlib.makeMultipleAlignment()
        
        query_nid = neighbours.mQueryToken
        
        sequence = self.mFasta.getSequence( query_nid )

        mali.add( alignlib.makeAlignatum( sequence ) )

        for n in neighbours.mMatches[:self.mMaxNumNeighbours]:
            if n.mSbjctToken == query_nid: continue

            map_query2sbjct = n.getAlignment()

            if map_query2sbjct.getLength() == 0:
                self.warn( "empty alignment: %s" % str( n ) )
                continue
            
            sequence = self.mFasta.getSequence( n.mSbjctToken )
            mali.add( alignlib.makeAlignatum( sequence ),
                      map_query2sbjct,
                      mali_is_in_row = True, 
                      insert_gaps_mali = False,
                      insert_gaps_alignatum = True,
                      use_end_mali = True,
                      use_end_alignatum = False )
            
        if E.getLogLevel() >= 6:
            x = 1
            outfile = open( "mali_%s" % query_nid, "w" )
            for line in str(mali).split("\n"):
                try:
                    a,b,c = line.split("\t")
                except ValueError:
                    continue
                outfile.write( ">%06i\n%s\n" % (x,b) )
                x += 1
            outfile.close()
    
        return mali
        
    def applyMethod(self, neighbours ):
        """output the graph.

        If mMergeRepeats is set, consecutive links are merged.
        Links are consecutive if they are adjacent both in the
        query and in the sbjct.
        
        This ensures that 1:many repeats are not merged, but will
        cover alignments split by transmembrane regions.
        
        """
        if self.mContinueAt:
            if neighbours.mQueryToken == self.mContinueAt:
                self.info("continuing processing at %s" % str(self.mContinueAt ) )
                self.mContinueAt = None
            return

        mali = self.buildMali( neighbours )

        query_nid = neighbours.mQueryToken

        self.debug( "working on profile %s" % query_nid )
            
        profile = alignlib.makeProfile( mali )
        profile.setStorageType( alignlib.Sparse )
        if self.mPrepareProfile: profile.prepare()

        self.mProfileLibrary.add( query_nid, profile )
            
    def finish( self ):
        """finish processing.
        
        add entries for sequences who only appear in the sbjct field.
        """
        if not self.isSubset():
            nids = self.mFasta.getContigSizes().keys()
            nadded = 0

            for nid in sorted(nids):
                if nid not in self.mProfileLibrary:
                    self.applyMethod( AddaIO.NeighboursRecord( nid, [] ) )
                    nadded += 1
                
            self.mOutput += nadded
            self.info( "added %i profiles for sequences without neighbours" % nadded )

        self.mProfileLibrary.close()
        
        AddaModule.finish(self)        

    #--------------------------------------------------------------------------
    def merge(self):
        """merge runs from parallel computations.

        returns true if merging was succecss.
        """
        if self.isComplete(): return
        
        infiles = glob.glob( "%s*" % self.mFilenameProfile )
        # remove suffixes
        infiles = list(set([ x[:-4] for x in infiles if x != self.mFilenameProfile ]))
        infiles.sort()

        last_nid = None
        found = set()
        ninput, noutput, nfound, nunknown, nduplicate = 0, 0, 0, 0, 0
        tokens = set(self.mFasta.keys())

        self.mProfileLibrary = ProfileLibrary.ProfileLibrary( self.mFilenameProfile,
                                                              "w" )

        for filename in infiles:
            infile = ProfileLibrary.ProfileLibrary( filename, "r" )

            for nid, profile in infile.iteritems_sorted():
                ninput += 1
                
                if nid in found:
                    nduplicates += 1
                    self.warn("duplicate nid: %i in file %s" % (nid, filename))
                if nid not in tokens:
                    nunknown += 1
                    self.warn("unknown nid: %i in file %s" % (nid, filename))
                found.add(nid)
                nfound += 1
                self.mProfileLibrary.add( nid, profile )
                noutput += 1

        missing = tokens.difference( found ) 
        if len(missing) > 0:
            self.warn( "the following nids were missing: %s" % str(missing) )
            
        self.info( "adding %i missing nids" % len(missing))
        
        for nid in missing:
            self.applyMethod( AddaIO.NeighboursRecord( nid, [] ) )

        self.info( "merging: parts=%i, ninput=%i, noutput=%i, nfound=%i, nmissing=%i, nduplicate=%i, nunknown=%i" %\
                       (len(infiles), ninput, noutput, nfound, len(missing), nduplicate, nunknown ) )

        self.info( "deleting %i parts" % len(infiles) )
        for infile in infiles:
            fn, fi = ProfileLibrary.getFileNames( infile )
            os.remove( fn )
            os.remove( fi )
        
        return len(missing) == 0 and nduplicate == 0 and nunknown == 0
