import sys, os, re, time

import alignlib
import ProfileLibrary

from AddaModule import AddaModule
import AddaIO

class AddaProfiles( AddaModule ):
    """write a links table for adda processing."""
    
    mName = "Profiles"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )
                
        self.mProfileLibrary = ProfileLibrary.ProfileLibrary( self.mConfig.get( "files", "output_profiles"), 
                                                              "w",
                                                              force=self.mOptions.force )
        
        self.mScaleFactor  = self.mConfig.get( "profiles", "scale_factor", 0.3 )
        self.mMaxNumNeighbours = self.mConfig.get( "profiles", "max_neighbours", 1000)
        self.mPrepareProfile = self.mConfig.get( "profiles", "prepare_profile", False ) 
        
        self.mLogOddor    = alignlib.makeLogOddorDirichlet( self.mScaleFactor )
        self.mRegularizor = alignlib.makeRegularizorDirichletPrecomputed()
        self.mWeightor    = alignlib.makeWeightor()
        
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
            
        if self.mOptions.loglevel >= 6:
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

        mali = self.buildMali( neighbours )

        query_nid = neighbours.mQueryToken

        self.debug( "working on profile %s" % query_nid )
            
        profile = alignlib.makeProfile( mali,
                                        alignlib.getDefaultEncoder(),
                                        self.mWeightor,
                                        self.mRegularizor, 
                                        self.mLogOddor )
        
        profile.setStorageType( alignlib.Sparse )
        if self.mPrepareProfile: profile.prepare()

        self.mProfileLibrary.add( query_nid, profile )
            
    def finish( self ):
        """finish processing.
        
        add entries for sequences who only appear in the sbjct field.
        """
        nids = self.mFasta.getContigSizes().keys()
        
        nadded = 0

        for nid in sorted(nids):
            if nid not in self.mProfileLibrary:
                self.applyMethod( AddaIO.NeighboursRecord( nid, [] ) )
                nadded += 1
                
        self.mOutput += nadded
        self.info( "added %i profiles for sequences without neighbours" % nadded )
        
        AddaModule.finish(self)        
