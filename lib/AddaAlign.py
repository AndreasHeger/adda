import sys, os, re, time, math, copy
import alignlib
import ProfileLibrary
from AddaModule import AddaModuleRecord
import SegmentedFile

class AddaAlign( AddaModuleRecord ):
    """align domains."""
    
    mName = "Align"
    
    def __init__(self, *args, **kwargs ):

        AddaModuleRecord.__init__( self, *args, **kwargs )

        self.mFilenameAlignments = self.mConfig.get("files","output_align", "adda.align" )

        self.mFilenameProfiles = self.mConfig.get( "files", "output_profiles", "adda.profiles")
        self.mFilenameMst = self.mConfig.get( "files", "output_mst", "adda.mst" )              
        
        self.mScaleFactor  = self.mConfig.get( "profiles", "scale_factor", 0.3)
        self.mMaxNumNeighbours = self.mConfig.get( "profiles", "max_neighbours", 1000 )
        self.mPrepareProfile = self.mConfig.get( "profiles", "prepare_profile", False ) 
        
        self.mMinOverlapResidues = self.mConfig.get( "align", "min_overlap_residues", 20 )
        self.mMinCoverage = self.mConfig.get( "align", "min_coverage", 0.2 )
        self.mMinOverlap = self.mConfig.get( "align", "min_overlap", 0.2 )
        self.mMask = self.mConfig.get( "align", "mask", True )
        self.mMethodsMask = map(int, self.mConfig.get( "align", "masks", "3,4" ).split(","))
        
        self.mUseCache = self.mConfig.get( "align", "use_cache", True )
        self.mCacheSize = self.mConfig.get( "align", "cache_size", 100 ) 

        ###############################################
        # options for zscore check
        self.mMinZScore = self.mConfig.get( "align", "min_zscore", 5.0 )
        self.mNumIterationsZScore = self.mConfig.get( "align", "num_iterations_zscore", 50 )

        # if score is 5 times the minimum score, do not compute zscore
        self.mSafetyThreshold = self.mConfig.get( "align", "safety_threshold", 5 )

        ###############################################
        # alignment parameters
        self.mGop = self.mConfig.get( "align", "gop", -10.0 )
        self.mGep = self.mConfig.get( "align", "gep", -1.0 )
        
        # minimum size for using a profile for alignments
        self.mMinProfileSize = self.mConfig.get( "align", "min_profile_size", 0 )

        # threshold parameters for significance check
        self.mMinAlignmentScore  = self.mConfig.get( "align", "min_alignment_score", 83.0 )
        self.mMinAlignmentMotifLength = self.mConfig.get( "align", "min_motif_length", 10 )

        self.mFilenames = (self.mFilenameAlignments, )
                
    #--------------------------------------------------------------------------------
    def startUp( self ):

        if self.isComplete(): return

        ###############################################
        # create objects for algorithm 
        self.mProfileLibrary = ProfileLibrary.ProfileLibrary( self.mFilenameProfiles, "r" )
        
        self.mChecker = self.checkLinkZScore
        self.mHeader = ("passed",
                        "qdomain",
                        "sdomain",
                        "qstart",
                        "qend",
                        "qali",
                        "sstart",
                        "send",
                        "sali",
                        "score",
                        "naligned",
                        "ngaps",
                        "zscore" )

        self.mAlignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, 
                                                        self.mGop,
                                                        self.mGep )

        # the cache to store alignandum objects
        self.mCache = {}        
        
        self.mLogOddor    = alignlib.makeLogOddorDirichlet( self.mScaleFactor )
        self.mRegularizor = alignlib.makeRegularizorDirichletPrecomputed()
        self.mWeightor    = alignlib.makeWeightor()
        self.mProfileLibrary.setWeightor( self.mWeightor )
        self.mProfileLibrary.setLogOddor( self.mLogOddor )
        self.mProfileLibrary.setRegularizor( self.mRegularizor )
        
        alignlib.setDefaultEncoder( alignlib.getEncoder( alignlib.Protein20 ) )

        ## initialize counters
        self.mNPassed, self.mNFailed, self.mNNotFound = 0, 0, 0

        self.mOutfile = self.openOutputStream( self.mFilenameAlignments )

        if self.mContinueAt == None:
            self.mOutfile.write( "\t".join( self.mHeader ) + "\n" ) 
            self.mOutfile.flush()

    #--------------------------------------------------------------------------------
    def mask( self, nid, alignandum):
        """mask a sequence or profile with nid.
        (do not mask membrane regions)
        """

        masks = self.mTableMasks.GetMasks( nid, self.mMethodsMask )
        for first_res, last_res, info, method in masks:
            for x in range(first_res, last_res+1):
                    alignandum.MaskResidue( x )

        # mask bias residue wise (otherwise to restrictive)
        masks = self.mTableMasks.GetMasks( nid, (1,) )
        if masks:
            sequence = self.mTableNrdb.GetSequence(nid)
            
            for first_res, last_res, info, method in masks:        
                for s in range(first_res, last_res+1):
                    if sequence[s-1] == info:
                        alignandum.MaskResidue( s )
                        
    #--------------------------------------------------------------------------------
    def getAlignandum( self, nid ):
        """get the alignandum object for an nid."""

        if self.mCache:
            if nid not in self.mCache:
                a = self.mProfileLibrary.getProfile(nid)
                self.mCache[nid] = a
                a.prepare()
                if self.mMask: self.mask( nid, a)
            else:
                a = self.mCache[nid]
        else:
            try:
                a = self.mProfileLibrary.getProfile(nid)
            except KeyError:
                self.warn( "profile for sequence %s not found." % str(nid))
                return None
            
            a.prepare()
            if self.mMask: self.mask( nid, a)

        if self.mLogLevel >= 5:
            E.debug( "alignandum for rep %s\n%s" % ( nid, str(a) ) )

        return a
    
    def registerExistingOutput(self, filename):
        """process existing output in filename to guess correct point to continue computation."""
        last_line = self.getLastLine( filename )
        if last_line:
            code, query_token, sbjct_token = last_line.split( "\t" )[:3]
            self.mContinueAt = (query_token, sbjct_token)
            self.info("processing will continue after pair %s" % str( self.mContinueAt) )    

    def applyMethod( self, line ):
        """output the graph."""

        t_start = time.time()
            
        self.mInput += 1
              
        if self.mContinueAt and (query_token,sbjct_token) == self.mContinueAt:
            self.mContinueAt = None
            t_start = time.time() 
            self.info("continuing processing after iteration %i" % self.mInput )
            return

        try:
            (query_token, sbjct_token) = line[:-1].split("\t")[:2]
        except ValueError:
            raise ValueError( "parsing error in line %s" % line )

        if self.mInput % self.mReportStep == 0:
            t = time.time() 
            self.info( "iteration=%i, passed=%i, failed=%i, notfound=%i, total time=%i, time per step=%f" %\
                           (self.mInput, self.mNPassed, self.mNFailed, self.mNNotFound,
                            t - t_start,
                            float(self.mReportStep * ( t - t_start )) / self.mInput, 
                            ) )

        query_nid, query_from, query_to = query_token.split("_")
        sbjct_nid, sbjct_from, sbjct_to = sbjct_token.split("_")
        query_from, query_to = map(int, (query_from, query_to) )
        sbjct_from, sbjct_to = map(int, (sbjct_from, sbjct_to) )

        self.debug( "checking link between %s (%i-%i) and %s (%i-%i)" %\
                    (query_nid, query_from, query_to,
                     sbjct_nid, sbjct_from, sbjct_to) )

        passed, alignment, extra_info = self.mChecker( query_nid, query_from, query_to,
                                                       sbjct_nid, sbjct_from, sbjct_to)

        if passed: 
            code = "+"
            self.mNPassed += 1
        else:
            code = "-"
            self.mNFailed += 1

        self.mOutfile.write( "\t".join( ( code,
                                    line[:-1],
                                    str(alignlib.AlignmentFormatEmissions( alignment )),
                                    str(alignment.getScore()), 
                                    str(alignment.getNumAligned()), 
                                    str(alignment.getNumGaps())) + extra_info ) + "\n" )                    
        self.mOutfile.flush()

        self.mOutput += 1
           
    def finish(self):    
        
        self.mOutfile.close()
        
        self.info( "aligned: %i links input, %i links passed, %i links failed, %i links not found" %\
                       (self.mInput, self.mNPassed, self.mNFailed, self.mNNotFound ) )
        
        AddaModuleRecord.finish( self )
        
    def checkLinkThreshold( self,
                   query_nid, query_from, query_to,
                   sbjct_nid, sbjct_from, sbjct_to):
        """check, whether two domains are homologous.
        
        The check is done whether the alignment store between the two
        domains is above a score threshold.
        """

        query_profile = self.getAlignandum( query_nid )
        query_profile.useSegment( query_from, query_to )

        sbjct_profile = self.getAlignandum( sbjct_nid )
        sbjct_profile.useSegment( sbjct_from, sbjct_to )        
        
        result = alignlib.makeAlignmentVector()

        alignator.align( result, query_profile, sbjct_profile )

        self.debug( "--> %i vs %i: score=%5.2f, length=%i, numgaps=%i, row_from=%i, row_to=%i, col_from=%i, col_to=%i" %\
                  (query_nid, sbjct_nid,
                   result.getScore(),
                   result.getLength(),
                   result.getNumGaps(),
                   result.getRowFrom(), result.getRowTo(),
                   result.getColFrom(), result.getColTo()) )

        query_profile.useFullLength()
        sbjct_profile.useFullLength()
        
        if result.getScore() > self.mMinAlignmentScore:
            return True,result, ()
        else:
            return False,result, ()

        
    def checkLinkZScore( self,
                         query_nid, query_from, query_to,
                         sbjct_nid, sbjct_from, sbjct_to):
        """check, whether two domains are homologous.
        
        The check is done using a zscore calculation.
        """

        result = alignlib.makeAlignmentVector()
        
        query_profile = self.getAlignandum( query_nid )
        sbjct_profile = self.getAlignandum( sbjct_nid )

        if not query_profile or not sbjct_profile:
            self.warn( "could not compute link %s_%i_%i - %s_%i_%i\n" % \
                       (query_nid, query_from, query_to,
                        sbjct_nid, sbjct_from, sbjct_to) )
            self.mNNotFound += 1
            return False, result, ("na",)
        
        query_profile.useSegment( query_from, query_to )
        sbjct_profile.useSegment( sbjct_from, sbjct_to )        
        
        self.mAlignator.align( result, query_profile, sbjct_profile )
        
        self.debug( "# --> %s vs %s: score=%5.2f, length=%i, numgaps=%i, row_from=%i, row_to=%i, col_from=%i, col_to=%i" %\
                    (query_nid, sbjct_nid,
                     result.getScore(),
                     result.getLength(),
                     result.getNumGaps(),
                     result.getRowFrom(), result.getRowTo(),
                     result.getColFrom(), result.getColTo()))
        
        if result.getLength() == 0:
            query_profile.useSegment()
            sbjct_profile.useSegment()
            return False, result, ("na",)
        
        elif result.getScore() < self.mMinAlignmentScore:
            query_profile.useSegment()
            sbjct_profile.useSegment()
            return False, result, ("na",)

        elif result.getScore() > self.mSafetyThreshold * self.mMinAlignmentScore:
            query_profile.useSegment()
            sbjct_profile.useSegment()
            return True,result, ("na",)
        
        z_params = alignlib.makeNormalDistributionParameters()
        alignlib.calculateZScoreParameters( z_params,
                                            query_profile,
                                            sbjct_profile,
                                            self.mAlignator,
                                            self.mNumIterationsZScore)
        
        mean   = z_params.getMean()
        stddev = z_params.getStandardDeviation()
        if stddev == 0: stddev = 1
        
        zscore = (result.getScore() - mean) / stddev
        
        self.debug( "--> mean=%f, stdev=%f, zscore=%f" % (mean, stddev, zscore) )
        
        query_profile.useSegment()
        sbjct_profile.useSegment()
        
        if zscore > self.mMinZScore:
            return True, result, ( "%5.2f" % zscore,)
        else:
            return False, result, ( "%5.2f" % zscore,)

