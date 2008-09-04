import sys, os, re, time, math, copy

from AddaModule import AddaModule

class AddaGraph( AddaModule ):
    """write a links table for adda processing."""
    
    mName = "Graph"
    
    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )

        self.mOutfile = self.openOutputStream( self.mConfig.get("files", "output_graph") )
                
        self.mMergeRepeats = self.mConfig.get( "graph", "merge_repeats") == "True"
        self.mMinDomainSize = int(self.mConfig.get('adda','min_domain_size'))
                
        self.mNJoined = 0
        self.mNLinksInput = 0
        self.mNLinksOutput = 0        
        
    def applyMethod(self, neighbours ):
        """output the graph.

        If mMergeRepeats is set, consecutive links are merged.
        Links are consecutive if they are adjacent both in the
        query and in the sbjct.
        
        This ensures that 1:many repeats are not merged, but will
        cover alignments split by transmembrane regions.
        
        """

        if len(neighbours.mMatches) == 0:
            return
        
        self.mNLinksInput += len(neighbours.mMatches)
        
        if self.mMergeRepeats:
            matches = []
            neighbours.mMatches.sort( lambda x,y: cmp( \
                    (x.mSbjctToken, x.mQueryFrom), (y.mSbjctToken, y.mQueryFrom )))
            last = neighbours.mMatches[0]
            
            for match in neighbours.mMatches[1:]:
                if match.mSbjctToken == last.mSbjctToken and \
                    0 < match.mQueryFrom - last.mQueryTo <= self.mMinDomainSize and \
                    0 < match.mSbjctFrom - last.mSbjctTo <= self.mMinDomainSize: 
                    self.mNJoined += 1
                    last.mEvalue = min( match.mEvalue, last.mEvalue )
                    self.debug( "joining: %s:%s-%s %s:%s-%s with %s:%s-%s %s:%s-%s" %\
                                (last.mQueryToken, last.mQueryFrom, last.mQueryTo,
                                 last.mSbjctToken, last.mSbjctFrom, last.mSbjctTo,
                                 match.mQueryToken, match.mQueryFrom, match.mQueryTo,
                                 match.mSbjctToken, match.mSbjctFrom, match.mSbjctTo))
                    
                else:   
                    matches.append(last)
                    last = match
                    continue
                    
                last.mQueryTo = max(last.mQueryTo, match.mQueryTo)
                last.mSbjctTo = max(last.mSbjctTo, match.mSbjctTo)                
                
            matches.append(last)
                 
        else:
            neighbours.mMatches.sort( lambda x,y: cmp( x.mSbjctToken, y.mSbjctToken ))
            matches = neighbours.mMatches()
            
        self.mNLinksOutput += len(matches)
        for m in matches:
            self.mOutfile.write( "%s\t%s\t%f\t%i\t%i\t%i\t%i\n" % \
                                 ( m.mQueryToken, 
                                   m.mSbjctToken, 
                                   m.mEvalue,
                                   m.mQueryFrom, m.mQueryTo,
                                   m.mSbjctFrom, m.mSbjctTo ) )    
            self.mOutfile.flush()
            
    def finish(self):
        
        self.mOutfile.close()
        
        self.info( "graph: %i links input, %i links output, %i links merged" %\
                   (self.mNLinksInput, self.mNLinksOutput, self.mNJoined ) )
        
        AddaModule.finish( self )
