import sys, os, re, time, math, copy, glob, optparse, logging
import Numeric

from AddaModule import AddaModule
import CorrespondenceAnalysis
import MatlabTools
import SegmentedFile
import Experiment as E

class AddaSegment( AddaModule ):

    mName = "Segment"

    def __init__(self, *args, **kwargs ):

        AddaModule.__init__( self, *args, **kwargs )

        self.mFilenameSegments = self.mConfig.get("files","output_segments", "adda.segments" ) 
        self.mFilenames = (self.mFilenameSegments, )

    #--------------------------------------------------------------------------
    def startUp( self ):

        self.mHeaders = ("nid","node","parent","level","start","end") 
        if not self.isComplete():
            self.mOutfile = self.openOutputStream( self.mFilenameSegments )
            self.mMinDomainSize = int(self.mConfig.get('adda','min_domain_size'), 30)
            self.mResolution = float( self.mConfig.get('segments','resolution', 10.0) )
            
            if self.mContinueAt == None:
                self.mOutfile.write( "\t".join( self.mHeaders) + "\n" )
                self.mOutfile.flush()

    #--------------------------------------------------------------------------
    def finish(self):
        self.mOutfile.close()        
        AddaModule.finish( self )

    #--------------------------------------------------------------------------
    def validate(self):
        
        infile = SegmentedFile.fileopen( self.mFilenameSegments )

        last_nid = None
        found = set()
        nfound, nunknown, nduplicate = 0, 0, 0
        for line in infile:
            ninput += 1
            nid = line[:line.index("\t")]
            if nid != last_nid:
                if nid in found:
                    nduplicates += 1
                    self.warn("duplicate nid: %i in file %s" % (nid, filename))
                if nid not in tokens:
                    nunknown += 1
                    self.warn("unknown nid: %i in file %s" % (nid, filename))
                found.add(nid)
                nfound += 1
                last_nid = nid
            noutput += 1

        missing = set(self.mFasta.getTokens()).difference( found ) 
        if len(missing) > 0:
            self.warn( "the following nids were missing: %s" % str(missing) )

        self.info( "merging: ninput=%i, noutput=%i, nfound=%i, nmissing=%i, nduplicate=%i, nunknown=%i" %\
                       (ninput, noutput, nfound, len(missing), nduplicate, nunknown ) )
        
        return len(missing) == 0 and nduplicate == 0 and nunknown == 0
        
    #--------------------------------------------------------------------------
    def readPreviousData(self, filename ):
        """process existing output in filename to guess correct point to continue computation."""
        
        self.info( "reading previous data from %s" % filename )
        
        infile = open( filename, "r" )
        
        def iterate_per_query(infile):
            
            last_nid = None
            
            for line in infile:

                try: 
                    (nid,node,parent,level,xfrom,xto) = line[:-1].split("\t")
                except ValueError:
                    self.warn( "parsing error in line %s\n" % line[:-1] )
                    continue
                
                if nid != last_nid:
                    if last_nid: yield last_nid, values
                    values = []
                    last_nid = nid
                    
                values.append( map( int, (nid,node,parent,level,xfrom,xto) ) )
                    
            if last_nid: yield last_nid, values
            
            raise StopIteration
            
        nnids = 0
        for nid, values in iterate_per_query(infile):
            nnids += 1
            self.mContinueAt = nid
        
        self.info("found %i nids" % (nnids,)  )
            
        infile.close()

    #--------------------------------------------------------------------------
    def applyMethod(self, neighbours ):
        """split a sequence into putative domains."""
        
        nid = neighbours.mQueryToken

        if self.mContinueAt:
            if nid == self.mContinueAt:
                self.info("continuing processing at %s" % str(self.mContinueAt ) )
                self.mContinueAt = None
            return

        length = self.mFasta.getLength( nid )
        
        if length > self.mConfig.get("segments", "max_sequence_length"):
            self.warn( "skipped: nid=%s, length=%i -> too long" % (nid, length) )
            return False
        
        self.debug( "starting nid=%s, length=%i" % (nid, length) )

        tree = self.getTree( nid, length, neighbours ) 

        if self.mConfig.get("segments", "covering_trees"):
            tree = self.convertTreeToCoveringTree( tree )
        
        max_depth = 0
        if tree:
            for node in range(len(tree)):
                (level, parent, left_child, right_child, ranges) = tree[node]
                max_depth = max(max_depth, level)
                for xfrom, xto in ranges:
                    self.mOutfile.write("\t".join( map(str,(nid,node,parent,level,xfrom,xto)))+ "\n")
                self.mOutfile.flush()

            self.debug( "finished: nid=%s, length=%i, size=%i, depth=%i" % (nid, length,
                                                                          len(tree), max_depth ) )
            return True
        else:
            self.warn( "failed: nid=%s, length=%i, time=%i" % (job_id, nid, length,
                                                              len(tree), max_depth ) )
                                                                                        
            return False
        
    #--------------------------------------------------------------------------
    def normalizeMatrix( self, old_matrix, sums ):
        """normalize matrix. Each element will by x2 / row / col.
        
        The matrix is converted to floats to prevent integer
        overflows and returned as percent.
        """
        matrix = old_matrix.astype(Numeric.Float)
        Numeric.multiply( matrix, matrix, matrix )
        Numeric.divide( matrix, sums, matrix )        
        matrix = Numeric.transpose( matrix )
        Numeric.divide( matrix, sums, matrix )
        matrix *= 100
        return matrix.astype( Numeric.Int )

    #--------------------------------------------------------------------------
    def addLocalBiasToMatrix( self, matrix ):
        """adds a local bias to the correlation matrix.
        Adds diagonal elements.
        """
    
        r,c=matrix.shape
        
        bias = float( self.mConfig.get('segments','matrix_bias_strength'))
        width = int(self.mConfig.get('segments','matrix_bias_width'))
        
        for x in range(0,r):
            for y in range(max(x-width,0),
                           min(x+width,c)):
                matrix[x,y] = max(matrix[x,y], bias )

        
    #--------------------------------------------------------------------------
    def addChildren( self, tree ):
        """Add children to a tree.
        """
    
        children = []
        for x in range(0,len(tree)):
            children.append([])
    
        for  node, parent, level, ranges in tree:
            if node != 0:
                children[parent].append(node)
    
        new_tree = []
        for x in range(0,len(tree)):
            new_tree.append( [list(tree[x]), list(children[x])] )
    
        return new_tree
    
    #--------------------------------------------------------------------------
    def convertRanges2ExpandedRanges( self, ranges, max_length ):
        """remove the resolution factor from the ranges.
        """
        
        new_ranges = []
                
        for r in ranges:
            new_ranges.append( ( max(int(r[0] * self.mResolution), 0),
                                 min(int(r[1] * self.mResolution), max_length)) )
        return new_ranges
    #--------------------------------------------------------------------------
    def convertResidues2Ranges( self, residues, max_length ):
    
        ## skip over gaps
        residues = list(residues)
        residues.sort()
    
        last_r = residues[0]
        first_from = last_r 
        ranges = []
        for r in residues[1:]:
            if r - last_r > 1:
                ranges.append( (first_from, last_r) )
                first_from = r
            last_r = r 
            
        ranges.append( (first_from, last_r) )        
    
        return self.convertRanges2ExpandedRanges(ranges, max_length)
    
    #--------------------------------------------------------------------------
    def convertTreeToList( self, result, parent, tree, max_length ):
        """convert a tree to a list.
        """
        
        for t in tree:
    
            (level, residues, children) = t
    
            new_node = len(result)
    
            ranges = self.convertResidues2Ranges( residues, max_length )
    
            if ranges:
                result.append( (new_node, parent, level, ranges) )
                self.convertTreeToList( result, new_node, children, max_length)
    
    #--------------------------------------------------------------------------
    def resolveRanges( self, left_ranges, right_ranges):
        """removes small ranges by switching them between the two lists.
        """
        new_left_ranges  = []
        new_right_ranges = []
    
        ranges = map( lambda x: (x[0], x[1], 0), left_ranges)
        ranges += map( lambda x: (x[0], x[1], 1), right_ranges)
    
        ranges.sort()
        
        last_left, last_right, last_is_right = ranges[0]
        for this_left, this_right, this_is_right in ranges[1:]:
    
            ## if segment is the same type, just combine
            if (last_is_right and this_is_right) or (not last_is_right and not this_is_right):
                last_right = this_right
                continue
                
            ## write if not consecutive and there is a small gap
            if this_left - last_right > self.mConfig['segments']['min_segment_length']:
                if last_is_right:
                    new_right_ranges.append((last_left, last_right))
                else:
                    new_left_ranges.append((last_left, last_right))
    
                last_left, last_right, last_is_right = this_left, this_right, this_is_right
                continue
    
            ## if current segment is too small: add to current type
            if (this_right - this_left) < self.mConfig['segments']['min_segment_length']:
                last_right = this_right
                continue
    
            ## if previous segment is too small to be output: add to next type
            if (last_right - last_left) < self.mConfig['segments']['min_segment_length']:
                last_right = this_right
                last_is_right = this_is_right
                continue
    
            ## otherwise: output
            if last_is_right:
                new_right_ranges.append((last_left, last_right))
            else:
                new_left_ranges.append((last_left, last_right))
                
            last_left, last_right, last_is_right = this_left, this_right, this_is_right            
                
        if last_is_right:
            new_right_ranges.append((last_left, last_right))
        else:
            new_left_ranges.append((last_left, last_right))
    
        self.debug( "ranges=%s" % str(ranges), 4 )
        self.debug( "new_left_ranges=%s" % str(new_left_ranges), 4)
        self.debug( "new_right_ranges=%s" % str(new_right_ranges), 4 )
            
        return new_left_ranges, new_right_ranges
            
    #--------------------------------------------------------------------------
    def getCoveringRanges( self, left_ranges, right_ranges, parent_ranges ):
        """make sure, that left_ranges and right_ranges both cover the
        range given by parent_ranges completely.
    
        Removes small fragments as well.
        """
    
        child_ranges = map( lambda x: (x[0], x[1], 0), left_ranges)
        child_ranges += map( lambda x: (x[0], x[1], 1), right_ranges)
    
        child_ranges.sort()
        parent_ranges.sort()
    
        new_left_ranges = []
        new_right_ranges = []
    
        parent_index = 0
        last_to = 0
    
        parent_left, parent_right = parent_ranges[parent_index]

        self.debug( "child_ranges=%s" % str(child_ranges) )
        self.debug( "parent_ranges=%s" % str(parent_ranges))
    
        last_left, last_right, last_is_right = child_ranges[0]
        
        for this_left, this_right, this_is_right in child_ranges[1:]:
    
            ## look at previous segment last_left to last_right:
            ## find matching parent_index:
            old_parent_index = parent_index
            while (min(parent_right, last_right) - max(parent_left, last_left)) < 0:
                parent_index += 1
                if parent_index == len(parent_ranges): break
                parent_left, parent_right = parent_ranges[parent_index]
    
            ## skip fragments that do not overlap
            if parent_index == len(parent_ranges):
                parent_index = old_parent_index
                last_left, last_right, last_is_right = this_left, this_right, this_is_right
                continue
                
            ## firstly: make segment covering
            new_left  = min(parent_left, last_left)
            new_right = min(max(parent_right, last_right), this_left)
    
            if last_is_right:
                new_right_ranges.append((new_left, new_right))
            else:
                new_left_ranges.append((new_left, new_right))
    
            ## reduce parent on left side
            parent_left=max(new_right, parent_left)
    
            last_left, last_right, last_is_right = this_left, this_right, this_is_right
    
        ## process last segment
        while (min(parent_right, last_right) - max(parent_left, last_left)) < 0:
            parent_index += 1
            if parent_index >= len(parent_ranges): break        
            parent_left, parent_right = parent_ranges[parent_index]
    
        new_left = min(parent_left, last_left)
        new_right = max(parent_right, last_right)
            
        if last_is_right:
            new_right_ranges.append((new_left, new_right))
        else:
            new_left_ranges.append((new_left, new_right))
    
        self.debug( "old left ranges=%s" % str(left_ranges))
        self.debug( "new left ranges=%s" % str(new_left_ranges))
        self.debug( "old right ranges=%s" % str(right_ranges))
        self.debug( "new right ranges=%s" % str(new_right_ranges))
        
        return new_left_ranges, new_right_ranges
    
    #----------------------------------------------------------------
    def removeSmallRanges( self, ranges, min_segment_length, max_separation = 1):
        """resolve ranges.
    
        ranges are defined as tuples: (from, to, type)
    
        Small ranges are deleted and added to neighbouring
        domains.
        """
        
        ranges.sort()
    
        new_ranges = []
    
        last_left, last_right, last_type = ranges[0]
        
        for this_left, this_right, this_type in ranges[1:]:
    
            # print this_left, this_right, this_type,
            ## write if not consecutive and there is a small gap, but only
            ## if segment is long enough, otherwise: discard
            if this_left - last_right > max_separation:
                if (last_right - last_left) >= min_segment_length:
                    new_ranges.append((last_left, last_right, last_type))
                last_left, last_right, last_type = this_left, this_right, this_type
                # print "a"
                continue
    
            
            ## if segment is the same type as last type, just combine
            if last_type == this_type:
                last_right = this_right
                # print "b"
                continue
                
            ## if current segment is too small: add to last type
            if (this_right - this_left) < min_segment_length:
                last_right = this_right
                # print "c"
                continue
    
            ## if previous segment is too small to be output: add to current type
            if (last_right - last_left) < min_segment_length:
                last_right = this_right
                last_type = this_type
                # print "d"
                continue
    
            ## otherwise: output
            new_ranges.append((last_left, last_right, last_type))
                
            last_left, last_right, last_type = this_left, this_right, this_type
            
        new_ranges.append((last_left, last_right, last_type))
    
        self.debug( "# old ranges=%s" % ranges)
        self.debug( "# new ranges=%s" % new_ranges)
            
        return new_ranges
    
    ##--------------------------------------------------------
    def setChildren( self, tree ):
        """sets children correctly."""

        ## set all children to zero
        for node in range(len(tree)):
        
            level, parent, left_child, right_child, ranges = tree[node]
            
            tree[node][2] = 0
            tree[node][3] = 0
        
            if tree[parent][2]:
                tree[parent][3] = node
            else:
                tree[parent][2] = node

    ##--------------------------------------------------------
    def collapseTree(self, tree):
        """remove all empty levels in the tree and renumber nodes
        so that they are continuous.
        """
    
        self.setChildren(tree)
        map_old2new = {}
    
        index = 1
        new_tree = []
    
        ## write root
        new_tree.append( [0, 0, 0, 0, tree[0][4]] )
        map_old2new[0] = 0
    
        ## PrintTree(tree)
    
        for old_node in range(1, len(tree)):
            level, parent, left_child, right_child, ranges = tree[old_node]
    
            ## if only a single child of a parent, skip this node
            ## if ranges are empty: skip this node
            if tree[parent][2] == 0 or tree[parent][3] == 0 or len(ranges) == 0:
                map_old2new[old_node] = map_old2new[parent]
                continue
    
            map_old2new[old_node] = index
    
            ## add to new tree
            new_tree.append( [new_tree[map_old2new[parent]][0] + 1, map_old2new[parent], 0, 0, ranges] )
            
            index += 1
            
            ## PrintTree( new_tree )
    
        ## PrintTree( new_tree )
        ## print "#########"
        
        return new_tree

    #--------------------------------------------------------------------------    
    def printTree(self, tree):
        """output a tree."""

        for node in range(0, len(tree)):
            print "%i\t" % node + "\t".join( map(str, tree[node]))

    #--------------------------------------------------------------------------    
    def convertTreeToCoveringTree( self, tree ):
        """make a covering tree out of a splitting tree.
        Shortening of domains is not allowed.
        """
    
        self.debug( "making covering trees" )
    
        ntree = self.addChildren( tree )
    
        #######
        # descend tree and add new domains
        # if domain has only a single child: delete the child and
        # rewire
        for t in ntree:
            info, children = t
    
            if info:
                node, parent, level, ranges = info
    
            if len(children) == 1:
                ntree[children[0]][0] = None
                ntree[node][1] = ntree[children[0]][1]
                
        #######
        # build new tree with new node identifiers
        current_node = 0
        covering_tree = []
    
        levels = map( lambda x: [], [0] * len(tree))
        
        for t in ntree:
            info, children = t
    
            if not info: continue
            node, parent, level, ranges = info
            
            if len(children) == 2:
    
                # add new node to tree, rename parent in children and
                # set borders
                leftchild = children[0]
                rightchild = children[1]                
    
                # change left child
                lnode, lparent, llevel, lranges = ntree[leftchild][0]
                rnode, rparent, rlevel, rranges = ntree[rightchild][0]            
    
                if ranges:
                    lranges, rranges = self.getCoveringRanges( lranges, rranges, ranges )
                else:
                    continue
    
                # change left child
                ntree[leftchild][0]= (None, current_node, level + 1, lranges) 
    
                # change right child            
                # cnode, cparent, clevel, cranges = ntree[rightchild][0]
                ntree[rightchild][0]= (None, current_node, level + 1, rranges )
    
            covering_tree.append( [level, parent, 0, 0, ranges] )
            levels[level].append( current_node )
                
            current_node += 1
    
        max_range = covering_tree[0][4][0][1]
    
        self.debug( "resulting tree" )
        if E.getLogLevel() >= 2:
            self.printTree( covering_tree )
        
        ###################################
        ## remove small fragments
        ## has to be done per level in order to be consistent
        ## done here and not during matrix decomposition, so that
        ## matrix needs not to be permuted more than once.
        for l in range(0, len(levels)):
            if len(levels[l]) == 0: break
            # collect all domains per level in a list of the form
            # (from, to, node)
            ranges = []
            for node in levels[l]:
                ranges += map(lambda x: (x[0], x[1], node), covering_tree[node][4])
                covering_tree[node][4] = []
                
            # and remove small fragments
            new_ranges = self.removeSmallRanges( ranges, 
                                                 int(self.mConfig.get('segments','min_segment_size')), 
                                                 int(self.mConfig.get('segments','min_segment_size')))
    
            # and put back into tree if there is more than one range
            for (xfrom, xto, node) in new_ranges:
                covering_tree[node][4].append( (xfrom, xto) )
    
        ###################################
        ## delete nodes with empty ranges or only a single child.
        ## renumber nodes so that there are no gaps
    
        if E.getLogLevel() >= 2:
            self.printTree( covering_tree )
        
        return self.collapseTree( covering_tree )
    
    #-----------------------------------------------------------------------------------------------------
    def buildMatrix( self,
                     query_nid, lsequence, neighbours ):
        
        """build matrix based on BLAST alignments to query_nid.

        matrix of size N*M
        N: number of neighbours
        M: length of query (scaled with resolution)

        alignments are truncated.

        the query is included in the matrix.
    
        if combine_repeats is set, multiple alignments between the query and a sbjct will
        be entered into the same row.

        if residue_level is set, entries are added on the residue level. The resolution parameter
        is ignored.
        """

        combine_repeats = self.mConfig.get( "segments", "combine_repeats")
                    
        query_length = int( math.floor( float(lsequence )) / float(self.mResolution))

        nindex = {}
    
        nneighbours = 0
        if combine_repeats:
            for neighbour in neighbours.mMatches:
                if not nindex.has_key(neighbour.mSbjctToken):
                    nindex[neighbour.mSbjctToken] = nneighbours
                    nneighbours += 1
        else:
            nneighbours = len(neighbours)

        # build matrix and add query sequence
        nneighbours += 1
        matrix = Numeric.zeros( (nneighbours, query_length), Numeric.Int)    
        matrix[0, 0:query_length] = 1
        row = 1
        
        for n in neighbours.mMatches:

            if combine_repeats:
                use_row = nindex[n.mSbjctToken]
            else:
                use_row = row
                row += 1
    
            yfrom = int(math.floor(n.mQueryFrom/self.mResolution))
            yto   = int(math.floor(n.mQueryTo/self.mResolution)) 
            matrix[use_row, yfrom:yto] = 1
                
        return matrix
    
    
    #--------------------------------------------------------------------------
    def getTree( self, nid, lsequence, neighbours):
        
        ## retrieve blast matrix
        ## make sure, add_self is 1, so that there are no empty columns
        ## todo
        blast_matrix = self.buildMatrix( nid, lsequence, neighbours )
        
        self.debug( "blast matrix for %s: %s" % (nid, str(blast_matrix.shape)))

        if E.getLogLevel() >= 3:
            MatlabTools.WriteMatrix(blast_matrix, outfile=open("blast_%s.matrix" % nid, "w"))
    
        nneighbours, lmatrix = blast_matrix.shape

        self.debug( "rows in blast matrix for %s: %s" % (nid, str(nneighbours))) 
    
        if nneighbours < int(self.mConfig.get("segments", "min_neighbours")):
            return []
        
        ## calculate dot product of the matrix
        dot_matrix = Numeric.matrixmultiply( Numeric.transpose( blast_matrix ), blast_matrix, )
    
        ## perform some matrix magic
        if int(self.mConfig.get('segments','multiply')) > 0:
            for x in range(0, int(self.mConfig.get('segments','multiply'))):
                dot_matrix = Numeric.matrixmultiply( dot_matrix, dot_matrix)

                if E.getLogLevel() >= 3:
                    MatlabTools.WriteMatrix(dot_matrix, outfile=open("correlation_%s_%i.matrix" % (nid, x), "w"))
    
        if self.mConfig.get('segments','matrix_add_local_bias'):
            self.addLocalBiasToMatrix( dot_matrix )
    
        if self.mConfig.get('segments','normalize'):
            dot_matrix = self.normalizeMatrix( dot_matrix, Numeric.sum( blast_matrix ))
            
        if E.getLogLevel() >= 3:
            self.debug( "correlation matrix for %s: %s" % (nid, str(dot_matrix.shape)))            
            MatlabTools.WriteMatrix(dot_matrix, outfile=open("correlation_%s.matrix" % nid, "w"))
            
        ## rearrange matrix if necessary
        map_row_new2old = range(0, lmatrix)        
        if self.mConfig.get('segments','permute') == "True":
            row_indices, col_indices =  CorrespondenceAnalysis.GetIndices( dot_matrix )
            
            if not row_indices:
                print "# error for %i: correspondence analysis did not converge" % nid
            else:
                map_row_new2old = Numeric.argsort(row_indices)
            
            dot_matrix = CorrespondenceAnalysis.GetPermutatedMatrix( dot_matrix, map_row_new2old, map_row_new2old)
                
            if E.getLogLevel() >= 3:
                self.debug( "permuted matrix for %s: %s" % (nid, str(dot_matrix.shape)))            
                MatlabTools.WriteMatrix(dot_matrix, outfile=open("permuted_%s.matrix" % nid, "w"))
        
        full_range = dot_matrix.shape
    
    
        xtree = self.splitMatrix( nid,
                                  dot_matrix,
                                  (0, lmatrix),
                                  0,
                                  int(self.mMinDomainSize / self.mResolution),
                                  int(self.mConfig.get('segments','min_distance_border')),
                                  map_row_new2old)
    
        self.debug( "xtree=%s" % xtree )
    
        tree = [(0, 0, 0, [(0, lsequence)])]
        self.convertTreeToList( tree, 0, xtree, lsequence )
    
        self.debug( "tree=%s" % str(tree) )
    
        return tree


    #-------------------------------------------------------------------------
    def splitMatrix( self,
                     nid,
                     matrix,
                     intervall,
                     level,
                     min_length = 30,
                     min_distance_border = 0,
                     map_row_new2old = None ):
        """
        1. calculate objective function for matrix in intervall.
    
    
        \           <- xfrom
         \
          \
        c1 \
        ----\       <- x
           | \
        cc |c2\
                    <- xto (one past last)
    
        l1 = x - xfrom
        l2 = xto - x
    
        chi-squared: 
    
        I[x] = (i11*i22-i21*i12)**2 * total / row&col-sums
    
    
        I[x] = mu[x]/F[x]
    
        """
        xfrom, xto = intervall
        l = xto - xfrom 
    
        if l < min_length: return []
    
        self.debug( "# splitting matrix on level %i in intervall %s %i" % (level, intervall, min_length))
    
        ## 1. build Interfaces
        I = Numeric.zeros( l, Numeric.Float )
    
        for x in range(xfrom+1, xto-1):
    
            i11 = float(Numeric.sum(Numeric.sum(matrix[xfrom:x,xfrom:x])))
            i22 = float(Numeric.sum(Numeric.sum(matrix[x:xto,x:xto])))
            i12 = float(Numeric.sum(Numeric.sum(matrix[xfrom:x,x:xto] )))
            i21 = i12
    
            row1 = i11 + i12
            row2 = i21 + i22
            col1 = i11 + i21
            col2 = i12 + i22
    
            l1 = x-xfrom
            l2 = xto - x
    
            a = i11 * i22 - i21 * i12
    
            n = row1 * row2 * col1 * col2
            if n > 0.0:
                I[l1] = a * a / n
            else:
                I[l1] = 0.0

            self.debug( "# %i\t%i\t%i\t%i\t%i\t%i\t%f" % (x, l1, l2, i11, i22, i12, I[l1]))

        ## 2. split at maximum
        if min_distance_border and l > 2 * min_distance_border:
            xmax = Numeric.argmax( I[min_distance_border:(-min_distance_border)] ) + min_distance_border
        else:
            xmax = Numeric.argmax( I )
    
        val  = I[xmax]
        pos = xmax + xfrom
    
        if xmax == 0: return []

        self.debug( "splitting at position %i with value %f %f" % (pos, val, xmax))
    
        result = []
    
        if xmax > min_length:
            result += [(level+1,
                        map_row_new2old[xfrom:pos],
                        self.splitMatrix( nid,
                                     matrix,
                                     (xfrom, pos),
                                     level+1,
                                     min_length,
                                     min_distance_border,
                                     map_row_new2old))]
    
    
        if (l - xmax) > min_length:
            result += [(level+1,
                        map_row_new2old[pos:xto],
                        self.splitMatrix( nid,
                                     matrix,
                                     (pos,xto),
                                     level+1,
                                     min_length,
                                     min_distance_border,
                                     map_row_new2old))]
            
        return result
    
