import sys, os, re, time, math, copy, random, glob  
import numpy
import scipy, scipy.stats, scipy.optimize
import matplotlib, pylab

from AddaModule import AddaModuleRecord
import AddaIO
import SegmentedFile
import multiprocessing

class Parameter:
    def __init__(self, value):
        self.value = value
    def set(self, value):
        self.value = value
    def __call__(self):
        return self.value

def fit(function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)
        
    if x is None: x = numpy.arange(y.shape[0])
    p = [param() for param in parameters]
    return scipy.optimize.leastsq(f, p)

class AddaFit( AddaModuleRecord ):
    """fit domains of a reference set to alignments and compute 
    parameters for ADDA's objective function.

    Briefly, each alignment between a pair of sequences is evaluated
    with respect of the domains in the sequences computing ``overhang``
    and ``transfer``.

    A domain might have ``overhang``, if it overlaps an alignment incompletely. 
    ``overhang`` is measured as the number of residues that are left uncovered. 

    ``transfer`` are the number of residues that the alignment links between any
    pair of domains of the same family in the two sequences. 

    input
       ``files:input_graph``: the pairwise alignment graph

       ``files:input_reference``: a reference domain definition

    output
       ``files:output_fit``: a config file with the estimated parameters

       ``files:output_fit_transfer``: a tab-separated histogram of transfer
          values.

       ``files:output_fit_overhang``: a tab-separated histogram of overhang
          values.

       ``files:output_fit_details``: details of the fitting procedure. This tab-separated
           table reports the transfer and overhang for combination of domains and alignment.

           class 
              domain family
           nid1
              sequence nid1
           dfrom1  
              domain start on nid1
           dto1  
              domain end on nid1
           afrom
              alignment start on nid1
           ato1 
              alignment end on nid1
           nid2
              sequence nid2
           dfrom
              domain start on nid2
           dto2
              domain end on nid2
           afrom2  
              alignment start on nid2
           ato2
              alignment end on nid2
           lali
              alignment length
           lx 
              length of domain on nid1
           ly 
              length of domain on nid2
           trans
              transfer value
           ptran
              percentage transfer (transfer/lali)
           atran
              average percentage transfer (transfer/ sqrt( lx * ly))
           score 
              alignment score (ln(evalue))
    """

    mName = "Fit"
    
    def __init__(self, *args, **kwargs ):

        AddaModuleRecord.__init__( self, *args, **kwargs )

        self.mFilenameFit = self.mConfig.get("files","output_fit", "adda.fit" )
        self.mFilenameOverhang = self.mConfig.get( "files", "output_fit_overhang" )
        self.mFilenameTransfer = self.mConfig.get( "files", "output_fit_transfer" )
        self.mFilenameDetails = self.mConfig.get( "files", "output_fit_details" )
        self.mMinTransfer = float(self.mConfig.get( "fit", "min_transfer" ))
        self.mMinOverhang = float(self.mConfig.get( "fit", "min_overhang" ))
        self.mFilenameNids = self.mConfig.get( "files", "output_nids", "adda.nids" )
        self.mFilenames = (self.mFilenameFit, self.mFilenameTransfer, self.mFilenameDetails, self.mFilenameOverhang )

        self.mOutfileDetails = None
    #--------------------------------------------------------------------------        
    def startUp( self ):

        if self.isComplete(): return

        self.mMapId2Nid = AddaIO.readMapId2Nid( open( self.mFilenameNids, "r") )

        self.mDomainBoundaries = {}

        rx_include = re.compile( self.mConfig.get( "fit", "family_include", "") )

        self.info( "reading domains from %s" % self.mConfig.get( "files", "input_reference") )
        infile = AddaIO.openStream( self.mConfig.get( "files", "input_reference") )

        ninput, nskipped_nid, nskipped_family, ndomains = 0, 0, 0, 0

        for line in infile:
            if line[0] == "#": continue
            ninput += 1
            token, start, end, family = line[:-1].split( "\t" )[:4]

            try:
                token = self.mMapId2Nid[token]
            except KeyError:
                nskipped_nid += 1
                continue

            if not rx_include.search( family): 
                nskipped_family += 1
                continue

            start, end = int(start), int(end)
            if token not in self.mDomainBoundaries:
                a = { family : [ (start, end) ] }
                self.mDomainBoundaries[token] = a
            else:
                a = self.mDomainBoundaries[token]
                if family not in a:
                    a[family] = [ (start, end) ]
                else:
                    a[family].append( (start,end) )
            ndomains += 1

        self.info( "read domain information: nsequences=%i, ndomains=%i, ninput=%i, nskipped_nid=%i, nskipped_family=%i" %\
                       (len(self.mDomainBoundaries), ndomains, ninput, nskipped_nid, nskipped_family))

        # result containers
        self.mTransferValues = []
        self.mOverhangValues = []

        # used to store data from an aborted run
        self.mValues = []
        self.mContinueAt = None

        # extract options
        self.mOutfileDetails  = self.openOutputStream( self.mFilenameDetails, register = True )

        if not self.mContinueAt:
            self.mOutfileDetails.write( """# FAMILY:          domain family
# NID1:         sequence nid1
# DFROM1:       domain start on nid1
# DTO1:         domain end on nid1
# AFROM1:       ali start on nid1
# ATO1:         ali end on nid1
# NID2:         sequence nid2
# DFROM2:       domain start on nid2
# DTO2:         domain end on nid2
# AFROM2:       ali start on nid2
# ATO2:         ali end on nid2
# LALI:         alignment length
# LX:           length of domain on nid1
# LY:           length of domain on nid2
# TRANS:        transfer value
# PTRAN:        percentage transfer (transfer/lali)
# ATRAN:        average percentage transfer (transfer/ sqrt( LX * LY))
# SCORE:        score of alignment
class\tnid1\tdfrom1\tdto1\tafrom1\tato1\tdnid2\tdfrom2\tdto2\tafrom2\tato2\tlali\tlx\tly\ttrans\tptran\tatran\tscore\n""")
            ## flushing is important with multiprocessing - why?
            ## if not flushed, the header and the EOF token appear twice.
            self.mOutfileDetails.flush()

    #--------------------------------------------------------------------------        
    def registerExistingOutput(self, filename):    

        if os.path.exists(filename):
            self.readPreviousData( filename )
            self.info("processing will continue after %s" % (str( self.mContinueAt ) ) )

    #--------------------------------------------------------------------------
    def processValues(self, values):
        """process data for a single query.
        
        The results are appended to self.mTransferValues and self.mOverhangValues
        """

        values.sort()

        def __iterate( values):
            t, nt, o, no = 0, 0, 0, 0
            last_f, last_q, last_s = None, None, None
            
            for f, q, s, transfer, overhang1, overhang2 in values:
                
                if f != last_f or \
                    q != last_q or \
                    s != last_s:
                    if nt > 0:
                        a = float(t)/nt
                    else:
                        a = None
                    
                    if no > 0:
                        b = float(o)/no
                    else:
                        b = None
                        
                    yield a,b 
            
                    t, nt, o, no = 0, 0, 0, 0
    
                    last_f, last_q, last_s = f, q, s
                    
                if transfer >= self.mMinTransfer:                                         
                    t += transfer
                    nt += 1
                    
                if overhang1 >= self.mMinOverhang:
                    o += overhang1
                    no += 1     
    
                if overhang2 >= self.mMinOverhang:
                    o += overhang2
                    no += 1     
                
            if nt > 0:
                a = float(t)/nt
            else:
                a = None
            
            if no > 0:
                b = float(o)/no
            else:
                b = None
                
            yield a,b 
                
            raise StopIteration

        for transfer, overhang in __iterate( values):                   

            if transfer:
                self.mTransferValues.append( transfer )
            if overhang:
                self.mOverhangValues.append( overhang )

    #--------------------------------------------------------------------------    
    def readPreviousData(self, filename = None):
        """process existing output in filename to guess correct point to continue computation."""
        
        if filename == None: filename = self.mFilenameDetails

        self.info( "reading previous data from %s" % filename )
        
        infile = open( filename, "r" )

        self.mTransferValues = []
        self.mOverhangValues = []
        
        def iterate_per_query(infile):
            
            last_query = None
            
            for line in infile:
                if line.startswith("#"): continue
                if line.startswith("class"): continue

                try: 
                    (family,
                     query_token, xfrom, xto, query_from, query_to,
                     sbjct_token, yfrom, yto, sbjct_from, sbjct_to,
                     lali, lx, ly, transfer, A, B, evalue) = line[:-1].split("\t")
                except ValueError:
                    self.warn( "parsing error in line %s\n" % line[:-1] )
                    continue
                
                if query_token != last_query:
                    if last_query: yield values
                    values = []
                    last_query = query_token
                    
                transfer, lx, ly = map( int, (transfer, lx, ly) )
                
                if transfer >= 0:
                    values.append( (family, query_token, sbjct_token, transfer, lx-transfer, ly-transfer) ) 
                    
            if last_query: 
                yield values
                self.mContinueAt = (query_token, sbjct_token)

            raise StopIteration
            
        for values in iterate_per_query(infile):
            self.processValues(values)
            
        self.info("read previous data from %s: transfer=%i, overhang=%i" % \
                      (filename, len(self.mTransferValues), len(self.mOverhangValues) ))
            
        infile.close()

    #--------------------------------------------------------------------------                    
    def applyMethod(self, neighbours ):
        """estimate ADDA penalties.

        This method calculates the distribution of::

           lali / sqrt( d1 * d2 ) 

        The computation only concerns domains of the same class.
            
        For each class:
        get all nids that have that class
        for each pair of nids, check if there is a link

        Repeats cause some counts to be inflated (most pronounced with the immunoglobulins)

        For example:
          * nid1: 3 domains
          * nid2: 2 domains
        
        Alignments: depending on domain arrangement 1 to 3 * 2, or:
          * nid1: 2 domains
          * nid2: 1 domain
          * alignments: 1 or 2

        If you want to eliminate repeats: which one?
            
        This method normalizes per family and per sequence pair.
        """

        values = []

        for n in neighbours.mMatches:
            if n.mQueryToken not in self.mDomainBoundaries or \
                n.mSbjctToken not in self.mDomainBoundaries: continue
 
            if self.mContinueAt:
                if (n.mQueryToken,n.mSbjctToken) == self.mContinueAt:
                    self.info("continuing processing at pair %s" % str(self.mContinueAt ) )
                    self.mContinueAt = None
                continue

            qdomains = self.mDomainBoundaries[n.mQueryToken]
            sdomains = self.mDomainBoundaries[n.mSbjctToken]
            
            for family in set(qdomains.keys()).intersection( set(sdomains.keys())):
                xdomains = qdomains[family]
                ydomains = sdomains[family]
                
                total_transfer = 0
                ntransfer = 0
                total_overhang = 0
                noverhang = 0
                
                for xfrom, xto in xdomains:
                    ovlx = min(xto,n.mQueryTo) - max(xfrom,n.mQueryFrom)
                    if ovlx < 0: continue                            
                    lx = xto - xfrom
                    for yfrom, yto in ydomains:
                        ovly = min(yto,n.mQueryTo) - max(yfrom,n.mQueryFrom)
                        if ovly < 0: continue

                        lali = min(n.mSbjctTo - n.mSbjctFrom, n.mQueryTo - n.mQueryFrom)
                        ly = yto - yfrom

                        zfrom = max(xfrom - n.mQueryFrom + n.mSbjctFrom, n.mSbjctFrom)
                        zto   = min(xto   - n.mQueryFrom + n.mSbjctFrom, n.mSbjctTo)                                
                        transfer = min(zto, yto) - max(zfrom, yfrom)

                        A = float(transfer) / float( lali )
                        B = float(transfer) / math.sqrt( float(lx * ly))
                        
                        self.mOutfileDetails.write( "\t".join( \
                                map(str, (family,
                                          n.mQueryToken, xfrom, xto, n.mQueryFrom, n.mQueryTo,
                                          n.mSbjctToken, yfrom, yto, n.mSbjctFrom, n.mSbjctTo,
                                          lali, lx, ly, transfer, A, B, n.mEvalue) )) + "\n" )
                        self.mOutfileDetails.flush()

                        if transfer >= 0:
                            values.append( (family, n.mQueryToken, n.mSbjctToken, transfer, lx-transfer, ly-transfer) ) 
                                         
        values.sort()
        self.processValues( values )
                                        
    #--------------------------------------------------------------------------
    def writeHistogram(self, outfile, intervals, frequencies ):
        
        for bin, value in zip(intervals, frequencies ):
            outfile.write( "%i\t%f\n" % (bin, value) )

    #--------------------------------------------------------------------------
    def cumulateHistogram( self, values ):
        t = 0
        v = []
        for x in values:
            t += x
            v.append( t )
        return v

    #--------------------------------------------------------------------------    
    def getCumulativeHistogram(self, values, reverse = False):

        x = numpy.arange( int(min(values)), int(max(values)) )
        h = scipy.stats.histogram2( values, x)
        if reverse:
            c = numpy.add.accumulate( numpy.array( h[::-1], numpy.float) )
        else: 
            c = numpy.add.accumulate( numpy.array( h, numpy.float) )
        total = max(c)
        y = c / total
        if reverse: y = y[::-1].copy()
        return x,y
    
    #--------------------------------------------------------------------------
    def plotHistogram(self, bins, vals,
                      title = None, 
                      filename = None,
                      f = None):
            
        # do not plot if called in subprocess. The first time this
        # function is called in a subprocess, it is fine, but called
        # by another, the error
        #
        # adda.py: Fatal IO error 0 (Success) on X server :0.0.
        #
        # appears.
        #
        # The test is not pretty and maintainable, but I could not
        # find how best to test if within MainProcess or not.
        if not re.search( "MainProcess", str(multiprocessing.current_process())):
            return

        pylab.plot( bins, vals )
        if f: 
            xstart, xend = pylab.gca().get_xlim()
            increment = (xend - xstart) / 100.0
            xvals = numpy.arange( xstart, xend, increment )
            yvals = f( xvals )
            pylab.plot( xvals, yvals )

        if title: pylab.title( title )
        pylab.xlabel( "residue" )
        pylab.ylabel( "relative cumulative frequency" )

        if filename:
            pylab.savefig( os.path.expanduser(filename) )
        else:
            pylab.show()
            
        pylab.clf()

    #--------------------------------------------------------------------------    
    def finish(self):

        self.mOutfile = self.openOutputStream( self.mFilenameFit, register = False )
        self.mOutfileTransfer = self.openOutputStream( self.mFilenameTransfer, register = False )
        self.mOutfileOverhang = self.openOutputStream( self.mFilenameOverhang, register = False )        
        
        self.info( "number of values: transfer=%i, overhang=%i" % (len(self.mTransferValues),
                                                                   len(self.mOverhangValues)) )
        
        x,y = self.getCumulativeHistogram(self.mTransferValues, reverse = True )
        self.writeHistogram( self.mOutfileTransfer, x, y )

        def f(x):
            return A() + B() * numpy.exp ( - numpy.exp( -(x - C()) / K() ) - (x - C()) / K() + 1 )

        A = Parameter(0.0)
        B = Parameter(1.0)
        C = Parameter(76.0)
        K = Parameter(8.6)

        result = fit(f,[A,B,C,K], y=y, x=x)

        self.plotHistogram(x, y, 
                           f = f,
                           filename = self.mFilenameTransfer + ".png",
                           title = "transfer" )


        # fit an exponential function to overhangs
        x,y = self.getCumulativeHistogram(self.mOverhangValues, reverse = True)
        self.writeHistogram( self.mOutfileOverhang, x, y )
        
        def f(x): return F() * numpy.exp ( -(x)* E() )
        
        E = Parameter(0.05)
        F = Parameter(1.0)
          
        result = fit(f,[E,F], y=y, x=x)
        
        self.plotHistogram(x, y, 
                           f = f,
                           filename = self.mFilenameOverhang + ".png",
                           title = "overhang" )
        
        self.mOutfile.write( "[optimise]\n" )
        self.mOutfile.write( "sigmoid_min=%f\n" % A() )
        self.mOutfile.write( "sigmoid_max=%f\n" % B() )
        self.mOutfile.write( "sigmoid_k=%f\n" % K() )
        self.mOutfile.write( "sigmoid_c=%f\n" % C() )
        self.mOutfile.write( "exponential_E=%f\n" % E() )                                
        self.mOutfile.write( "exponential_F=%f\n" % F() )

        self.mOutfile.close()
        self.mOutfileTransfer.close()
        self.mOutfileOverhang.close()

        ## close here, so that all is flushed before merge is called
        if self.mOutfileDetails: self.mOutfileDetails.close()

        AddaModuleRecord.finish( self )
        
    #--------------------------------------------------------------------------
    def merge(self):
        """merge runs from parallel computations.
        """

        # remove unwanted results
        for x in (self.mFilenameTransfer, self.mFilenameOverhang, self.mFilenameFit):
            for fn in glob.glob( "%s.0*" % x ):
                os.remove(fn)

        # merge the details file and compute stats
        if SegmentedFile.merge( self.mFilenameDetails ):
            self.readPreviousData( self.mFilenameDetails )
            self.finish()
            
