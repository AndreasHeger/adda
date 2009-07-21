import sys, os, re, time, types, gzip
import alignlib
from ConfigParser import ConfigParser as PyConfigParser
import Experiment as E

class ConfigParser( PyConfigParser ):
    """config parser with defaults."""
    def __init__(self, *args, **kwargs):
        PyConfigParser.__init__(self, *args, **kwargs)
    
    def get( self, section, option, default = None):
        """Get an option value for the named section.
        
        If default is given and the option does not exist,
        default is returned. The default value determines
        the type.
        """
        if default != None:
            if not self.has_option( section, option ): 
                return default
                    
            t = type(default)
            if t is types.StringType:
                return PyConfigParser.get( self, section, option )
            elif t is types.IntType:
                return self.getint( section, option )
            elif t is types.FloatType:
                return self.getfloat( section, option )
            elif t is types.BooleanType:
                return self.getboolean( section, option )        
            raise TypeError, "unknown type %s" % t
        else:
            return PyConfigParser.get( self, section, option )

def openStream( filename ):
    """open an input stream.
    """
    
    if filename[-3:] == ".gz":
        return gzip.open( filename, "r" )
    else:
        return open( filename, "r" )

class NeighbourRecordPairsdb:
    """a pairwise alignment.

    The alignment is parsed from the input line.

    The input format is tab-separated columns:

    ``query_token`` the query
    ``sbjct_token`` the sbjct
    ``evalue`` : the E-Value
    ``query_from``: the first aligned residue in query
    ``query_to``: the last aligned residue + 1 in query
    ``query_ali``: the aligned query in compressed form
    ``sbjct_from``: the first aligned residue in sbjct
    ``sbjct_to``: the last aligned residue + 1 in sbjct
    ``sbjct_ali``: the aligned sbjct in compressed form

    Additional columns are ignored.
    """

    def __init__(self, line ): 
        (self.mQueryToken, self.mSbjctToken, self.mEvalue,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.mEvalue = float(self.mEvalue)

    def __str__( self ):

        return "\t".join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.mEvalue,
            self.mQueryFrom, self.mQueryTo, self.mQueryAli,
            self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli)))
    
    def getAlignment(self ):
        """parse alignment into a AlignmentVector object."""
        r = alignlib.makeAlignmentVector()
        f = alignlib.AlignmentFormatEmissions()
        f.mRowFrom, f.mRowTo, f.mRowAlignment = self.mQueryFrom, self.mQueryTo, self.mQueryAli
        f.mColFrom, f.mColTo, f.mColAlignment = self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli     
        f.copy( r )
        return r   

class NeighbourRecordPairsdbRealign:
    """a pairwise alignment.

    The input format is tab-separated columns:

    ``query_token`` the query
    ``sbjct_token`` the sbjct
    ``evalue`` : the E-Value
    ``query_from``: the first aligned residue in query
    ``query_to``: the last aligned residue + 1 in query
    ``query_ali``: the aligned query in compressed form
    ``sbjct_from``: the first aligned residue in sbjct
    ``sbjct_to``: the last aligned residue + 1 in sbjct
    ``sbjct_ali``: the aligned sbjct in compressed form

    Additional columns are ignored and the alignment itself
    are ignored forcing ADDA to realign.
    """

    def __init__(self, line ): 
        (self.mQueryToken, self.mSbjctToken, self.mEvalue,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.mEvalue = float(self.mEvalue)

    def __str__( self ):

        return "\t".join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.mEvalue,
            self.mQueryFrom, self.mQueryTo, self.mQueryAli,
            self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli)))
    
    def getAlignment(self ):
        """parse alignment into a AlignmentVector object."""
        return None

class NeighbourRecordSimap:
    """a pairwise alignment.

    The alignment is built on demand from the coordinates
    by re-alignment.

    The input format is tab-separated columns:

    ``query_token`` the query
    ``sbjct_token`` the sbjct
    ``evalue`` : the E-Value
    ``query_from``: the first aligned residue in query
    ``query_to``: the last aligned residue + 1 in query
    ``sbjct_from``: the first aligned residue in sbjct
    ``sbjct_to``: the last aligned residue + 1 in sbjct

    Additional columns are ignored.
    """

    def __init__(self, line ): 
        (self.mQueryToken, self.mSbjctToken, 
         self.mEvalue,
         self.mQueryFrom, self.mQueryTo,
         self.mSbjctFrom, self.mSbjctTo) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.mEvalue = float(self.mEvalue)
        self.mAlignment = None

    def __str__( self ):

        return "\t".join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.mEvalue,
            self.mQueryFrom, self.mQueryTo,
            self.mSbjctFrom, self.mSbjctTo )))
    
    def getAlignment(self, fasta ):
        """parse alignment into a AlignmentVector object."""
        return None

class NeighbourRecordPairsdbOld(NeighbourRecordPairsdb):
    """a pairwise alignment in old pairsdb format.

    The old pairsdb format used one-based coordinates.
    """
    def __init__(self, line ): 
        NeighbourRecordPairsdb.__init__( self, line )
        self.mQueryFrom -= 1
        self.mSbjctFrom -= 1

class NeighbourIterator:

    def _iterate( self, infile, record = NeighbourRecordPairsdb ):
    
        for line in infile:

            if line[0] == "#": continue
            if not line.strip(): continue
            yield record ( line )
        
        raise StopIteration

    def __init__(self, f, *args, **kwargs):
        self.mIterator = self._iterate(f)

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.mIterator.next()
        except StopIteration:
            return None

class NeighboursRecord:
    def __init__(self, token, matches):
        self.mQueryToken = token
        self.mMatches = matches

class NeighboursIterator:

    def __init__(self, f, map_id2nid = None, *args, **kwargs):
        """
        f: the input file object.
        tokens: a collection of tokens to filter with.
        """
        
        self.mIterator = self._iterate(f)
        self.mMapId2Nid = map_id2nid

    def _iterate( self, infile ):
    
        last_nid = None
    
        iterator = NeighbourIterator( infile, self.mRecord )
        last_token = None
        
        while 1:

            r = iterator.next()
            if not r: break
            if self.mMapId2Nid:
                if (r.mQueryToken not in self.mMapId2Nid or \
                        r.mSbjctToken not in self.mMapId2Nid ):
                    continue 
                r.mQueryToken = self.mMapId2Nid[r.mQueryToken]
                r.mSbjctToken = self.mMapId2Nid[r.mSbjctToken]
            if r.mQueryToken != last_token:
                if last_token:
                    yield NeighboursRecord( last_token, matches )
                matches = []
                last_token = r.mQueryToken
                
            matches.append( r )
            
        if last_token:
            yield NeighboursRecord( last_token, matches )
        raise StopIteration

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.mIterator.next()
        except StopIteration:
            return None

class NeighboursIteratorPairsdb( NeighboursIterator ):
    """iterate over Pairsdb formatted file."""
    mRecord = NeighbourRecordPairsdb

class NeighboursIteratorSimap(  NeighboursIterator ):
    """iterate over SIMAP formatted file."""
    mRecord = NeighbourRecordSimap

    def _iterate( self, infile ):
    
        last_nid = None
    
        iterator = NeighbourIterator( infile, self.mRecord )
        last_token = None
        
        alignator = alignlib.makeAlignatorDPFull( alignment.ALIGNMENT_LOCAL, 
                                                  -10, -2)

        q,s = None, None

        while 1:

            r = iterator.next()
            if not r: break
            if self.mMapId2Nid:
                if (r.mQueryToken not in self.mMapId2Nid or \
                        r.mSbjctToken not in self.mMapId2Nid ):
                    continue 
                r.mQueryToken = self.mMapId2Nid[r.mQueryToken]
                r.mSbjctToken = self.mMapId2Nid[r.mSbjctToken]

            if r.mQueryToken != last_token:
                if last_token:
                    yield NeighboursRecord( last_token, matches )
                matches = []
                last_token = r.mQueryToken
                q = alignlib.makeSequence( fasta.getSequence( r.mQueryToken ) )

            # do a re-alignment
            s = alignlib.makeSequence( fasta.getSequence( r.mSbjctToken ) )
            q.useSegment( r.mQueryFrom, r.mQueryTo )
            s.useSegment( r.mSbjctFrom, r.mSbjctTo )
            ali = alignlib.makeAlignmentVector()

            alignator.align( ali, q, s )
            r.mAlignment = ali
            matches.append( r )
            
        if last_token:
            yield NeighboursRecord( last_token, matches )
        raise StopIteration


def readMapId2Nid( infile ):
    """read map from adda.nids file."""
    
    m = {}
    for line in infile:
        if line.startswith("#"): continue
        data = line[:-1].split("\t")[:2]
        m[data[1]] = data[0]
    return m

def readMapNid2Domains( infile, map_id2nid, rx_include ):
    """read reference domain file.
    
    Only include families matching the regulare expression rx_include.
    """

    domain_boundaries = {}

    rx_include = re.compile( rx_include )

    ninput, nskipped_nid, nskipped_family, ndomains = 0, 0, 0, 0

    for line in infile:
        if line[0] == "#": continue
        ninput += 1
        token, start, end, family = line[:-1].split( "\t" )[:4]

        try:
            token = map_id2nid[token]
        except KeyError:
            nskipped_nid += 1
            continue

        if not rx_include.search( family): 
            nskipped_family += 1
            continue

        start, end = int(start), int(end)
        if token not in domain_boundaries:
            a = { family : [ (start, end) ] }
            domain_boundaries[token] = a
        else:
            a = domain_boundaries[token]
            if family not in a:
                a[family] = [ (start, end) ]
            else:
                a[family].append( (start,end) )
        ndomains += 1

    E.info( "read domain information: nsequences=%i, ndomains=%i, ninput=%i, nskipped_nid=%i, nskipped_family=%i" %\
                (len(domain_boundaries), ndomains, ninput, nskipped_nid, nskipped_family))

    return domain_boundaries
