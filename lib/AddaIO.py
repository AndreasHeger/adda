import sys, os, re, time, types, gzip
import alignlib
from ConfigParser import ConfigParser as PyConfigParser

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

class NeighbourRecord:

    def __init__(self, line):
    
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
    
    def getAlignment(self):
        """parse alignment into a AlignmentVector object."""
        r = alignlib.makeAlignmentVector()
        f = alignlib.AlignmentFormatEmissions()
        f.mRowFrom, f.mRowTo, f.mRowAlignment = self.mQueryFrom, self.mQueryTo, self.mQueryAli
        f.mColFrom, f.mColTo, f.mColAlignment = self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli     
        f.copy( r )
        return r   
        
class NeighbourIterator:

    def _iterate( self, infile ):
    
        for line in infile:

            if line[0] == "#": continue
    
            if not line.strip(): continue
            
            yield NeighbourRecord( line )
        
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

    def _iterate( self, infile ):
    
        last_nid = None
    
        iterator = NeighbourIterator( infile )
        last_token = None
        
        while 1:

            r = iterator.next()
            if not r: break
            if self.mTokens and \
                (r.mQueryToken not in self.mTokens or \
                 r.mSbjctToken not in self.mTokens ):
                continue 
            if r.mQueryToken != last_token:
                if last_token:
                    yield NeighboursRecord( last_token, matches )
                matches = []
                last_token = r.mQueryToken
                
            matches.append( r )
            
        if last_token:
            yield NeighboursRecord( last_token, matches )
        raise StopIteration

    def __init__(self, f, tokens = None, *args, **kwargs):
        """
        f: the input file object.
        tokens: a collection of tokens to filter with.
        """
        
        self.mIterator = self._iterate(f)
        if tokens: self.mTokens = tokens

    def __iter__(self):
        return self

    def next(self):
        try:
            return self.mIterator.next()
        except StopIteration:
            return None

        
        
        