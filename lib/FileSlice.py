import math, gzip, struct, itertools
import Experiment as E

def getZipSize(gzipfile):
    f = open(gzipfile, "rb")
    if f.read(2) != "\x1f\x8b":
        raise IOError("not a gzip file")
    f.seek(-4, 2)
    return struct.unpack("<i", f.read())[0]


class Iterator:

    def __init__(self, filename, nchunks, chunk, iterator, *args, **kwargs ):

        if filename.endswith(".gz"):
            self.mFileSize = getZipSize( filename )
            self.mInfile = gzip.open( filename, "r" )
        else:
            self.mInfile = open( filename, "r" )
            self.mInfile.seek(0, 2)
            self.mFileSize = self.mInfile.tell()

        chunk_size = int(math.ceil( float(self.mFileSize) / nchunks ) )
        start_pos = chunk_size * chunk
        self.mEndPos = start_pos + chunk_size

        self.mInfile.seek( start_pos )
        self.mIterator = iterator(self.mInfile, *args, **kwargs )

    def __iter__(self):
        return self

    def __del__(self):
        self.mInfile.close()

    def next( self ):

        if self.mInfile.tell() > self.mEndPos:
            raise StopIteration
        return self.mIterator.next()

class IteratorMultiline:

    def __init__(self, filename, nchunks, chunk, iterator, *args, **kwargs ):

        if filename.endswith(".gz"):
            self.mFileSize = getZipSize( filename )
            self.mInfile = gzip.open( filename, "r" )
        else:
            self.mInfile = open( filename, "r" )
            self.mInfile.seek(0, 2)
            self.mFileSize = self.mInfile.tell()

        if nchunks > 1:
            self.mChunkSize = int(math.ceil( float(self.mFileSize) / nchunks ) )
            start_pos = self.mChunkSize * chunk
            self.mEndPos = start_pos + self.mChunkSize
            # seek is expensive in compressed files
            self.mInfile.seek( start_pos + self.mChunkSize )
            self.mInfile.readline()
            self.mEndPos = self.mInfile.tell()

            self.mInfile.seek( start_pos )
            ## position yourself at a newline
            if start_pos > 0: self.mInfile.readline()
            self.mStartPos = self.mInfile.tell()
            self.mIterator = iterator(self.mInfile, *args, **kwargs )
            self.mLastPos = self.mStartPos
        else:
            self.mInfile.seek( 0 )
            self.mEndPos = self.mFileSize 
            self.mIterator = iterator(self.mInfile, *args, **kwargs )

        E.info( "nchunks=%i, chunk=%i, start=%i, filesize=%i, filename=%s" %\
                    (nchunks,  chunk, self.mInfile.tell(), self.mFileSize, filename ) )

    def __iter__(self):
        return self

    def __del__(self):
        self.mInfile.close()

    def next( self ):

        pos, record = self.mIterator.next()
        # print "next:", "start of record=", pos, record, "last=", self.mLastPos, "start=", self.mStartPos, "end=",self.mEndPos
        if pos > self.mEndPos: raise StopIteration
        return record

def iterator( infile ):

    if infile.tell() != 0: infile.readline()
    while 1:
        line = infile.readline()
        if not line: break
        yield line

def groupby( infile, key  ):

    # skip first putatively incomplete entry
    last_k, last_p, data = None, 0, []
    if infile.tell() != 0: 
        # read first line
        last_p = infile.tell()
        line = infile.readline()
        if not line: raise StopIteration
        last_k = key(line)
        # advance to next key change
        while 1:
            if key(line) != last_k:
                last_k = key(line)
                data.append( line )
                break
            last_p = infile.tell()
            line = infile.readline()
            if not line: break
    else:
        line = infile.readline() 
        # print "first", line,
        last_k = key( line )
        data.append( line )

    # iterate over full records
    while 1:
        p = infile.tell()
        line = infile.readline()
        # print "loop", "start=", p, "end=", infile.tell(), line,
        if not line: break
        k = key(line)
        if k != last_k:
            yield last_p, data
            # print "Reset", k, last_k
            last_k, last_p, data=k, p, []
        data.append(line)
        
    if data: yield last_p, data

class Iterator3:

    def __init__(self, filename, nchunks, chunk, iterator, *args, **kwargs ):

        if filename.endswith(".gz"):
            self.mFileSize = getZipSize( filename )
            self.mInfile = gzip.open( filename, "r" )
        else:
            self.mInfile = open( filename, "r" )
            self.mInfile.seek(0, 2)
            self.mFileSize = self.mInfile.tell()

        self.mChunkSize = int(math.ceil( float(self.mFileSize) / nchunks ) )
        start_pos = self.mChunkSize * chunk
        self.mInfile.seek( start_pos + self.mChunkSize )
        self.mInfile.readline()
        self.mEndPos = self.mInfile.tell()

        self.mInfile.seek( start_pos )
        ## position yourself at a newline
        if start_pos > 0: self.mInfile.readline()
        self.mStartPos = self.mInfile.tell()
        self.mIterator = iterator(self.mInfile, *args, **kwargs )
        self.mLastPos = self.mStartPos

    def __iter__(self):
        return self

    def __del__(self):
        self.mInfile.close()

    def next( self ):

        if self.mInfile.tell() > self.mEndPos:
            raise StopIteration
        return self.mIterator.next()

def groupby2( infile, key  ):

    # ignore first  
    print "start", infile.tell()
    x = itertools.groupby( infile, key )
    print "define", infile.tell()

    if infile.tell() > 0: 
        r = x.next()
        print "next", infile.tell()
        print "skipped=",r

    for k,i in x:
        print "before", infile.tell()
        l = list(i)
        print "after", infile.tell()
        yield l

