import math

class Iterator:

    def __init__(self, filename, nchunks, chunk, iterator, *args ):
        
        self.mInfile = open( filename, "r" )
        
        self.mFileSize = self.mInfile.seek( 2, 0) 
        self.mInfile.seek(0, 2)
        self.mFileSize = self.mInfile.tell()
        
        chunk_size = int(math.ceil( float(self.mFileSize) / nchunks ) )
        start_pos = chunk_size * chunk
        self.mEndPos = start_pos + chunk_size

        self.mInfile.seek( start_pos )
        self.mIterator = iterator(self.mInfile, *args )

    def __iter__(self):
        return self

    def next( self ):

        # print "tell", self.mInfile.tell(), self.mEndPos, self.mFileSize
        
        if self.mInfile.tell() > self.mEndPos:
            raise StopIteration
        return self.mIterator.next()
            
def iterator( infile ):

    # print "starting iterator at", infile.tell()

    if infile.tell() != 0: infile.readline()
    while 1:
        line = infile.readline()
        if not line: break
        yield line

def iterator_group_by( infile, num_columns, group_by ):

    l, d = "#", []
    while l.startswith("#") and len(d) != num_columns:
        l = infile.readline()
        d = l.split( "\t" )

    while 1:
        yield d
        l = infile.readline()
        d = l.split( "\t" )
