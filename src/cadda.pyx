"""
Pyrex extension classes used by `cadda.py`.
"""

cdef extern from "string.h":
    ctypedef int size_t
    void *memcpy(void *dst,void *src,size_t len)
    void *memmove(void *dst,void *src,size_t len)
    void *memset(void *b,int c,size_t len)
    size_t strlen(char *s)
    char *strncpy(char *dest, char *src, size_t n)

cdef extern from "stdlib.h":
    void free(void *)
    void *malloc(size_t)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    int c_abs "abs" (int)
    void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "stdio.h":
    ctypedef struct FILE:
        pass
    ctypedef struct fpos_t:
        pass
    FILE *fopen(char *,char *)
    int fclose(FILE *)
    int sscanf(char *str,char *fmt,...)
    int sprintf(char *str,char *fmt,...)
    int fprintf(FILE *ifile,char *fmt,...)
    int ferror(FILE *stream)
    size_t fwrite( void *ptr, 
                   size_t size, 
                   size_t nmemb,
                   FILE *stream)
    size_t fread(void *ptr, 
                 size_t size, 
                 size_t nmemb, 
                 FILE *stream)
    char *fgets(char *str,int size,FILE *ifile)
    int fgetpos(FILE *stream, fpos_t *pos)
    int fsetpos(FILE *stream, fpos_t *pos)
    int printf(char *format, ...)

cdef extern from "string.h":
    int strcmp(char *s1, char *s2)
    int strncmp(char *s1,char *s2,size_t len)
    char *strcpy(char *dest,char *src)
    char *strdup(char *)
    char *strcat(char *,char *)

#cdef extern from "zlib.h":
#    ctypedef struct z_stream:
#        unsigned char * next_in
#        int avail_in
#        void * next_out
#        int avail_out
#    int deflate(z_stream * strm, int flush)
#    int deflateEnd( z_stream * strm )
#    int deflateInit( z_stream * strm, int level )

cdef extern from "cadda.h":
    int cadda_optimise_initialise()
    int cadda_optimise_destroy()
    double cadda_optimise_iteration()
    int cadda_optimise_save_partitions( char * )
    int cadda_optimise_load_partitions( char * )
    long cadda_optimise_get_num_partitions()
    int cadda_convert( char * )
    long cadda_build_mst( char *, char * )
    int cadda_build_index()
    int cadda_check_index()
    void cadda_dump_parameters()
    void cadda_setFilenameSegments( char *)
    void cadda_setFilenameGraph( char *)
    void cadda_setFilenameIndex( char *)
    void cadda_setFilenameTransfers( char *)
    void cadda_setFilenameNids( char *)
    void cadda_setFilenameDomains( char *)
    void cadda_setFilenameDomainGraph( char *)
    void cadda_setFilenameMst( char *)
    void cadda_setLogLevel(int)
    void cadda_setResolution( int )
    void cadda_setReportStep( int )
    void cadda_setK(double)
    void cadda_setC(double)
    void cadda_setMax(double)
    void cadda_setMin(double)
    void cadda_setE(double)
    void cadda_setF(double)
    void cadda_setRelativeOverhang( int )
    void cadda_setOnlyQuery( int )
    void cadda_setDescend( int )
    void cadda_setDisallowShortening( int )
    void cadda_setMaxIterations( int )
    void cadda_setEvalueThresholdTrustedLinks( double) 
    int toCompressedFile( unsigned char *, size_t, FILE * )
    int fromCompressedFile( unsigned char *, size_t, FILE * )
    
def optimise_iteration():
    return cadda_optimise_iteration()

def optimise_initialise():
    return cadda_optimise_initialise()

def optimise_destroy():
    return cadda_optimise_destroy()

def optimise_get_num_partitions():
    return cadda_optimise_get_num_partitions()

def optimise_load_partitions( filename ):
    return cadda_optimise_load_partitions( filename )

def optimise_save_partitions( filename ):
    return cadda_optimise_save_partitions( filename )

def convert( filename):
    return cadda_convert( filename )

def build_mst( out_filename, in_filename ):
    return cadda_build_mst( out_filename, in_filename )

def build_index():
    return cadda_build_index()

def check_index():
    return cadda_check_index()

def dump_parameters():
    cadda_dump_parameters()

def setFilenameSegments(v):
    """set input filename with segments."""
    cadda_setFilenameSegments(v)

def setFilenameGraph(v):
    """set input filename with graph."""
    cadda_setFilenameGraph(v)

def setFilenameIndex(v):
    """set input filename with index."""
    cadda_setFilenameIndex(v)

def setFilenameMst(v):
    """set filename with mst."""
    cadda_setFilenameMst(v)

def setFilenameNids(v):
    """set input filename with nids."""
    cadda_setFilenameNids(v)

def setFilenameDomains(v):
    """set input filename with domains."""
    cadda_setFilenameDomains(v)

def setFilenameDomainGraph(v):
    """set filename of domain graph."""
    cadda_setFilenameDomainGraph(v)

def setFilenameTransfers(v):
    """set input filename with transfer data."""
    cadda_setFilenameTransfers(v)

def setLogLevel(v):
    """set the logging level."""
    cadda_setLogLevel(v)

def setSigmoidK(v):
    """set parameter K: smoothness of sigmoid."""
    cadda_setK(v)

def setSigmoidC(v):
    """set parameter C: inflection point of sigmoid."""
    cadda_setC(v)

def setSigmoidMax(v):
    """set parameter M: maximum sigmoid."""
    cadda_setMax(v)

def setSigmoidMin(v):
    """set parameter N: minimum sigmoid."""
    cadda_setMin(v)

def setExponentialE(v):
    """set parameter E: exponential decay rate."""
    cadda_setE(v)

def setExponentialF(v):
    """set parameter F: exponential decay scale."""
    cadda_setF(v)

def setRelativeOverhang(v):  
    """if true, use relative overhang."""
    cadda_setRelativeOverhang( v )

def setOnlyQuery(v):  
    """if true, use only the score based on the query."""
    cadda_setOnlyQuery( v )

def setResolution(v):  
    """resolution to use."""
    cadda_setResolution( v )

def setDescend(v):  
    """if true, descend."""
    cadda_setDescend( v )

def setDisallowShortening(v):  
    """if true, disallow shortening."""
    cadda_setDisallowShortening( v )

def setMaxIterations(v):  
    """set the maximum number of iterations."""
    cadda_setMaxIterations( v )

def setReportStep(v):
    """set reporting interval."""
    cadda_setReportStep(v)
    
def setEvalueThresholdTrustedLinks( v ):
    """set evalue threshold for trusted links.""" 
    cadda_setEvalueThresholdTrustedLinks( v )
    
import alignlib

DEF MAX_BUFFER_SIZE = 10000000

cdef extern from "adda.h":

    ctypedef int Nid
    ctypedef int Residue
    ctypedef fpos_t FileIndex
    ctypedef int Length
    ctypedef int uResidue

## todo: convert to class
ctypedef struct Neighbour:
    Nid sbjct_nid
    float evalue
    uResidue query_start
    uResidue query_end
    uResidue sbjct_start
    uResidue sbjct_end
    Length query_alen
    Length sbjct_alen
    char * query_ali
    char * sbjct_ali
    
cdef fromNeighbour( Neighbour * n, neighbour ):
    '''load data from neighbour.'''
    if n.query_ali != NULL: free( n.query_ali)
    if n.sbjct_ali != NULL: free( n.sbjct_ali)
    n.sbjct_nid = neighbour.mSbjctToken
    n.query_start = neighbour.mQueryFrom
    n.query_end = neighbour.mQueryTo
    n.sbjct_start = neighbour.mSbjctFrom
    n.sbjct_end = neighbour.mSbjctTo
    n.evalue = neighbour.mEvalue
    n.query_alen = len( neighbour.mQueryAli )
    n.sbjct_alen = len( neighbour.mSbjctAli )

    # copy string explicitely, as lifetime of python object neighbours
    # is no guaranteed.
    n.query_ali = <char*>calloc( n.query_alen + 1, sizeof(char) )
    n.sbjct_ali = <char*>calloc( n.sbjct_alen + 1, sizeof(char) )

    strncpy( n.query_ali, neighbour.mQueryAli, n.query_alen + 1)
    strncpy( n.sbjct_ali, neighbour.mSbjctAli, n.sbjct_alen + 1)

class NeighbourRecord(object):

    def __str__( self ):

        return "\t".join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.mEvalue,
            self.mQueryFrom, self.mQueryTo,
            self.mSbjctFrom, self.mSbjctTo )))

cdef toNeighbour( Nid query_nid, Neighbour * n ):
    '''load data from neighbour.'''

    dest = NeighbourRecord()

    dest.mQueryToken = query_nid
    dest.mSbjctToken = n.sbjct_nid
    dest.mEvalue = n.evalue
    dest.mQueryFrom = n.query_start
    dest.mQueryTo = n.query_end
    dest.mQueryAli = n.query_ali
    dest.mSbjctFrom = n.sbjct_start
    dest.mSbjctTo = n.sbjct_end
    dest.mSbjctAli = n.sbjct_ali
    
    return dest

cdef toFile( Neighbour * n, FILE * output_f ):
    '''write neighbour to file'''
    fwrite( &n.sbjct_nid, sizeof(Nid), 1, output_f )
    fwrite( &n.evalue, sizeof(float), 1, output_f )
    fwrite( &n.query_start, sizeof(uResidue), 1, output_f )
    fwrite( &n.query_end, sizeof(uResidue), 1, output_f )
    fwrite( &n.sbjct_start, sizeof(uResidue), 1, output_f )
    fwrite( &n.sbjct_end, sizeof(uResidue), 1, output_f )
    fwrite( &n.query_alen, sizeof(Length), 1, output_f )
    fwrite( &n.sbjct_alen, sizeof(Length), 1, output_f )
    fwrite( n.query_ali, sizeof( char ), n.query_alen, output_f )
    fwrite( n.sbjct_ali, sizeof( char ), n.sbjct_alen, output_f )

cdef fromFile( Neighbour * n, FILE * input_f ):
    '''read neighbour from file'''
    if n.query_ali != NULL: free( n.query_ali)
    if n.sbjct_ali != NULL: free( n.sbjct_ali)
    fread( &n.sbjct_nid, sizeof(Nid), 1, input_f )
    fread( &n.evalue, sizeof(float), 1, input_f )
    fread( &n.query_start, sizeof(uResidue), 1, input_f )
    fread( &n.query_end, sizeof(uResidue), 1, input_f )
    fread( &n.sbjct_start, sizeof(uResidue), 1, input_f )
    fread( &n.sbjct_end, sizeof(uResidue), 1, input_f )
    fread( &n.query_alen, sizeof(Length), 1, input_f )
    fread( &n.sbjct_alen, sizeof(Length), 1, input_f )
    n.query_ali = <char*>calloc( n.query_alen + 1, sizeof(char) )
    n.sbjct_ali = <char*>calloc( n.sbjct_alen + 1, sizeof(char) )
    fread( n.query_ali, sizeof( char ), n.query_alen, input_f )
    fread( n.sbjct_ali, sizeof( char ), n.sbjct_alen, input_f )

cdef toStdout( Nid query_nid, Neighbour * n ):
    '''print neighbour to stdout in pairsdb format'''
    printf("%i\t%i\t%f\t%i\t%i\t%s\t%i\t%i\t%s\n",
           query_nid,
           n.sbjct_nid,
           n.evalue,
           n.query_start,
           n.query_end,
           n.query_ali,
           n.sbjct_start,
           n.sbjct_end,
           n.sbjct_ali)

cdef unsigned char * toBuffer( Neighbour * n, unsigned char * buffer):
    '''copy information in *n* into buffer.

    returns pointer to position in buffer after writing all data
    '''
    cdef size_t s
    s = sizeof( Neighbour ) - 2 * sizeof( char * )
    memcpy( buffer, n, s )
    buffer += s

    s = sizeof( char ) * (n.query_alen + 1)
    memcpy( buffer, n.query_ali, s )
    buffer += s

    s = sizeof( char ) * (n.sbjct_alen + 1)
    memcpy( buffer, n.sbjct_ali, s )
    buffer += s

    return buffer

cdef unsigned char * fromBuffer( Neighbour * n, unsigned char * buffer):
    '''copy data in *n* into buffer.

    returns pointer to position in buffer after reading all data
    '''
    if n.query_ali != NULL: free( n.query_ali)
    if n.sbjct_ali != NULL: free( n.sbjct_ali)

    cdef size_t s
    s = sizeof( Neighbour ) - 2 * sizeof( char * )
    memcpy( n, buffer, s )
    buffer += s

    n.query_ali = <char*>calloc( n.query_alen + 1, sizeof(char) )
    n.sbjct_ali = <char*>calloc( n.sbjct_alen + 1, sizeof(char) )

    s = sizeof( char ) * ( n.query_alen + 1)
    memcpy( n.query_ali, buffer, s )
    buffer += s

    s = sizeof( char ) * (n.sbjct_alen + 1)
    memcpy( n.sbjct_ali, buffer, s )
    buffer += s

    return buffer

DEF Z_OK           = 0
DEF Z_STREAM_END   = 1
DEF Z_NEED_DICT    = 2
DEF Z_ERRNO        = (-1)
DEF Z_STREAM_ERROR = (-2)
DEF Z_DATA_ERROR   = (-3)
DEF Z_MEM_ERROR    = (-4)
DEF Z_BUF_ERROR    = (-5)
DEF Z_VERSION_ERROR = (-6)

# cdef abcToCompressed( unsigned char * buffer, size_t size, FILE * output_f ):
#     '''compress buffer of size.'''


#     cdef z_stream strm
#     cdef int level = 9
#     cdef size_t have

#     print "initing"
#     ret = deflateInit(&strm, level);
#     if ret != Z_OK: return ret
                                                                                                                                                                                   
#     cdef unsigned char * compressed
#     cdef size_t max_compressed
#     max_compressed = size * 2    
#     print "allocaitng"
    
#     compressed = <unsigned char *>calloc( max_compressed, sizeof(unsigned char) )

#     # compress in one go
#     strm.next_in = <unsigned char*>buffer
#     strm.avail_in = size
#     strm.avail_out = max_compressed
#     strm.next_out = compressed

#     print "before deflating"
#     ret = deflate(&strm, 1)
    
#     assert ret != Z_STREAM_ERROR, "compression error"

#     have = max_compressed - strm.avail_out

#     if (fwrite(compressed, 1, have, output_f) != have or ferror(output_f)):
#         deflateEnd(&strm)
#         free( compressed )
#         return Z_ERRNO

#     deflateEnd(&strm)
#     free( compressed )
#     return Z_OK
    
def indexGraph( graph_iterator, num_nids, output_filename_graph, output_filename_index ):
    """translate the pairsdb input graph into an ADDA formatted graph.
    
    This method reformats and indexes a neighbourhood graph.

    The number of nids must be known beforehand and the nids are assumed
    to be contiguous from 1 to num_nids.

    The ADDA graph format is binary and consists of records of
    neighbourhood lists. Each record starts with:

    Nid query_nid
    size_t nneighbours
    Neigbour [] neighbours
    
    where Neighbour is a struct of:

    Nid sbjct_nid 
    float evalue
    uResidue query_start
    uResidue query_end
    uResidues bjct_start
    uResidue sbjct_end
    Length query_ali_len
    Length sbjct_ali_lon
    char [] query_ali  
    char [] sbjc_ali
    
    query_ali and sbjct_ali are `\0` terminated strings.
    
    Each neighbour-record is gzipped.

    The index format is:
    Nid number of nids
    FileIndex [] index
    """

    # allocate index
    cdef FileIndex * index
    cdef Nid nnids 
    # add 1 for nid=0
    nnids = num_nids + 1
    # sets file positions for unknown ids to 0.
    index = <FileIndex*>calloc( nnids, sizeof( FileIndex ) )
    if index == NULL:
        raise ValueError( "memory allocation for index failed" )

    # open output file
    cdef FILE * output_f
    
    output_f = fopen( output_filename_graph, "wb" );
    if output_f == NULL:
        raise ValueError( "opening of file %s failed" % output_filename_graph )

    # iterate over graph
    cdef Nid query_nid
    cdef FileIndex pos
    cdef Neighbour neighbour
    neighbour.query_ali = NULL
    neighbour.sbjct_ali = NULL
    cdef size_t nneighbours
    cdef unsigned char * buffer
    buffer = <unsigned char *>calloc( MAX_BUFFER_SIZE, sizeof(unsigned char) )
    cdef unsigned char * p1 
    cdef size_t used

    # write empty entry with nid 0. This is a place-holder
    # for entries without neighbours
    nneighbours = 0
    query_nid = 0
    fgetpos( output_f, &pos )
    index[query_nid] = pos
    fwrite( &query_nid, sizeof( Nid ), 1, output_f )
    fwrite( &nneighbours, sizeof( size_t ), 1, output_f )

    for neighbours in graph_iterator:
        
        if neighbours == None: break

        query_nid = neighbours.mQueryToken

        # save index position
        fgetpos( output_f, &pos )
        index[query_nid] = pos

        # convert neighbours
        nneighbours = len(neighbours.mMatches)
        
        # write record to file
        fwrite( &query_nid, sizeof( Nid ), 1, output_f )
        fwrite( &nneighbours, sizeof( size_t), 1, output_f )

        p1 = buffer
        for n in neighbours.mMatches:
            fromNeighbour( &neighbour, n )
            p1 = toBuffer( &neighbour, p1 )
            
        used = p1 - buffer
        err = toCompressedFile( buffer, used, output_f )
        if err: raise ValueError( "error %i while writing compressed buffer to file" )

    fclose( output_f )

    # save index
    output_f = fopen( output_filename_index, "wb" );
    if output_f == NULL:
        raise ValueError( "opening of file %s failed" % output_filename_index )
    fwrite( &nnids, sizeof( Nid ), 1, output_f )
    fwrite( index, sizeof( FileIndex ), nnids, output_f )
    fclose(output_f)


cdef class IndexedNeighbours:
    """access to indexed ADDA graph."""

    cdef FILE * mFile
    cdef FileIndex * mIndex
    cdef Nid mNids

    def __init__(self, filename_graph, filename_index ):

        cdef FILE * index_f
        index_f = fopen( filename_index, "rb" )
        if index_f == NULL:
            raise OSError( "could not open index %s " % filename_index )
        cdef Nid nnids
        if fread( &nnids, sizeof(Nid), 1, index_f ) != 1 or ferror( index_f):
            raise OSError( "could not read index from %s" % filename_index )

        self.mIndex = <FileIndex *>calloc( sizeof(FileIndex), nnids)
        if self.mIndex == NULL:
            raise MemoryError( "out of memory when allocating index for %i nids" % nnids )

        if fread( self.mIndex, sizeof(FileIndex), nnids, index_f ) != nnids or ferror(index_f):
            raise OSError( "could not read index from %s" % filename_index )
        
        self.mFile = fopen( filename_graph, "rb" )
        
        self.mNids = nnids

    def getNeighbours( self, nid ):
        '''retrieve neighbours for *nid*'''
        
        assert 0 < nid < self.mNids, "nid %i out of range, maximum is %i" % (nid, self.mNids - 1)
        
        cdef int r
        r = fsetpos( self.mFile, &self.mIndex[nid] )

        if r != 0:
            raise OSError( "Could not go to file position for nid %i" % nid )
        
        cdef int n 
        cdef size_t nneighbours
        cdef Nid query_nid 

        n = fread( &query_nid, sizeof(Nid), 1, self.mFile )
        n += fread( &nneighbours, sizeof(size_t), 1, self.mFile )
    
        assert n == 2, "wrong item count while reading from graph"
        if nid != query_nid and query_nid != 0:
            raise ValueError( "index returned wrong nid: %i instead of %i" % (query_nid, nid) )
        
        if query_nid == 0: return []

        cdef unsigned char * buffer
        buffer = <unsigned char *>calloc( MAX_BUFFER_SIZE, sizeof(unsigned char) )

        cdef int retval
        retval = fromCompressedFile( buffer, MAX_BUFFER_SIZE, self.mFile )
        if retval != 0: raise ValueError("error while reading data for %i" % nid )

        # create neighbours
        cdef Neighbour neighbour
        neighbour.query_ali = NULL
        neighbour.sbjct_ali = NULL

        cdef unsigned char * p
        cdef int i

        result = []

        p = buffer
        for i from 0 <= i < nneighbours:
            p = fromBuffer( &neighbour, p )
            result.append( toNeighbour( query_nid, &neighbour) )

        free( buffer )
        return result

