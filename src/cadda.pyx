"""
Pyrex extension classes used by `cadda.py`.
"""

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
    cadda_setQueryOverhang( v )

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
     
