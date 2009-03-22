####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Experiment.py,v 1.1 2002/11/18 13:00:34 heger Exp $
##
##
####
####

"""Module for record keeping of experiments. The module
provides functions for

   * argument parsing
   * record keeping (logging)
   * benchmarking

"""

import string,re,sys,time,inspect,os,optparse,logging,random, uuid, collections, types, getopt

global_starting_time = time.time()
global_options = None
global_args    = None
global_id      = uuid.uuid4()
global_benchmark = collections.defaultdict( int )

def GetHeader():
    """return a header string with command line options and
    timestamp"""
    system, host, release, version, machine = os.uname()
    return "# output generated by %s\n# job started at %s on %s -- %s\n# pid: %i, system: %s %s %s %s" %\
           (string.join(sys.argv, " "),
            time.asctime(time.localtime(time.time())),
            host,
            global_id,
            os.getpid(),
            system, release, version, machine)

def GetParams( options = None):
    """return a string containing script parameters.
    Parameters are all variables that start with "param_".
    """
    result = []
    if options:
        members = options.__dict__
        for k, v in sorted(members.items()):
            result.append("# %-40s: %s" % (k, str(v).encode("string_escape")))
    else:
        vars = inspect.currentframe().f_back.f_locals
        for var in filter(lambda x: re.match("param_",x), vars.keys()):
            result.append("# %-40s: %s" % (var, str(vars[var]).encode("string_escape")))

    if result:
        return string.join(result, "\n")
    else:
        return "# no parameters."


def GetFooter():
    """return a header string with command line options and
    timestamp."""
    return "# job finished in %i seconds at %s -- %s -- %s" %\
           (time.time() - global_starting_time,
            time.asctime(time.localtime(time.time())),
            string.join( map( lambda x: "%5.2f" % x, os.times()[:4]), " "),
            global_id)

def Start( parser = None,
           quiet = False,
           add_csv_options = False,
           add_mysql_options = False,
           add_psql_options = False,
           add_pipe_options = True,
           add_cluster_options = False ):
    """set up an experiment.

    returns a tuple containing (options, args).

    The options class is extended with a logger module.
    """
    
    if not parser:
        parser = optparse.OptionParser( version = "%prog version: $Id$" )

    global global_options, global_args, global_starting_time

    global_starting_time = time.time()

    parser.add_option("-v", "--verbose", dest="loglevel", type="int",
                      help="loglevel [%default]. The higher, the more output." )

    parser.add_option( "--timeit", dest='timeit_file', type="string",
                       help="store timeing information in file [%default]." )        
    parser.add_option( "--timeit-name", dest='timeit_name', type="string",
                       help="name in timing file for this class of jobs [%default]." )
    parser.add_option( "--timeit-header", dest='timeit_header', action="store_true",
                       help="add header for timing information [%default]." )        
                                     
    if quiet:
        parser.set_defaults( loglevel = 0 )
    else:
        parser.set_defaults( loglevel = 1 )

    parser.set_defaults(
        timeit_file = None,
        timeit_name = 'all',
        timeit_header = None,
        )

    if add_csv_options:
        parser.add_option("--dialect", dest="csv_dialect", type="string",
                          help="csv dialect to use [%default]." )

        parser.set_defaults(
            csv_dialect = "excel-tab",
            csv_lineterminator = "\n",
            )

    if add_psql_options:
        parser.add_option("-C", "--connection", dest="psql_connection", type="string",
                          help="psql connection string [%default]."  )
        parser.add_option("-U", "--user", dest="user", type="string",
                          help="database user name [%default]."  )
        
        parser.set_defaults( psql_connection = "db:andreas" )
        parser.set_defaults( user = "" )
        
    if add_cluster_options:
        parser.add_option( "--use-cluster", dest="use_cluster", action = "store_true",
                          help="use cluster [%default]."  )
        parser.add_option( "--cluster-priority", dest="cluster_priority", type="int",
                          help="set job priority on cluster [%default]."  )
        parser.add_option( "--cluster-queue", dest="cluster_queue", type="string",
                          help="set cluster queue [%default]."  )
        parser.add_option( "--cluster-num-jobs", dest="cluster_num_jobs", type="int",
                          help="number of jobs to submit to the queue execute in parallel [%default]."  )
        parser.add_option( "--cluster-options", dest="cluster_options", type="string",
                          help="additional options for cluster jobs, passed on to qrsh [%default]."  )
        
        parser.set_defaults( use_cluster = False,
                             cluster_queue = "medium_jobs.q",
                             cluster_priority = -10,
                             cluster_num_jobs = 100,
                             cluster_options = "")

    if add_pipe_options:
        parser.add_option("-I", "--stdin", dest="stdin", type="string",
                          help="file to read stdin from [default = stdin].",
                          metavar = "FILE"  )
        parser.add_option("-L", "--log", dest="stdlog", type="string",
                          help="file with logging information [default = stdout].",
                          metavar = "FILE"  )
        parser.add_option("-E", "--error", dest="stderr", type="string",
                          help="file with error information [default = stderr].",
                          metavar = "FILE"  )
        parser.add_option("-S", "--stdout", dest="stdout", type="string",
                          help="file where output is to go [default = stdout].",
                          metavar = "FILE"  )
        
        parser.set_defaults( stderr = sys.stderr )
        parser.set_defaults( stdout = sys.stdout )
        parser.set_defaults( stdlog = sys.stdout )
        parser.set_defaults( stdin = sys.stdin )
        
    if add_mysql_options:
        if not parser.has_option( "--host"):
            parser.add_option("-H", "--host", dest="host", type="string",
                              help="mysql host [%default]."  )
        if not parser.has_option( "--database"):
            parser.add_option("-D", "--database", dest="database", type="string",
                              help="mysql database [%default]."  )
        if not parser.has_option( "--user"):
            parser.add_option("-U", "--user", dest="user", type="string",
                              help="mysql username [%default]."  )
        if not parser.has_option( "--password"):
            parser.add_option("-P", "--password", dest="password", type="string",
                              help="mysql password [%default]."  )
        if not parser.has_option( "--port"):
            parser.add_option("-O", "--port", dest="port", type="int",
                              help="mysql port [%default]."  )
        
        parser.set_defaults( host = "db",
                             port = 3306,
                             user = "",
                             password = "",
                             database = "" )

    (global_options, global_args) = parser.parse_args()

    if add_pipe_options:
        if global_options.stdout != sys.stdout: 
            global_options.stdout = open(global_options.stdout, "w")
        if global_options.stderr != sys.stderr:
            if global_options.stderr == "stderr":
                global_options.stderr = global_options.stderr
            else:
                global_options.stderr = open(global_options.stderr, "w")
        if global_options.stdlog != sys.stdout:
            global_options.stdlog = open(global_options.stdlog, "a")            
        if global_options.stdin != sys.stdin: 
            global_options.stdin = open(global_options.stdin, "r")
    else:
        global_options.stderr = sys.stderr
        global_options.stdout = sys.stdout
        global_options.stdlog = sys.stdout
        global_options.stdin = sys.stdin
    
    if global_options.loglevel >= 1:
        global_options.stdlog.write(GetHeader() + "\n" )
        global_options.stdlog.write(GetParams( global_options) + "\n")
        global_options.stdlog.flush()

    ## configure logging
    ## map from 0-10 to logging scale
    ## 0: quiet
    ## 1: little verbositiy
    ## >1: increased verbosity
    if global_options.loglevel == 0:
        lvl = logging.ERROR
    elif global_options.loglevel == 1:
        lvl = logging.INFO
    else:
        lvl = logging.DEBUG

    if global_options.stdout == global_options.stdlog:
        logging.basicConfig(
            level=lvl,
            format='# %(asctime)s %(name)s %(levelname)s %(message)s',
            stream = global_options.stdlog )
    else:
        logging.basicConfig(
            level=lvl,
            format='%(asctime)s %(name)s %(levelname)s %(message)s',
            stream = global_options.stdlog )
        
    return global_options, global_args
    
def Stop():
    """stop the experiment."""

    if global_options.loglevel >= 1 and global_benchmark:
        t = time.time() - global_starting_time
        global_options.stdlog.write( "######### Time spent in benchmarked functions ###################\n" )
        global_options.stdlog.write( "# function\tseconds\tpercent\n" )
        for key, value in global_benchmark.items():
            global_options.stdlog.write( "# %s\t%6i\t%5.2f%%\n" % (key, value, (100.0 * float(value) / t)))
        global_options.stdlog.write( "#################################################################\n" )

    if global_options.loglevel >= 1:
        global_options.stdlog.write(GetFooter() + "\n")

    if global_options.timeit_file:

        outfile = open(global_options.timeit_file, "a")

        if global_options.timeit_header:
            outfile.write( "\t".join( ("name", "wall", "user", "sys", "cuser", "csys",
                                       "host", "system", "release", "machine",
                                       "start", "end", "path", "cmd" ) ) + "\n" )

        csystem, host, release, version, machine = map(str, os.uname())
        uusr, usys, c_usr, c_sys = map( lambda x: "%5.2f" % x, os.times()[:4])
        t_end = time.time()
        c_wall = "%5.2f" % (t_end - global_starting_time)
        
        if sys.argv[0] == "run.py":
            cmd = global_args[0]
            if len(global_args) > 1:
                cmd += " '" + "' '".join(global_args[1:]) + "'" 
        else:
            cmd = sys.argv[0]

        result = "\t".join( (global_options.timeit_name,
                             c_wall, uusr, usys, c_usr, c_sys,
                             host, csystem, release, machine,
                             time.asctime(time.localtime( global_starting_time )),
                             time.asctime(time.localtime( t_end )),
                             os.path.abspath( os.getcwd() ),
                             cmd ) ) + "\n"

        outfile.write( result )
        outfile.close()


class Experiment2( optparse.Values ):
    """add an interface to the logging module.

    (work in progress).
    """

    def __init__(self, options ):

        for a, b in options.__dict__.items():
            setattr( self, a, b)
        
        logging.basicConfig(
            level=self.loglevel,
            format='%(asctime)s %(levelname)s %(message)s',
            stream = self.stdlog )

        #def __del__(self):
        #Stop()

    def log( self, message ):
        logging.log( message )

    def info( self, message ):
        logging.info( message )

    def warning( self, message ):
        logging.warning( message )

    def warn( self, message ):
        logging.warning( message )

    def debug( self, message ):
        logging.debug( message )

    def error( self, message ):
        logging.error( message )
        
    def critical( self, message):
        logging.critical( message )

class Experiment:

    mShortOptions = ""
    mLongOptions  = []

    mLogLevel = 0
    mTest = 0
    mDatabaseName = None

    mName = sys.argv[0]


    def __init__(self ):

        # process command-line arguments
        (self.mOptlist, self.mArgs) = self.ParseCommandLine()

        # set options now
        self.ProcessOptions(self.mOptlist)

    def DumpParameters( self ):
        """dump parameters of this object. All parameters start with a lower-case m."""

        members = self.__dict__

	print "#--------------------------------------------------------------------------------------------"
        print "#" + string.join(sys.argv)
        print "# pid: %i, system:" % os.getpid(), string.join(os.uname(), ",") 
	print "#--------------------------------------------------------------------------------------------"
        print "# Parameters for instance of <" + self.mName + "> on " + time.asctime(time.localtime(time.time()))

        member_keys = list(members.keys())
        member_keys.sort()
        for member in member_keys:
            if member[0] == 'm':
                print "# %-40s:" % member, members[member]

	print "#--------------------------------------------------------------------------------------------"
        sys.stdout.flush()

    #-----------------------------> Control functions <--------------------------------

    ##------------------------------------------------------------------------------------
    def ProcessOptions( self, optlist ):
        """Sets options in this module. Please overload as necessary."""

        for o,a in optlist:
            if o in ( "-V", "--Verbose" ):
                self.mLogLevel = string.atoi(a)
            elif o in ( "-T", "--test") :
                self.mTest = 1

    ##------------------------------------------------------------------------------------        
    def ProcessArguments( self, args ):
        """Perform actions as given in command line arguments."""

        if self.mLogLevel >= 1:
            self.DumpParameters()
            
        for arg in args:
            if arg[-1] == ")":
                statement = "self.%s" % arg
            else:
                statement = "self.%s()" % arg
            exec statement

            if self.mLogLevel >= 1:
                print "--------------------------------------------------------------------------------------------"
                print statement + " finished at " + time.asctime(time.localtime(time.time()))
                print "--------------------------------------------------------------------------------------------"

    ##------------------------------------------------------------------------------------            
    def ParseCommandLine( self ):
        """Call subroutine with command line arguments."""

        self.mShortOptions = self.mShortOptions + "V:D:T"
        self.mLongOptions.append( "Verbose=" )
        self.mLongOptions.append( "Database=" )
        self.mLongOptions.append( "Test" )

        try:
            optlist, args = getopt.getopt(sys.argv[1:],
                                          self.mShortOptions,
                                          self.mLongOptions)
        except getopt.error, msg:
            self.PrintUsage()
            print msg
            sys.exit(2)

        return optlist, args

    #--------------------------------------------------------------------------------        
    def Process( self ):
        self.ProcessArguments(self.mArgs)

    #--------------------------------------------------------------------------------
    def PrintUsage( self ):
        """print usage information."""

        print "# valid short options are:", self.mShortOptions
        print "# valid long options are:", str(self.mLongOptions)
    
def benchmark(func):
    """Decorator adding printing of time taken to run methods."""

    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        key = "%s:%i" % (func.func_name, func.func_code.co_firstlineno)
        global_benchmark[key] += t2-t1
        global_options.stdlog.write( '## benchmark: %s completed in %6.4f s\n' % (key, (t2-t1)))
        global_options.stdlog.flush()
        return res
    return wrapper

# there are differences whether you cache a function or
# an objects method
def cachedmethod(function):
    return Memoize(function)

class Memoize(object):
    def __init__(self, fn):
        self.cache={}
        self.fn=fn

    def __get__(self, instance, cls=None):
        self.instance = instance
        return self

    def __call__(self,*args):
        if self.cache.has_key(args):
            return self.cache[args]
        else:
            object = self.cache[args] = self.fn(self.instance, *args)
            return object

def log( loglevel, message ):
    """log message at loglevel."""
    logging.log( loglevel, message )

def info( message ):
    logging.info( message )

def warning( message ):
    logging.warning( message )

def warn( message ):
    logging.warning( message )

def debug( message ):
    logging.debug( message )

def error( message ):
    logging.error( message )
        
def critical( message):
    logging.critical( message )

def getLogLevel():
    return global_options.loglevel
