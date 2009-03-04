import sys, os, re, time, math, copy, glob, optparse, gzip, types, subprocess
import hashlib, base64, string
import fcntl
import logging

class Error(Exception):
    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)


class AddaError(Error):
    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line = None):
        if line:
            self.message = message + " at line " + line
        else:
            self.message = message
        
class AddaModule:
    
    mName = "Adda"
    
    def __init__(self, options, config, fasta = None):

        self.mOptions = options
        self.mConfig = config
        self.mFasta = fasta
        self.mRequirements = []
        
        self.mInput = 0
        self.mOutput = 0
        self.mTime = 0
        
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename='adda.log',
                            filemode='a')
        
        self.mLogger = logging.getLogger( self.mName )
        
        self.mReportStep = self.mConfig.get( "adda", "report_step", 1000 )
        self.mFilenamePersistence = self.mConfig.get("adda", "filename_persistence", "adda.private" )
        
    def apply(self, neighbours ):
        """perform action on neighbours."""

        self.debug( "started: token=%s" % (neighbours.mQueryToken ))

        self.mInput += 1
        t1 = time.time()
        self.applyMethod( neighbours )
        t2 = time.time()
        self.mTime += t2 - t1
        self.mOutput += 1

        self.debug( "finished: token=%s, time=%i" % (neighbours.mQueryToken,
                                                     t2-t1 ) )
    
    def run(self):
        self.mInput = 0
        self.mOutput = 0
        t1 = time.time()
        self.applyMethod()
        t2 = time.time()
        self.mTime += t2 - t1

    def execute( self, cmd ):
        """execute a shell command."""

        self.info( "executing command: %s" % cmd )
        try:
            retcode = subprocess.call( cmd, shell=True)
            if retcode < 0:
                self.error( "command failed with code %i: %s" % ( -retcode, cmd ) )
                raise AddaError( "child command %s failed, see log for details." % cmd )
        except OSError, e:
            self.error( "command failed: %s: %s" % ( -retcode, str(e), cmd ) )
            raise AddaError( "child command %s failed, see log for details." % cmd )

    def finish(self):
        """do aggregate computations (if needed)."""

        self.info( "finished: ninput=%i, noutput=%i, time=%i" % (self.mInput, 
                                                                 self.mOutput,
                                                                 self.mTime ) )
    
        
    def getHID ( self, sequence ):
        """returns a hash identifier for a sequence.
        """
        # do the encryption
        h = hashlib.md5(sequence).digest()
        # map to printable letters: hid has length 22
        r = base64.encodestring(h)[0:22]

        # finally substitute some characters:
        # '/' for '_', so we have legal file names
        # '[' for '+' and ']' for '=' for internet-applications
        hid = string.replace(r  , '/', '_')
        hid = string.replace(hid, '+', '[')
        hid = string.replace(hid, '=', ']')
        return hid

    def registerExistingOutput(self, filename ):
        """process existing output in infile to guess correct point to continue computation."""
        raise AddaError( "module does not support appending.""")
    
    def getLastLine( self, filename, read_size = 1024 ):
        """return last line of a file.
        """
        f = open(filename, 'rU' )  # U is to open it with Universal newline support
        offset = read_size
        f.seek(0, 2)
        file_size = f.tell()
        if file_size == 0:
            f.close() 
            return None
        while 1:
            if file_size < offset:
                offset = file_size
            f.seek(-1*offset, 2)
            read_str = f.read(offset)
            # Remove newline at the end
            if read_str[offset - 1] == '\n':
                read_str = read_str[:-1]
            lines = read_str.split('\n')
            if len(lines) >= 2:
                f.close()
                return lines[-1]
            if offset == file_size:   # reached the beginning
                f.close()
                return read_str
            offset += read_size
        f.close()
        return None
    
    def openOutputStream(self, filename ):
        """opens an output stream.
        
        If the output filename exists an error is raised
        unless
        1. --force is set: the existing file will be overwritten
        2. --append is set: data will be appended. The registerExistingOutput
            method is called to give the module the chance to advance 
            the input stream to the appropriate point for continuation.
        """

        if self.mOptions.suffix:
            fn = filename + self.mOptions.suffix
        else:
            fn = filename

        
        if os.path.exists( fn ):
            if self.mOptions.append:            
                self.registerExistingOutput( fn )
                self.log( "opening output file %s in append mode" % fn )
                return open( fn, "a" )
            elif self.mOptions.force:
                self.log( "opening output file %s: overwriting existing file" % fn )
                return open( fn, "w" )
            else:
                raise AddaError( "output file %s already exists" % fn )

        self.debug( "opening output file %s" % fn )

        return open( fn, "w" )
        
    def log(self, msg, loglevel = 1 ):
        self.mLogger.log( loglevel, msg )
        
    def warn(self, msg ):
        self.mLogger.warn( msg )

    def debug(self, msg ):
        self.mLogger.debug( msg )

    def info( self, msg ):
        self.mLogger.info( msg )

    def critical( self, msg ):
        self.mLogger.critical( msg )

    def error( self, msg ):
        self.mLogger.error( msg )
        
    def saveValue(self, key, value ):
        """save a key/value pair in the persistence file."""
        fd = os.open( self.mFilenamePersistence, "a" )
        fcntl.flock( fd, fcntl.LOCK_EX )
        v = PyConfigParser.read( fd )
        if not v.has_section("data"):
            v.add_section("data")
        v.set( "data", key, value )
        fd.seek( 0 )
        PyConfigParser.write( fd )
        fcntl.flock( fd, fcntl.LOCK_UN)
        os.close( fd )
        
    def getValue(self, key):
        """get a value for a given key in the persistence file."""
        fd = os.open( self.mFilenamePersistence, "r" )
        fcntl.flock( fd, fcntl.LOCK_SH )
        v = PyConfigParser.read( fd )
        k = v.get( "data", key, value )
        fcntl.flock( fd, fcntl.LOCK_UN)
        os.close( fd )
        
        
        
