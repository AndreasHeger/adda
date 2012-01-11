####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Pairsdb.py,v 1.3 2002/11/18 13:03:28 heger Exp $
##
##
####
####


#----------------------------------------------------------------
# Name:		Pairsdb
#--------------General Information-------------------------------
# File:		Pairsdb.py
# Version:      1.0
# Description:	basic definitions used by all python modules
# Author:	Andreas Heger (heger@ebi.ac.uk)
#--------------Documentation-------------------------------------
#
#
#--------------Change History------------------------------------
# 20 Jan 2000	Created
# 
#
#----------------------------------------------------------------

import sys
import os
from Database import Database

# in the near future make this file read a configuration file.

# note, that symbolic links break at a few cases. In particular, the routines os.path can not
# follow links across different partitions?

# Hierarchy of loglevels, uncomment only one
#LOGLEVEL       = 0		# completely silent output
#LOGLEVEL       = 1		# write process-indicators
#LOGLEVEL       = 2		# write information
#LOGLEVEL       = 3      # dump out every SQL-statement

LOGLEVEL       = 3		# completely silent output

PATH_HOME = "/production"

PATH_TEMP= "/production/tmp"

DEFAULT_SOCKET          = None # "/tmp/mysql.sock"
DEFAULT_HOST            = 'localhost'
DEFAULT_PORT            = 3306

# path to where temporary files are stored
PATH_LOAD       = '/tmp'
PATH_DUMP       = '/tmp'                 

# path to helper executables, platform dependent
PATH_EXEC      = '/data/bin'

# path to scripts
PATH_SCRIPTS   = '/data/scripts/pairsdb'

PATH_DATA      = '/data/data'                             # path to data used by helpers

PATH_BLAST     = '/data/databases/blast'			# path, where to put blast-related data-stuff

PATH_LOCK      = PATH_HOME + '/lock'                            # path, where all the locks reside
PATH_SHARE     = PATH_HOME + '/share'                           # path to shared directory

DEFINITIONFILE = PATH_HOME + '/share/tables.txt'
USERSFILE      = PATH_HOME + '/share/users.txt'

# Environment fuer Pairsdb
FILE_LOG                = PATH_HOME + '/log/logfile'

REFERENCE_DATABASE_SIZE = 65000000                                               # for normalizing E-Values
EVALUE_CUTOFF           = 1.0

#----------------------------------------
# Function:    Connect
# Description: Connect to mysql-database
# Author     : Andreas Heger
# Created    : 22.1.2000
#----------------------------------------
class Pairsdb (Database):

    def __init__( self ):
        Database.__init__( self )
        
    def Connect( self,
		 host = DEFAULT_HOST,
		 user = 'heger',
		 passwd = 'HeinBloed',
		 unix_socket = DEFAULT_SOCKET,
                 port = DEFAULT_PORT,
                 dbname = 'pairsdb' ):
        """connect to a database on in mysql.
        """
        
	self.dbname = dbname
	self.host   = host
	self.user   = user
	self.passwd = passwd
	self.socket = unix_socket
        self.port   = port
	
	return Database.Connect( self, host = host, user = user, passwd = passwd, dbname = dbname, port = port, socket = unix_socket  )

    def Create( self,
		host = 'localhost',
		user = 'root',
		passwd = 'HeinBloed',
		unix_socket = DEFAULT_SOCKET ):

	self.host   = host
	self.user   = user
	self.passwd = passwd
	self.socket = unix_socket
	self.dbname = 'pairsdb'

	return Database.Create( self )





