//--------------------------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "adda.h"


//--------------------------------------------------------------------------------------------
bool fileExists (const std::string & filename)
{
	FILE * infile = fopen(filename.c_str(), "r");
	if (infile != NULL)
	{
		fclose(infile);
		return true;
	}
	else
	{
		return false;
	}
}

//--------------------------------------------------------------------------------------------
FILE * openFileForWrite( const std::string & filename )
{
	if (fileExists(filename))
	{
		std::cerr << "# file " << filename << " already exists, aborting." << std::endl;
		exit(EXIT_FAILURE);
	}

	FILE * file = fopen( filename.c_str(), "w");

	if (file == NULL)
	{
		std::cerr << "# error while opening " << filename << " for writing." << std::endl;
		exit(EXIT_FAILURE);
	}

	return file;
}

//--------------------------------------------------------------------------------------------
FILE * openFileForRead( const std::string & filename )
{
	FILE * file = fopen( filename.c_str(), "r");
	if (file == NULL)
	{
		std::cerr << "# file " << filename << " could not be opened for reading." << std::endl;
		exit(EXIT_FAILURE);
	}

	return file;
}
