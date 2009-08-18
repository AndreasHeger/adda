//--------------------------------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <zlib.h>
#include "adda.h"
#include <cassert>

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

// save buffer into file with zlib compression
int toCompressedFile( unsigned char * buffer, size_t uncompressed_size, FILE * output_f )
{

  uLongf compressed_size = uncompressed_size * 2;
  Bytef * compressed = (Bytef *)calloc( compressed_size, sizeof(Bytef) );
  int level = 9;
  int zok = compress2(compressed, &compressed_size, buffer, uncompressed_size, level);

  if ( zok != Z_OK || fwrite( &compressed_size, sizeof(uLongf), 1, output_f ) != 1 || ferror( output_f ))
    {
      free( compressed );
      return Z_ERRNO;
    }
        
  if ( fwrite(compressed, 1, compressed_size, output_f) != compressed_size || ferror(output_f))
    {
      free( compressed );
      return Z_ERRNO;
    }
  return zok;
}

// save compressed data into buffer. Buffer has to be large enough.
int fromCompressedFile( unsigned char * buffer, size_t uncompressed_size, FILE * input_f )
{
  uLongf compressed_size;
  
  if ( fread( &compressed_size, sizeof(uLongf), 1, input_f) != 1 || ferror(input_f))
    {
      return Z_ERRNO;
    }
  
  Bytef * compressed = (Bytef *)calloc( compressed_size, sizeof(Bytef) );

  if ( ( fread(compressed, 1, compressed_size, input_f) != compressed_size) || ferror(input_f) )
    {
      free( compressed );
      return Z_ERRNO;
    }

  int zok = uncompress (buffer, &uncompressed_size, compressed, compressed_size);

  free(compressed);
  return zok;
}
  
  /*
int toCompressed( unsigned char * buffer, size_t uncompressed_size, FILE * output_f )
{

  z_stream strm;

  // allocate deflate state 
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  
  int level = 9;
  size_t compressed_size = 0;

  int ret = deflateInit(&strm, level);
  if (ret != Z_OK)
    {
      return ret; 
    }

  size_t max_compressed = uncompressed_size * 2;
  unsigned char * compressed = (unsigned char *)calloc( max_compressed, sizeof(unsigned char) );

  // compress in one go
  strm.next_in = (unsigned char*)buffer;
  strm.avail_in = uncompressed_size;
  strm.avail_out = max_compressed;
  strm.next_out = compressed;
  
  ret = deflate(&strm, 1);

  assert( ret != Z_STREAM_ERROR);
  compressed_size = max_compressed - strm.avail_out;

  fwrite( &uncompressed_size, sizeof(size_t), 1, output_f );
  fwrite( &compressed_size, sizeof(size_t), 1, output_f );
  
  if (fwrite(compressed, 1, compressed_size, output_f) != compressed_size or ferror(output_f))
    {
      deflateEnd(&strm);
      free( compressed );
      return Z_ERRNO;
    }
  
  deflateEnd(&strm);
  free( compressed );
  return Z_OK;
}
*/


  
  

/*
// allocate memory in buffer
int fromCompressed( unsigned char ** buffer, FILE * input_f )
{
  int ret;
  size_t uncompressed_size;
  size_t compressed_size;
  z_stream strm;

  fread( &uncompressed_size, sizeof(size_t),1, input_f);
  fread( &compressed_size, sizeof(size_t),1, input_f);
  
  unsigned char * compressed = (unsigned char *)calloc( compressed_size, sizeof(unsigned char) );
  (*buffer) = (unsigned char *)calloc( uncompressed_size, sizeof(unsigned char) );
  
  // allocate inflate state 
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  
  ret = inflateInit(&strm);

  if (ret != Z_OK)
    return ret;

  // decompress until deflate stream ends or end of file 

  strm.avail_in = fread(compressed, 1, compressed_size, input_f);
  if (ferror(input_f)) 
    {
      (void)inflateEnd(&strm);
      return Z_ERRNO;
    }

  strm.next_in = compressed;
  strm.avail_out = uncompressed_size;
  strm.next_out = (*buffer);
  ret = inflate(&strm, Z_NO_FLUSH);
  
  assert(ret != Z_STREAM_ERROR);
  // state not clobbered 
  switch (ret) 
    {
    case Z_NEED_DICT:
      ret = Z_DATA_ERROR;
      // and fall through 
    case Z_DATA_ERROR:
    case Z_MEM_ERROR:
	    (void)inflateEnd(&strm);
	    return ret;
    }

  // clean up and return 
  (void)inflateEnd(&strm);

  return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}


*/
