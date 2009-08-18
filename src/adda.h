//--------------------------------------------------------------------------------
// Project adda
//
// Copyright (C) 2003 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: adda.h,v 1.3 2006/06/15 09:19:32 aheger Exp $
//--------------------------------------------------------------------------------    

#ifndef ADDA_H
#define ADDA_H 1

// enable large file system support
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif

#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif

#define _FILE_OFFSET_BITS 64

#include <map>
#include <cstdio>
#include <iostream>
#include <fstream>

//------------------------------------------------------------------------
FILE * openFileForRead( const std::string & filename );
FILE * openFileForWrite( const std::string & filename );

//------------------------------------------------------------------------
typedef long Nid;
typedef fpos_t FileIndex;
typedef int Index;
typedef int Residue;
typedef double Score;

// types for storing alignments in Adda graph
// length of an alignment string
typedef short unsigned int Length;
typedef short unsigned int uResidue;

typedef std::map< Nid, Index > MapNid2Index;

#define MAX_LINE_LENGTH 10000000
#define SEPARATOR '\t'
#define SEPARATOR_RANGE '_'
// end of file token
#define TOKEN "#//\n"

// minimum probability to avoid taking the log of 0
#define SMALL_PROBABILITY 1e-20

//------------------------------------------------------------------------
template< class Array >
void fillParameterArrays( std::ifstream & infile,
			  Array & array)
{
  
  unsigned int x;
  double y;

  while (!infile.eof()) 
  {
    // skip comments
    char c = infile.peek();
    if (c == '#')
      {
	infile.ignore(10000, '\n');
	continue;
      }
    
    infile >> x >> y;
    infile.ignore(10000, '\n');
    if (infile.eof()) break;
    if (array.size() <= x) array.resize( x + 1, 0);
    array[x] = y;
  }
}

//------------------------------------------------------------------------
/** fill values from 0 to first value with first value.
    fill values in between with linear interpolation of adjacent values.
 */
template< class Array >
void interpolateParameterArrays( Array & array) 
{
  
  // set left side
  unsigned int x=0;
  while (x < array.size() && array[x] == 0) { ++x; };

  for (unsigned int y = 0; y < x; ++y) { array[y] = array[x]; };

  double last_value = array[x];

  for ( x+=1 ; x < array.size() - 1; ++x) {
    if (array[x] == 0) {
      unsigned int y = x + 1;
      while (y < array.size() && array[y] == 0) ++y;
      array[x] = last_value + (array[y] - last_value) / (double)(y - x);
    }
    last_value = array[x];
  }
}

#endif

