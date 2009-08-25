//--------------------------------------------------------------------------------
// Project adda
//
// Copyright (C) 2003 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id$
//--------------------------------------------------------------------------------    

/* Convert a sequence graph to a domain graph.
 */

#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>

#include <cstdio>
#include <cassert>
#include <cstring>

#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>

#include "adda.h"

//---------------------------------------------------------------
extern std::string param_file_name_domains;
extern std::string param_file_name_graph;
extern int param_min_residues_overlap;
extern double param_min_coverage;
extern double param_min_relative_overlap;
extern double param_evalue_threshold_trusted_links;

extern int param_loglevel;
extern int param_mode;

struct Domain 
{
	Domain( Residue start = 0, Residue end = 0 ) :
		mStart( start ), mEnd( end )
		{};
		Residue mStart;
		Residue mEnd;
};

typedef std::vector< Domain > Domains;

// Parser class
class LinkProcessor
{
protected:
	struct Link 
	{
		Link( std::string query_token,
				std::string sbjct_token,
				double score,
				double overlap,
				double coverage_query,
				double coverage_sbjct,
				double relative_overlap
		) : mQueryToken( query_token ),
		mSbjctToken( sbjct_token ),
		mScore( score ),
		mOverlap( overlap ),
		mCoverageQuery( coverage_query),
		mCoverageSbjct( coverage_sbjct),
		mRelativeOverlap( relative_overlap )
		{
		}

		std::string mQueryToken;
		std::string mSbjctToken;
		Score mScore;
		double mOverlap;
		double mCoverageQuery;
		double mCoverageSbjct;
		double mRelativeOverlap;
	};


public:
	LinkProcessor( const std::vector<Domains> & domains, MapNid2Index & map_nid2index, int min_residues_overlap ) :
		mDomains(domains), 
		mMapNid2Index( map_nid2index ), 
		mMinResiduesOverlap( min_residues_overlap)
		{
		}

	virtual ~LinkProcessor()
	{
		mLinks.clear();
	};

	void process( std::ofstream & outstream, const char * line )
	{

		mLine = line;
		// skip comments and header
		if (line[0] == '#' || strncmp(line, "query_nid", 9 == 0) ) return;
		
		parse();

		if (mQueryNid == mSbjctNid)
		  return;

		calculateLinks();

		if (param_loglevel >= 3)
		{
			std::cout << "# query_nid=" << mQueryNid 
					<< " sbjct_nid=" << mSbjctNid 
					<< " links computed=" << mLinks.size() 
					<< std::endl;
		}

		selectLinks();

		if (param_loglevel >= 3)
		{
			std::cout << "# query_nid=" << mQueryNid 
					<< " sbjct_nid=" << mSbjctNid 
					<< " links selected=" << mLinks.size() 
					<< std::endl;
		}
		
		printLinks( outstream );

	}


	virtual void printHeader( std::ofstream & outstream ) const
	{
		outstream << "# QDOM:	query domain (nid_from_to)." << std::endl;
		outstream << "# SDOM:	sbjct domain (nid_from_to)." << std::endl;
		outstream << "# SCORE:	overlap score." << std::endl;
		outstream << "# EVAL:	evalue of link." << std::endl;
		outstream << "# QCOV:	coverage of query domain by alignment." << std::endl;
		outstream << "# SCOV:	coverage of sbjct domain by alignment." << std::endl;
		outstream << "# OVL:	common overlap in residues between QDOM and SDOM." << std::endl;
		outstream << "# ROVL:	relative overlap between QDOM and SDOM." << std::endl;
		outstream << "# QDOM\tSDOM\tSCORE\tEVAL\tQCOV\tSCOV\tOVL\tROVL" << std::endl;;
	}

protected:
	
	virtual void parse()
	{
		int query_nid, sbjct_nid;
		float score;

		sscanf( mLine,
			"%d\t%d\t%f\t%i\t%i\t%i\t%i", 
			&query_nid, 
			&sbjct_nid, 
			&score,
			&mQueryFrom, &mQueryTo,
			&mSbjctFrom, &mSbjctTo);

		mScore = score;
		mQueryNid = query_nid;
		mSbjctNid = sbjct_nid;
	}

	virtual void calculateLinks()
	{
		mLinks.clear();

		Residue offset = mSbjctFrom - mQueryFrom;

		Index query_index = mMapNid2Index[mQueryNid];
		Index sbjct_index = mMapNid2Index[mSbjctNid];    

		// iterate over query
		Domains::const_iterator query_it( mDomains[query_index].begin()), 
											query_end( mDomains[query_index].end() );

		for (;query_it != query_end; ++query_it)
		{
			Residue query_domain_from = query_it->mStart;
			Residue query_domain_to   = query_it->mEnd;	

			// check if overlap
			int query_overlap = std::min(mQueryTo, query_domain_to) - std::max(mQueryFrom, query_domain_from);
			int lquery = query_domain_to - query_domain_from;

			if (query_overlap < mMinResiduesOverlap) continue;

			// map to alignment
			int mapped_query_domain_from = std::max(query_domain_from + offset, mSbjctFrom);
			int mapped_query_domain_to   = std::min(query_domain_to + offset, mSbjctTo);

			// check for overlap with domains in sbjct
			// iterate over sbjct
			Domains::const_iterator sbjct_it( mDomains[sbjct_index].begin()), 
								sbjct_end( mDomains[sbjct_index].end() );

			for (;sbjct_it != sbjct_end; ++sbjct_it)
			{
				Residue sbjct_domain_from = sbjct_it->mStart;
				Residue sbjct_domain_to   = sbjct_it->mEnd;	

				int lsbjct = sbjct_domain_to - sbjct_domain_from;

				// check if overlap
				int mapped_overlap = std::min(mapped_query_domain_to, sbjct_domain_to) - std::max(mapped_query_domain_from, sbjct_domain_from);
				int mapped_union = std::max(mapped_query_domain_to, sbjct_domain_to) - std::min(mapped_query_domain_from, sbjct_domain_from);

				if (param_loglevel > 4)
				{
					std::cout << "# link: query=" << mQueryFrom << "-" << mQueryTo << " sbjct=" << mSbjctFrom << "-" << mSbjctTo << std::endl;
					std::cout << "# query_domain="<< query_domain_from << "-" << query_domain_to << std::endl;
					std::cout << "# mapped_domain="<< mapped_query_domain_from << "-" << mapped_query_domain_to << std::endl;
					std::cout << "# sbjct_domain=" << sbjct_domain_from << "-" << sbjct_domain_to << " " << mapped_overlap << std::endl;
				}
				if (mapped_overlap < mMinResiduesOverlap) continue;

				double cov_query = double(mapped_overlap) / lquery;
				double cov_sbjct = double(mapped_overlap) / lsbjct;
				double ovl = double(mapped_overlap) / double(mapped_union);

				char buffer1[MAX_LINE_LENGTH];
				sprintf( buffer1, "%i_%i_%i", (int)mQueryNid, query_domain_from, query_domain_to);

				char buffer2[MAX_LINE_LENGTH];
				sprintf( buffer2, "%i_%i_%i", (int)mSbjctNid, sbjct_domain_from, sbjct_domain_to);

				mLinks.push_back( 
						Link( std::string(buffer1),
								std::string(buffer2),
								mScore,
								mapped_overlap,
								cov_query,
								cov_sbjct,
								ovl ) );
			}
		}
	}

	virtual void selectLinks( )
	{
	}

	virtual void printLinks( std::ofstream & outstream)
	{
		std::vector< Link >::const_iterator it(mLinks.begin()), end(mLinks.end());

		// The score is given by the relative overlap, but the evalue is taken into account: 
		// All links with an evalue of less than
		// param_evalue_threshold_trusted_links get a better score
		// than those with evalues above the threshold.

		for (;it != end;++it)
		{
			outstream
			<< it->mQueryToken << SEPARATOR
			<< it->mSbjctToken << SEPARATOR
			<< ( (it->mScore < param_evalue_threshold_trusted_links) ? 
					(200 - it->mRelativeOverlap * 100) : 
						(100 - it->mRelativeOverlap * 100) ) << SEPARATOR
			<< it->mScore << SEPARATOR
			<< it->mCoverageQuery << SEPARATOR
			<< it->mCoverageSbjct << SEPARATOR
			<< it->mOverlap << SEPARATOR
			<< it->mRelativeOverlap << SEPARATOR
			<< std::endl;
		}
	}

protected:
	const char * mLine;
	Nid mQueryNid;
	Nid mSbjctNid;
	Residue mQueryFrom;
	Residue mQueryTo;
	Residue mSbjctFrom;
	Residue mSbjctTo;
	Score mScore;
	
	std::vector< Link >  mLinks;

	const std::vector< Domains > & mDomains;
	MapNid2Index & mMapNid2Index;

	int mMinResiduesOverlap;

};


class LinkProcessorMax : public LinkProcessor
{
public:
	LinkProcessorMax( const std::vector<Domains> & domains, MapNid2Index & map_nid2index,
			int min_residues_overlap, double min_coverage, double min_relative_overlap ) :
				LinkProcessor( domains, map_nid2index, min_residues_overlap ),
				mMinCoverage( min_coverage ), mMinRelativeOverlap( min_relative_overlap )
				{
				}

protected:

	class SorterByOverlap 
	{
	public:
		SorterByOverlap() {};
		bool operator()( const Link & a, const Link & b) const 
		{
			return a.mRelativeOverlap > b.mRelativeOverlap;
		}
	};

	virtual void selectLinks()
	{
		
		std::sort( mLinks.begin(), mLinks.end(), SorterByOverlap());

		std::vector <Link> old_links;
		std::copy( mLinks.begin(), mLinks.end(), std::back_inserter< std::vector<Link> >(old_links) );
		mLinks.clear();

		std::vector< Link >::const_iterator it(old_links.begin()), end(old_links.end());

		std::map< std::string, bool> written;

		for (;it != end;++it)
		{
	        
			if ( written.find( it->mQueryToken ) != written.end() ||
					written.find( it->mSbjctToken ) != written.end() )
				continue;

			written[it->mQueryToken] = true;
			written[it->mSbjctToken] = true;

			// print, if either domain is mostly covered
			if (it->mCoverageQuery >= mMinCoverage ||
					it->mCoverageSbjct >= mMinCoverage ||
					it->mRelativeOverlap >= mMinRelativeOverlap )
				mLinks.push_back( *it );
		}
		old_links.clear();
	}

	double mMinCoverage;
	double mMinRelativeOverlap;

};

/** read domains and build map_nid2index at the same time.
 */
void readDomains( FILE * file,
		  std::vector< Domains > & domains,
		  MapNid2Index & map_nid2index) 
{
  domains.clear();
  map_nid2index.clear();

  int nid = 0;
  int start = 0;
  int end = 0;

  size_t max_length = MAX_LINE_LENGTH;
  char * buffer = new char[MAX_LINE_LENGTH+1];
  
  while (!feof(file))
    {

      // skip over header
      ssize_t r = getline( &buffer, &max_length, file);
      if (r != -1)
	{

	  // skip comments or header (starting with query_nid)
	  if (buffer[0] == '#' || strncmp(buffer, "nid", 3) == 0 ) continue;

	  if (sscanf(buffer, "%i\t%i\t%i", &nid, &start, &end) != 3)
	    break;
	  
	  Index index;
	  if (map_nid2index.find( nid ) == map_nid2index.end())
	    {
	      index = map_nid2index.size();
	      map_nid2index[nid] = index;
	      domains.resize( index + 1 );
	    }
	  else
	    {
	      index = map_nid2index[nid];      	    
	    }
	  
	  domains[index].push_back( Domain(start, end)) ;
	}

    }
  delete [] buffer;
  
}


//--------------------------------------------------------------------------------
int cadda_convert( const char * filename )
{
  //------------------------------------------------------------
  if (param_loglevel >= 1)
    std::cout << "# reading domains from " << param_file_name_domains << std::endl;
  
  std::vector< Domains > domains;
  MapNid2Index map_nid2index;
  
  {
    FILE * file = openFileForRead( param_file_name_domains );
    readDomains( file, domains, map_nid2index );
    fclose(file );
  }
  
  if (param_loglevel >= 1)
    std::cout << "# read domains for " << domains.size() << " nids." << std::endl;

  LinkProcessor * parser = NULL;
  
  std::ofstream outstream( filename );
	
  switch (param_mode )
    {
    case 0: 
      parser = new LinkProcessor( domains,
				  map_nid2index,
				  param_min_residues_overlap );
		break;
    case 1:
      parser = new LinkProcessorMax( domains,
				     map_nid2index,
				     param_min_residues_overlap,
				     param_min_coverage,
				     param_min_relative_overlap );
      break;
    }

  // parser->printHeader( outstream );
  
  // main loop
  FILE * infile = openFileForRead( param_file_name_graph );
  {

    char * line = new char[MAX_LINE_LENGTH + 1];
    size_t max_length = MAX_LINE_LENGTH;
    
    // skip over header
    ssize_t r = getline( &line, &max_length, infile);
    assert (strncmp( line, "nid", 3) == 0);
    
    while (!feof(infile) && getline( &line, &max_length, infile) != -1)
      {
	if (line[0] == '#') continue;
	parser->process( outstream, line );
      }
    delete [] line;
  }
  
  delete parser;
  
  outstream << TOKEN;
  
  outstream.close();
  return 1;
}
