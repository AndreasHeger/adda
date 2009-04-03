//--------------------------------------------------------------------------------
// Project adda
//
// Copyright (C) 2003 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: adda.cpp,v 1.3 2006/06/15 09:19:32 aheger Exp $
//--------------------------------------------------------------------------------

#include "adda.h"

#include <math.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <cassert>
#include <cstdio>

#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>

// #define DEBUG

using namespace std;

/** working with partitions

    during the flow of the algorithm, partitions get deleted and added.
    Usually partitions are split, which means the full partition is deleted
    and two new partitions are created.

    score changes, if
    a partition is deleted:
	remove all pairwise scores between deleted partition
	and all other current partitions
    a partition is added:
	calculate pairwise scores between all new partitions.

    In order to do this in a fully incrementally manner, a matrix of all pairwise
    scores would have to be kept. The matrix would change in size continuously. I
    think it impractical.

    Instead the difference in score is calculated by summing over the differences
    between all pairwise scores for the deleted old and the added two new partitions.

    ranges are like python ranges, i.e. from:last_element+1

    - access links from file, do not keep in memory.
    - read tree from file

    have some troubles with streampos. Have to use unsigned long long, as
    there is no input (<<) operator for streampos in sstream anymore (there
    was one in strstream).

    -> alternatively: prepares links file and build index on the fly.
 */

// global variables
static vector<double> param_overhangs;
static vector<double> param_transfers;

extern unsigned int param_loglevel;

/* file names */
extern  std::string param_file_name_trees;
extern  std::string param_file_name_graph;
extern  std::string param_file_name_index;
extern  std::string param_file_name_nids;
extern  std::string param_file_name_transfers;
extern  std::string param_file_name_domains;

/* various options */
extern  bool param_disallow_shortening;
extern  bool param_descend;
extern  int param_max_iterations;
extern  bool param_use_file_nids;

/* function parameters for sigmoid */
extern  double param_resolution;

extern  double param_real_k;	// smooth sigmoid (has to incorporate the resolution (real k = resolution * k)
extern  double param_real_c;	// average domain size of 100
extern  double param_real_max;
extern  double param_real_min;

extern  bool param_relative_overhang;
extern  double param_real_e; // Probability of alignment ending directly in alignment is 5%
extern  double param_real_f;

static double param_e = 0;
static double param_c = 0;
static double param_k = 0;

// safety threshold for small overlaps (should be in the range of the resolution);
extern  int param_threshold_overlap;

extern  bool param_only_query;	// if true, calculate only score for query


//--------------------------------------------------------------------------------
// use int for Residue, as sometimes the difference is
// taken between two Residues and unsigned int gives underflow
typedef unsigned int Node;

typedef std::pair< Residue, Residue> Range;
typedef std::vector< Range > Ranges;

struct Partition {
	Partition( Node xnode,
			Residue xfrom,
			Residue xto) :
				node(xnode), from(xfrom), to(xto) {};
				Node node;
				Residue from;
				Residue to;
};

typedef std::list< Partition > PartitionList;
typedef std::vector< PartitionList > Partitions;

struct Link {
	Link( Residue xquery_from,
			Residue xquery_to,
			Index xsbjct_index,
			Residue xsbjct_from,
			Residue xsbjct_to) :
				query_from(xquery_from), query_to(xquery_to),
				sbjct_index(xsbjct_index),
				sbjct_from(xsbjct_from), sbjct_to(xsbjct_to) {};
				Residue query_from;
				Residue query_to;
				Index sbjct_index;
				Residue sbjct_from;
				Residue sbjct_to;
};

typedef std::list< Link > LinkList;
typedef std::vector< LinkList > Links;

typedef std::map<Nid,FileIndex> FileIndexMap;
typedef std::map<Nid,Index> IndexMap;
typedef std::vector<Nid> NidMap;

struct TreeNode
{
	TreeNode( int p, int l, int r, Residue f, Residue t) :
		mParent(p), mLeftChild(l), mRightChild(r), mFrom(f), mTo(t) {}
	int mParent;
	int mLeftChild;
	int mRightChild;
	Residue mFrom;
	Residue mTo;
};


typedef std::vector< TreeNode > Tree;
typedef std::vector<Tree> Trees;


std::ostream & operator<<( std::ostream & output, const Partition & src)
{
	output << " node=" << src.node << " from=" << src.from << " to=" << src.to;
	return output;
}

std::istream & operator>>( std::istream & input, Partition & target) {
	return input;
}

std::ostream & operator<<( std::ostream & output, const PartitionList & src) {
	std::copy( src.begin(), src.end(), std::ostream_iterator< Partition >( std::cout, "\n"));
	return output;
}

std::ostream & operator<<( std::ostream & output, const Link & src) {
	output << " sbjct_index=" << src.sbjct_index
	<< " query_from=" << src.query_from << " query_to=" << src.query_to
	<< " sbjct_from=" << src.sbjct_from << " sbjct_to=" << src.sbjct_to;
	return output;
}


std::ostream & operator<<( std::ostream & output, const LinkList & src) {
	std::copy( src.begin(), src.end(), std::ostream_iterator< Link >( std::cout, "\n"));
	return output;
}

std::ostream & operator<<( std::ostream & output, const Range & src) {
	output << "from=" << src.first << " to=" << src.second;
	return output;
}

std::ostream & operator<<( std::ostream & output, const TreeNode & src) {
	output << src.mParent << "\t" << src.mLeftChild << "\t" << src.mRightChild << "\t" << src.mFrom << "\t" << src.mTo;
	return output;
}

void printTrees( std::ostream & output,
		const Trees & trees, const NidMap & map_index2nid )
{
	Index index = 0;
	for (;index<(Index)trees.size();++index)
	{
		Nid nid = map_index2nid[index];
		Index ii = 0;
		for (; ii < (Index)trees[index].size(); ++ii)
			output << index << "\t" << nid << "\t" << ii << "\t" << trees[index][ii] << std::endl;
	}
}

void printPartitions( std::ostream & output ,
		const Partitions & partitions, const NidMap & map_index2nid )
{
	Index index = 0;
	for (;index<(Index)partitions.size();++index)
	{
		Nid nid = map_index2nid[index];
		PartitionList::const_iterator it(partitions[index].begin()), end(partitions[index].end());
		for (;it!=end;++it)
			output << nid << "\t" << it->from << "\t" << it->to << "\t" << endl;
	}
}

inline void printSection()
{
	std::cout << "##-------------------------------------------------------" << std::endl;
}

inline int convert(int x) { return (int)(floor( ( (x-1) / param_resolution) ) ); }

/** retrieve all links from nid using the table table_links starting at position
    index
 */
template< class OutputIter >
void fillLinks( FILE * infile,
		const FileIndex & index,
		const Nid & nid,
		const IndexMap & map_nid2index,
		OutputIter it)
{

	fsetpos( infile, &index );

#ifdef DEBUG
	// check if you are correctly positioned
	{
		Nid query_nid, sbjct_nid;
		float score;
		Residue query_from, query_to, sbjct_from, sbjct_to;
		fscanf( infile,
			"%ld\t%ld\t%f\t%i\t%i\t%i\t%i",
			&query_nid,
			&sbjct_nid,
			&score,
			&query_from, &query_to,
			&sbjct_from, &sbjct_to);

		if (query_nid != nid)
		{
		        std::cerr << "# positioning error in fillLinks: nid=" \
				  << nid << "\tquery_nid=" << query_nid << std::endl;
			exit(EXIT_FAILURE);
		}
		fsetpos( infile, &index );
	}
#endif

	char buffer[MAX_LINE_LENGTH+1];

	while (!feof(infile))
	{
		Nid query_nid, sbjct_nid;
		float score;
		Residue query_from, query_to, sbjct_from, sbjct_to;

		{
			int r = fscanf( infile,
				"%ld\t%ld\t%f\t%i\t%i\t%i\t%i",
				&query_nid,
				&sbjct_nid,
				&score,
				&query_from, &query_to,
				&sbjct_from, &sbjct_to);
			assert( r == 7);
		}

		char * r = fgets( buffer, MAX_LINE_LENGTH, infile );
		assert( r != NULL );

#ifdef DEBUG
		std::cout << "read the following: " << std::endl;
		std::cout << " " << query_nid << " " << sbjct_nid << " " << query_from << "-" << query_to \
			  << " " << sbjct_from << "-" << sbjct_to << " " << endl;
#endif

		if (query_nid != nid)
			break;

		if (query_nid == sbjct_nid)
			continue;

		IndexMap::const_iterator i;
		if ( (i = map_nid2index.find( sbjct_nid )) != map_nid2index.end())
		{
			Index sbjct_index = i->second;
			*it = Link(convert(query_from),
				   convert(query_to)+1,
				   sbjct_index,
				   convert(sbjct_from),
				   convert(sbjct_to)+1);
			++it;
		}
	}

	// reset stream, move away from eof.
	if (feof(infile))
		rewind(infile);

}

template< class T>
inline void printMap( const T & index, const char * title = "" )
{
	typedef typename T::const_iterator iterator;
	std::cout << "# ";
	iterator it(index.begin()), end(index.end());

	std::cout << title ;

	for (; it != end; ++ it)
	  std::cout << it->first << "->" << it->second << ":";
	std::cout << std::endl;
}

//--------------------------------------------------------------------------------
/** read trees from file
    nodes are sorted by nid and node!
    input format of trees file is
    nid	node	parent	level	xfrom	xto
 */
void fillTrees( ifstream & infile,
		const IndexMap & map_nid2index,
		Trees & trees)
{
    
  while (!infile.eof())
    {
      char c = infile.peek();

      // skip comments or header line (starting with 'nid')
      if (c == '#' || c == 'n')
	{
	  infile.ignore(10000, '\n');
	  continue;
	}
		
      Nid nid;
      int node, parent, level, xfrom, xto;
		
      infile >> nid >> node >> parent >> level >> xfrom >> xto;
      infile.ignore(10000, '\n');
      if (infile.eof()) break;
		
      IndexMap::const_iterator it = map_nid2index.find(nid);
      
      if (it != map_nid2index.end())
	{
	  Index index = it->second;
	  // std::cout << nid << "\t" << parent << "\t" << xfrom << "\t" << xto << "\t" << index <<endl;
	  trees[index].push_back( TreeNode(parent,0,0,xfrom,xto) );
	}
    }

	Trees::iterator it(trees.begin()), end(trees.end());
	for (;it!=end;++it)
	{
		Tree & t = *it;
		// skip root, child can not be 0 (unless in leaves)
		for (unsigned int i = 1; i < t.size(); ++i)
		{
			int parent = t[i].mParent;
			if (t[parent].mLeftChild)
			{
				t[parent].mRightChild = i;
			}
			else
			{
				t[parent].mLeftChild = i;
			}
		}
	}
}

//--------------------------------------------------------------------------------
/** get child partitions for a given partition
 */
template< class OutputIter >
void fillPartitionsWithChildren( const Trees & trees,
		const Index index,
		const Node parent_node,
		OutputIter it)
{

	int left = trees[index][parent_node].mLeftChild;
	int right = trees[index][parent_node].mRightChild;
	if (left)
	{
		*it = Partition( left, trees[index][left].mFrom, trees[index][left].mTo );
		++it;
	}

	if (right)
	{
		*it = Partition( right, trees[index][right].mFrom, trees[index][right].mTo);
		++it;
	}
}

//------------------------------------------------------------------------
// fill vector of nids with list of components

template< class OutputIter >
void fillNidsFromFile( ifstream & infile,
		OutputIter  it)
{  
#define MAX_SEQUENCE_LENGTH 100000
	// check that the header is "nid"
	std::string header("");
	infile >> header;
	assert( header == "nid" );
	infile.ignore(MAX_SEQUENCE_LENGTH, '\n');
	
	Nid nid;
	
	while (!infile.eof())
	{
	  // skip comments 
	  char c = infile.peek();
	  if (c == '#' ) 
	    {
	      infile.ignore(MAX_SEQUENCE_LENGTH, '\n');
	      continue;
	    }
	  
	  infile >> nid;
	  infile.ignore(MAX_SEQUENCE_LENGTH, '\n');
	  if (infile.eof()) break;
	  *it = nid;
	  ++it;
	}
}

inline double getTransfer( Residue & transfer )
{
	if (transfer < 0)
		return param_transfers[0];
	else if (transfer >= (Residue)param_transfers.size())
		return param_transfers[param_transfers.size()-1];
	else
		return param_transfers[transfer];
}

inline double getOverhang( Residue & overhang )
{
	return param_real_f * exp( -param_e * overhang);
}


//--------------------------------------------------------------------------------
/*
  Calculate transfer score between two domains on two sequences and an alignment
  between them.

  Mapping is necessary (see case three).


  1.
	   [----transfer-----]
   [-------------l1-----------------]
           -------------------
   [------------l2------------------]
  2.
                         [-transfer-]
   [-------------l1-----------------]
                         -------------------
   [------------l2------------------]
  3.
                         [-transfer-]
   [-------------l1-----------------]
              ----------------------------------------------------------
			 [------------l2------------------]


  The transfer score is:

  OF = log( (double)(l1-transfer) * (double)(l2-transfer))


 */
inline double calculateScore(
		const Range & part1, const Range & part2,
		const Range & link1, const Range & link2 )
{
	Residue l1 = part1.second - part1.first;
	Residue l2 = part2.second - part2.first;

	// take min, so that transfer will not be negative
	// make sure all values are treated as signed values, as the differences can
	// be negative
	Residue lali = std::min( link1.second - link1.first, link2.second - link2.first);
	Residue start = std::max( part1.first  - link1.first + link2.first, link2.first);
	Residue end   = std::min( part1.second - link1.first + link2.first, link2.second);

	Residue transfer = std::min( part2.second, end) - max( part2.first, start);

	// has to be one, as due to the scaling I might have small overlap
	// which would get penalized heavily!!!
	if (transfer <= param_threshold_overlap)
		return 0;

	// surprise score for splitting the alignment into a domain of size transfer
	double p = 1.0 - getTransfer( transfer );
	if (p <= 0) p = SMALL_PROBABILITY;
	
	double s = -log(p);

	// surprise score for missing out on residues:
	Residue o1, o2;
	if (param_relative_overhang)
	{
		o1 = l1 - transfer * 100 / l1;
		o2 = l2 - transfer * 100 / l2;
	} else
	{
		o1 = l1 - transfer;
		o2 = l2 - transfer;
	}

	// surprise score for missing out on residues:
	double p1 = getOverhang( o1 );
	double s1 = -log(p1);

	double p2 = getOverhang( o2 );
	double s2 = -log(p2);

	if (param_loglevel >= 3)
	{
		cout << endl << "part1=" << part1 << endl << "part2=" << part2 << endl << "link1=" << link1 << endl << "link2=" << link2 << endl;
		cout << " ( start=" << start << ", end=" << end
		<< ",l1=" << l1 << ",l2=" << l2
		<< ",lali=" << lali << ",t=" << transfer
		<< ")" << endl;
		cout << " (p=" << p  << ",s="  << s
		<< ",p1=" << p1 << ",s1=" << s1
		<< ",p2=" << p2 << ",s2=" << s2
		<< ")" << endl;
	}

	if (param_only_query)
		return s + s2;
	else
		return s + s1 + s2;
}


double calculatePartitionScore( LinkList & links,
		Partitions & partitions,
		Partition  & old_partition,
		PartitionList & new_partitions)
{

	double delta_score = 0;

	//------------------------
	// 1. loop over links
	LinkList::iterator it(links.begin()), it_end( links.end());

	for (; it != it_end; ++it)
	{

		Index sbjct_index = it->sbjct_index;

		//------------------------
		// 2. loop over partitions in sbjct sequence linked to query sequence
		PartitionList::iterator pit(partitions[sbjct_index].begin()), pend(partitions[sbjct_index].end());

		for (;pit!=pend;++pit)
		{

			if (param_loglevel >= 4)
				cout << "------> testing partition idx=" << sbjct_index << " from=" << pit->from << " to=" << pit->to << endl;

			// calculate score for original partition
			double old_score = calculateScore( Range( pit->from, pit->to),
					Range( old_partition.from, old_partition.to),
					Range( it->sbjct_from, it->sbjct_to),
					Range( it->query_from, it->query_to));


			// calculate score for new partitions
			double new_score = 0;
			PartitionList::iterator nit( new_partitions.begin()), nit_end( new_partitions.end());
			for (; nit != nit_end; ++nit)
			{
				new_score += calculateScore( Range( pit->from, pit->to),
						Range( nit->from, nit->to),
						Range( it->sbjct_from, it->sbjct_to),
						Range( it->query_from, it->query_to));

			}

			delta_score += new_score - old_score;

			if (param_loglevel >= 4)
				cout << endl;
			if (param_loglevel >= 3)
				cout << "------> " << sbjct_index
				<< " new_score=" << new_score
				<< " old_score=" << old_score
				<< " inc=" << new_score - old_score
				<< " delta=" << delta_score
				<< endl;
		}
	}
	return delta_score;
}

//--------------------------------------------------------------------------------
void fillFileIndexMap( FileIndexMap & map_nid2fileindex, std::string & file_name_index)
{

	FILE * file = fopen(file_name_index.c_str(), "r");

	if (file == NULL)
	{
		std::cerr << "could not open filename with indices: " << file_name_index << std::endl;
		exit(EXIT_FAILURE);
	}

	while(!feof(file))
	{

		Nid nid = 0;
		FileIndex index;

		{
		  int r = fread(&nid,sizeof(Nid), 1, file);
		  assert( r == 1 );
		}
		if (feof(file)) break;
		{
		  int r = fread(&index, sizeof(FileIndex), 1, file);
		  assert( r == 1 );
		}
		map_nid2fileindex[nid] = index;
	}

	fclose(file);

}


// global variables
IndexMap global_map_nid2index;
NidMap global_map_index2nid;
FILE * global_file_links = NULL;
FileIndexMap global_map_nid2fileindex;
Partitions global_partitions;
Trees global_trees;

//--------------------------------------------------------------------------------
/** optimize partitions.

    below is a greedy optimizer. Every current node in the tree is split. The split
    is retained only, if it is an improvement over a previous assignment.

    The order of commands is as follows:

    	0. let improvement = 0;
		1. get next partition p_i
		2. split partitition into p_i1 and p_i2
		3. calculate improvement
		4. if improvement > 0:
				delete old partition p_i
				add new partitions p_i1, p_i2
				set improvement > 0

		save partitions to file

	return improvement
 */
double cadda_optimise_iteration()
{

	if (param_loglevel >= 3)
	{
		std::cout << "# ------------------------------------------------------------------------" << endl;
		std::cout << "# partitions at start of iteration" << endl;
		printPartitions(
				std::cout,
				global_partitions,
				global_map_index2nid );
	}

	double improvement = 0;

	Index index = 0;
	for (; index < (Index)global_map_index2nid.size(); ++index)
	{
		Nid nid = global_map_index2nid[index];
		if (param_loglevel >= 2)
			cout << "--> checking split of sequence " << nid << "(index=" << index << ")" << endl;

		PartitionList::iterator it(global_partitions[index].begin()), end(global_partitions[index].end());

		while (it!=end)
		{

			Node node = it->node;

			PartitionList new_partitions;

			fillPartitionsWithChildren(
					global_trees,
					index,
					node,
					back_insert_iterator< PartitionList >(new_partitions));

			if (new_partitions.size() == 0)
			{
				++it;
				continue;
			}

			if (param_loglevel >= 2)
			{
			  cout << "# ----> new partitions for sequence " << nid << ": ";
			  std::copy(
				    new_partitions.begin(),
				    new_partitions.end(),
				    ostream_iterator< Partition >( std::cout, ";"));
			  cout << endl;
			}

			if (param_disallow_shortening && (new_partitions.size() != 2))
			{
				++it;
				continue;
			}

			LinkList links;

			fillLinks( global_file_links,
				   global_map_nid2fileindex.find(nid)->second,
				   nid,
				   global_map_nid2index,
				   back_insert_iterator< LinkList >(links));

			if (param_loglevel >= 3)
				cout << "# --> found " << links.size() << " links" << endl;

			double score = calculatePartitionScore( links, global_partitions, *it, new_partitions );

			if (score < 0)
			{
				if (param_loglevel >= 2)
					cout << "# ----> substituting partitition" << endl;

				global_partitions[index].insert( it, new_partitions.begin(), new_partitions.end());
				it = global_partitions[index].erase( it );
				improvement += score;

				if (param_descend)
					for (unsigned int i = 0; i < new_partitions.size(); ++i)
						--it;
			}
			else
			{
				if (param_loglevel >= 2)
					std::cout << "# ----> keeping partitition" << endl;
				++it;
			}

			if (param_loglevel >= 4)
				std::cout << "# ----------------------------------------------------------------" << endl;

		}

	}

	if (param_loglevel >= 1)
		std::cout << "# --> improvement=" << -improvement << std::endl;

	return -improvement;
}

//--------------------------------------------------------------------------------
int cadda_optimise_destroy()
{

	fclose( global_file_links );
	global_map_nid2index.clear();
	global_map_index2nid.clear();
	global_map_nid2fileindex.clear();
	global_partitions.clear();
	global_trees.clear();

	return 1;
}


long cadda_optimise_get_num_partitions()
{
	long total = 0;
	for (Index i = 0; i<(Index)global_partitions.size();++i)
	{
		total += global_partitions[i].size();
	}
	return total;
}

int cadda_optimise_load_partitions( const char * filename )
{
	std::cout << "# partitions loaded from " << filename << endl;
	std::cout << "# partitions loading not implemented!" << endl;
	return 0;
}

int cadda_optimise_save_partitions( const char * filename )
{
	std::ofstream outfile;

	outfile.open( filename );

	if (param_loglevel >= 2)
		std::cout << "# partitions written to " << filename << endl;

	outfile << "nid\tstart\tend" << std::endl;	
	printPartitions( outfile, global_partitions, global_map_index2nid);
	outfile << TOKEN;
	outfile.close();

	return 1;
}


//--------------------------------------------------------------------------------
int cadda_optimise_initialise()
{

	param_e = param_real_e * param_resolution;
	param_c = param_real_c / param_resolution;
	param_k = param_real_k / param_resolution;

	if (param_real_f == 0) param_real_f = param_real_e;

	// safety threshold for small overlaps (should be in the range of the resolution);
	param_threshold_overlap = (int)(12.0 / param_resolution);

	/*------------------------------------------------------------------*/
	// read distribution function from file
	//
	{
		if (param_loglevel >= 1)
		{
			cout << "# retrieving values for transfers from " << param_file_name_transfers << std::endl;
			std::cout.flush();
		}

		ifstream fin(param_file_name_transfers.c_str());
		if (!fin)
		{
			std::cerr << "could not open filename with transfers: " << param_file_name_transfers << std::endl;
			exit(EXIT_FAILURE);
		}
		fillParameterArrays( fin, param_transfers );
		fin.close();

		interpolateParameterArrays(param_transfers);

		if (param_loglevel >= 5)
			for (unsigned int x = 0; x < param_transfers.size(); ++x)
				std::cout << "# " << x << "\t" << param_transfers[x] << std::endl;

	}

	/*------------------------------------------------------------------*/
	// retrieve nids of component
	if (param_loglevel >= 1)
	{
		std::cout << "# retrieving nids from " << param_file_name_nids << std::endl;
		std::cout.flush();
	}

	{
		ifstream fin(param_file_name_nids.c_str());
		if (!fin) {
			std::cerr << "could not open filename with nids: " << param_file_name_nids << std::endl;
			exit(EXIT_FAILURE);
		}
		fillNidsFromFile( fin,
				back_insert_iterator< NidMap >(global_map_index2nid));
		fin.close();
	}

	if (param_loglevel >= 1)
		std::cout << global_map_index2nid.size() << " nids in partition " << endl;
	if (param_loglevel >= 5)
		std::copy( global_map_index2nid.begin(),
					global_map_index2nid.end(),
					ostream_iterator<Nid>( std::cout, "\n"));

	/*------------------------------------------------------------------*/
	// make map of nid->index
	{
		NidMap::iterator it(global_map_index2nid.begin()), end(global_map_index2nid.end());
		for (Index i=0;it!=end;++it,++i)
			global_map_nid2index[*it] = i;
	}

	if (param_loglevel >= 5)
	{
		NidMap::iterator it(global_map_index2nid.begin()),
			end(global_map_index2nid.end());
		std::cout << "# map_index2nid\tmap_index2nid" << std::endl;
		for (Index i=0;it!=end;++it,++i)
			cout << "# " << i << "\t" << *it << "\t" << global_map_nid2index[*it] << "\t" << (global_map_nid2index.find(*it) != global_map_nid2index.end()) << std::endl;
	}

	/*------------------------------------------------------------------*/
	// read tree
	if (param_loglevel >= 1)
	{
	  cout << "# retrieving trees: ";
	  std::cout.flush();
	}

	global_trees.resize(global_map_index2nid.size());
	{
		ifstream fin(param_file_name_trees.c_str());
		if (!fin) {
			std::cerr << "could not open filename with trees: " << param_file_name_trees << std::endl;
			exit(EXIT_FAILURE);
		}

		fillTrees( fin,
					global_map_nid2index,
					global_trees);

		fin.close();
	}

	if (param_loglevel >= 1)
	{
		Trees::iterator it(global_trees.begin()), end(global_trees.end());
		int count = 0;
		for (;it!=end;++it) count += (it->size() > 0);
		std::cout << count << " trees found " << endl;
	}

	if (param_loglevel >= 5)
	{
		printSection();
		printTrees( std::cout,
					global_trees,
					global_map_index2nid );
	}


	/*------------------------------------------------------------------*/
	// open links file: read indices
	{
		if (param_loglevel >= 1)
		{
			std::cout << "## retrieving indices from " << param_file_name_index << "." << std::endl;
		}

		fillFileIndexMap( global_map_nid2fileindex, param_file_name_index );

		if (param_loglevel >= 1)
		{
			std::cout << "## retrieved " << global_map_nid2fileindex.size() << " indices." << std::endl;
		}
	}

	/*------------------------------------------------------------------*/
	if (param_loglevel >= 1)
		cout << "# opening links file: " << std::endl;

	global_file_links = fopen(param_file_name_graph.c_str(),"r");

	if (global_file_links == NULL)
	{
		std::cerr << "could not open filename with links: " << param_file_name_graph << std::endl;
		exit(EXIT_FAILURE);
	}

	/*------------------------------------------------------------------*/
	// fill partitions with initial values
	global_partitions.resize(global_map_index2nid.size());
	{
		for (Index index = 0; index < (Index)global_map_index2nid.size(); ++index)
			if (global_trees[index].size())
				global_partitions[index].push_back(
						Partition( 0,
						global_trees[index][0].mFrom,
						global_trees[index][0].mTo) );
	}

	if (param_loglevel >= 5)
	{
		std::cout << "# partitions at beginning " << endl;
		std::copy( global_partitions.begin(),
					global_partitions.end(),
					ostream_iterator< PartitionList >( std::cout, ""));
	}

	return 1;
}





























