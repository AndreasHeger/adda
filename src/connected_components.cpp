#include "connected_components.h"
#include <cassert>
#include <vector>
#include <map>
#include <iostream>

template<class T, class C>
Components<T,C>::Components() 
{
}

template<class T, class C>
Components<T,C>::~Components() 
{
}

template<class T, class C>
int Components<T,C>::addToken( const T & a )
{
	Index index = 0;

	// make a persistent copy
	C q(a);
	
	MapToken2VertexIterator it (mMapToken2Vertex.find(q) );
	if ( it == mMapToken2Vertex.end() )
	{
		index = mMapVertex2Token.size();
		mMapToken2Vertex[q] = index;
		mMapVertex2Token.push_back(q);
	}
	else
	{
		index = (*it).second;
	}
	
	return index + 1;
	
}

template<class T, class C>
bool Components<T,C>::add( const T & a, const T & b )
{
  Index index1 = addToken(a);
  Index index2 = addToken(b);

  if (index1 > index2)
    {
      std::swap( index1, index2 );
    }
    
  if ((size_t)index2 + 1 >= mDad.size() )
    {
#ifdef DEBUG
      std::cout << "resizing mDad from " << mDad.size() << " to " << (index2 * 2) + 1 << std::endl;
#endif
      mDad.resize( (index2 * 2) + 1, 0);
    }
  
  if (index1 == index2)
    return false;
  
  bool retval = false;
  
  /* Sedgewick: Algo in C, p 507 */
  {
    int x,y,i,j,t;
    x = i = index1;
    y = j = index2;

#ifdef DEBUG
    std::cout << "before: " << a << "=" << index1-1 <<" " << b << "=" << index2-1
	      << " " << mDad[x] << " " << mDad[y] << std::endl;
#endif

    while (mDad[i] > 0) i = mDad[i];
    while (mDad[j] > 0) j = mDad[j];
    while (mDad[x] > 0)
    {
      t = x;
      x = mDad[x];
      mDad[t] = i;
    }

    while (mDad[y] > 0)
    {
      t = y;
      y = mDad[y];
      mDad[t] = j;
    }

    // nodes belong to two different components.
    if (i != j)
    {
      retval = true;
      if (mDad[j] < mDad[i])
	{
	  mDad[j] += mDad[i] - 1;
	  mDad[i] = j;
	}
      else
	{
	  mDad[i] += mDad[j] - 1;
	  mDad[j] = i;
	}
    }
    
#ifdef DEBUG
    std::cout << "after : " << a << "=" << index1-1 <<" " << b << "=" << index2-1
	      << " " << mDad[index1] << " " << mDad[index2] << " i=" << i << " j=" << j << std::endl;
#endif

  }
  return retval;
}

//------------------------------------------------------------------------
template<class T,class C>
int Components<T,C>::get( const T & v) 
{
  return getComponent( getIndex(v) );
}

//------------------------------------------------------------------------
// Components<T>::Index does not work:
//  expected constructor, destructor, or type conversion before "Components"
template<class T,class C>
int Components<T,C>::getIndex( const T & v) 
{
  MapToken2VertexIterator it;
  
  if ( (it = mMapToken2Vertex.find(v)) == mMapToken2Vertex.end() )
    return 0;
    
  return (*it).second + 1;
}

//------------------------------------------------------------------------
template<class T, class C>
int Components<T,C>::getComponent( Components<T,C>::Index in) 
{
  assert( in < mDad.size() );
  assert( in >= 0 );
  
  if (in == 0) return 0;

  int i,x;
  i = x = in;	
  while (mDad[x] > 0) x = mDad[x];
  if (i != x) mDad[i] = x;
  return x;
}

//------------------------------------------------------------------------
template<class T, class C>
int Components<T,C>::getNumNodes()
{
  return mMapToken2Vertex.size();
}


//------------------------------------------------------------------------
template<class T, class C>
const C & Components<T,C>::getToken(Index i)
{
  assert( i > 0);
  assert( i <= mMapVertex2Token.size() );
  return mMapVertex2Token[i-1];
}


//------------------------------------------------------------------------
template<class T, class C>
void Components<T,C>::reset()
{
  mDad.clear();
  mMapToken2Vertex.clear();
  mMapVertex2Token.clear();
}


// explicit instantiations
template class Components<const char *, std::string>;
template class Components<int,int>;
template class Components<std::string, std::string>;

 
