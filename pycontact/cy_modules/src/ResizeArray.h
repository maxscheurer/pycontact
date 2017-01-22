/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ResizeArray.h,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.45 $	$Date: 2010/12/16 04:08:38 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Automatically-adjusting single-dim array template.
 * 
 * LICENSE:
 *   UIUC Open Source License
 *   http://www.ks.uiuc.edu/Research/vmd/plugins/pluginlicense.html
 *
 ***************************************************************************/
#ifndef RESIZEARRAY_TEMPLATE_H
#define RESIZEARRAY_TEMPLATE_H

#include <string.h>

/// A template class which implements a dynamically-growing, automatically
/// resizing array of data of a given type.  Elements in the array may be
/// accessed via the [] operator.  When new data is added to the end of an
/// array, the size of the array is automatically increased if necessary.
///
/// XXX Do not parametrize this class with a datatype which cannot be
///     shallow-copied!  This class uses memcpy to resize, and therefore
///     classes which contain dynamically-allocated memory blocks will
///     crash and burn if the ResizeArray ever gets resized.
template<class T>
class ResizeArray {
private:
  T *allocate(size_t n) { return new T[n]; }
  void deallocate(T *p) { delete [] p; }

  T *data;      ///< list of items, and pointer to current item.
  int sz;       ///< max number of items that can be stored in the array
  int currSize; ///< largest index used + 1

public:
  /// Constructor
  /// The first argument is the initial internal size of the array, i.e. the
  /// initial number of elements for which to allocate memory (although the
  /// initial external size of the array will be zero).  
  ResizeArray(int s = 3) {
    currSize = 0;
    sz = (s > 0 ? s : 10);
    data = allocate(sz); 
  }

  ~ResizeArray() {
    deallocate(data);
  }
  
  int num(void) const { return currSize; } ///< current size of array 
  T& operator[](int N) { return data[N]; } ///< unchecked accessor, for speed
  T const& operator[](int N) const { return data[N]; } ///< a const version of above

  /// add a new element to the end of the array.  Return index of new item.
  void append(const T& val) {
    if (currSize == sz) {    // extend size of array if necessary
      int newsize = (int)((float)sz * 1.3f);

      // guarantee minimum required size increase, since the scaled value
      // may truncate back to the original size value when the initial number
      // of elements is very small.
      if (newsize == sz)
        newsize++;

      // shallow copy the data to a newly allocated block since we can't
      // do something better like realloc()
      T *newdata = allocate(newsize); 
      memcpy(newdata, data, currSize * sizeof(T));
      deallocate(data); 
    
      // save new values
      data = newdata;
      sz = newsize;
    }
    data[currSize++] = val;
  }

  /// remove an item from the array, shifting remaining items down by 1
  void remove(int n) {
    if (n < 0 || n >= currSize) return;
    for (int i=n; i<currSize-1; i++)
      data[i] = data[i+1];
    currSize--;
  }

  /// remove the last item from the array, unchecked for speed
  T& pop() {
    currSize--;
    return data[currSize]; 
  }

  /// delete entire array by defining size to be empty
  void clear() {
    currSize = 0;
  }

  /// truncate the array by defining the size to be N items less
  void truncatelastn(int N) {
    currSize -= N;
    if (currSize < 0) 
      currSize=0;
  }

  /// scan the array until the first item that matches in the array is
  /// found.  Return the index if found, (-1) otherwise.
  int find(const T& val) {
    int i;
  
    for(i=0; i < currSize; i++) {
      if(data[i] == val) 
        return i;
    }
  
    return -1;
  }
};

#endif

