// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiIndexedVector.hpp"
#define TINY_ELEMENT 1.0e-30
//#############################################################################

void
OsiIndexedVector::clear()
{
  int i;
  for (i=0;i<nElements_;i++) {
    elements_[indices_[i]]=0.0;
  }
  nElements_ = 0;
  clearBase();
  delete [] packedElements_;
  packedElements_=NULL;
}

//#############################################################################

void
OsiIndexedVector::empty()
{
  delete [] indices_;
  indices_=NULL;
  delete [] elements_;
  elements_=NULL;
  nElements_ = 0;
  capacity_=0;
  clearBase();
  delete [] packedElements_;
  packedElements_=NULL;
}

//#############################################################################

OsiIndexedVector &
OsiIndexedVector::operator=(const OsiIndexedVector & rhs)
{
  if (this != &rhs) {
    clear();
    gutsOfSetVector(rhs.capacity_,rhs.nElements_, rhs.indices_, rhs.elements_,
      OsiPackedVectorBase::testForDuplicateIndex(),
      "operator=");
  }
  return *this;
}

//#############################################################################

OsiIndexedVector &
OsiIndexedVector::operator=(const OsiPackedVectorBase & rhs)
{
  if (this != &rhs) {
    clear();
    gutsOfSetVector(rhs.getNumElements(), rhs.getIndices(), rhs.getElements(),
      OsiPackedVectorBase::testForDuplicateIndex(),
      "operator= from base");
  }
  return *this;
}

//#############################################################################

void
OsiIndexedVector::borrowVector(int size, int numberIndices, int* inds, double* elems)
{
  empty();
  capacity_=size;
  nElements_ = numberIndices;
  indices_ = inds;  
  elements_ = elems;
  
  // whole point about borrowvector is that it is lightweight so no testing is done
}

//#############################################################################

void
OsiIndexedVector::returnVector()
{
  indices_=NULL;
  elements_=NULL;
  nElements_ = 0;
  capacity_=0;
  clearBase();
  delete [] packedElements_;
  packedElements_=NULL;
  // and clear index set
  clearIndexSet();
}

//#############################################################################

void
OsiIndexedVector::setVector(int size, const int * inds, const double * elems,
                            bool testForDuplicateIndex)
{
  clear();
  gutsOfSetVector(size, inds, elems, testForDuplicateIndex, "setVector");
}
//#############################################################################


void 
OsiIndexedVector::setVector(int size, int numberIndices, const int * inds, const double * elems)
{
  clear();
  gutsOfSetVector(size, numberIndices, inds, elems, false, "setVector");
}
//#############################################################################

void
OsiIndexedVector::setConstant(int size, const int * inds, double value,
                              bool testForDuplicateIndex)
{
  clear();
  gutsOfSetConstant(size, inds, value, testForDuplicateIndex, "setConstant");
}

//#############################################################################

void
OsiIndexedVector::setFull(int size, const double * elems,
                          bool testForDuplicateIndex) 
{
  // Clear out any values presently stored
  clear();
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
  // and clear index set
  clearIndexSet();
  
  if (size<0)
    throw CoinError("negative number of indices", "setFull", "OsiIndexedVector");
  
  reserve(size);
  nElements_ = 0;
  // elements_ array is all zero
  int i;
  for (i=0;i<size;i++) {
    int indexValue=i;
    double value=elems[i];
    if (fabs(value)>=TINY_ELEMENT) {
      elements_[indexValue]=value;
      indices_[nElements_++]=indexValue;
    }
  }
}
//#############################################################################

/** Access the i'th element of the full storage vector.  */
double &
OsiIndexedVector::operator[](int index) const
{
  if ( index >= capacity_ ) 
    throw CoinError("index >= capacity()", "[]", "OsiIndexedVector");
  if ( index < 0 ) 
    throw CoinError("index < 0" , "[]", "OsiIndexedVector");
  double * where = elements_ + index;
  return *where;
  
}
//#############################################################################

void
OsiIndexedVector::setElement(int index, double element)
{
  if ( index >= nElements_ ) 
    throw CoinError("index >= size()", "setElement", "OsiIndexedVector");
  if ( index < 0 ) 
    throw CoinError("index < 0" , "setElement", "OsiIndexedVector");
  elements_[indices_[index]] = element;
}

//#############################################################################

void
OsiIndexedVector::insert( int index, double element )
{
  if ( index < 0 ) 
    throw CoinError("index < 0" , "setElement", "OsiIndexedVector");
  if (index >= capacity_)
    reserve(index+1);
  if (elements_[index])
    throw CoinError("Index already exists", "insert", "OsiIndexedVector");
  indices_[nElements_++] = index;
  elements_[index] = element;
  // and clear index set
  clearIndexSet();
  delete [] packedElements_;
  packedElements_=NULL;
}

//#############################################################################

void
OsiIndexedVector::add( int index, double element )
{
  if ( index < 0 ) 
    throw CoinError("index < 0" , "setElement", "OsiIndexedVector");
  if (index >= capacity_)
    reserve(index+1);
  if (elements_[index]) {
    element += elements_[index];
    if (fabs(element)>= TINY_ELEMENT) {
      elements_[index] = element;
    } else {
      elements_[index] = TINY_ELEMENT*1.0e-10;
    }
  } else if (fabs(element)>= TINY_ELEMENT) {
    indices_[nElements_++] = index;
    elements_[index] = element;
   }
  // and clear index set
  clearIndexSet();
  delete [] packedElements_;
  packedElements_=NULL;
}

//#############################################################################

int
OsiIndexedVector::clean( double tolerance )
{
  int number = nElements_;
  int i;
  nElements_=0;
  for (i=0;i<number;i++) {
    int indexValue = indices_[i];
    if (fabs(elements_[indexValue])>=tolerance) {
      indices_[nElements_++]=indexValue;
    } else {
      elements_[indexValue]=0.0;
    }
  }
  // and clear index set
  clearIndexSet();
  delete [] packedElements_;
  packedElements_=NULL;
  return nElements_;
}

//#############################################################################

void
OsiIndexedVector::append(const OsiPackedVectorBase & caboose) 
{
  const int cs = caboose.getNumElements();
  
  const int * cind = caboose.getIndices();
  const double * celem = caboose.getElements();
  int maxIndex=-1;
  int i;
  for (i=0;i<cs;i++) {
    int indexValue = cind[i];
    if (indexValue<0)
      throw CoinError("negative index", "append", "OsiIndexedVector");
    if (maxIndex<indexValue)
      maxIndex = indexValue;
  }
  reserve(maxIndex+1);
  bool needClean=false;
  int numberDuplicates=0;
  for (i=0;i<cs;i++) {
    int indexValue=cind[i];
    double value=celem[i];
    if (elements_[indexValue]) {
      numberDuplicates++;
      value += elements_[indexValue];
      if (fabs(value)<TINY_ELEMENT) 
	needClean=true; // need to go through again
      elements_[indexValue]=value;
    } else {
      if (fabs(value)>=TINY_ELEMENT) {
	elements_[indexValue]=value;
	indices_[nElements_++]=indexValue;
      }
    }
  }
  if (needClean) {
    // go through again
    int size=nElements_;
    nElements_=0;
    for (i=0;i<size;i++) {
      int indexValue=indices_[i];
      double value=elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	indices_[nElements_++]=indexValue;
      } else {
        elements_[indexValue]=0.0;
      }
    }
  }
  // and clear index set
  clearIndexSet();
  delete [] packedElements_;
  packedElements_=NULL;
  if (numberDuplicates &&  testForDuplicateIndex())
    throw CoinError("duplicate index", "append", "OsiIndexedVector");
}

//#############################################################################

void
OsiIndexedVector::swap(int i, int j) 
{
  if ( i >= nElements_ ) 
    throw CoinError("index i >= size()","swap","OsiIndexedVector");
  if ( i < 0 ) 
    throw CoinError("index i < 0" ,"swap","OsiIndexedVector");
  if ( j >= nElements_ ) 
    throw CoinError("index j >= size()","swap","OsiIndexedVector");
  if ( j < 0 ) 
    throw CoinError("index j < 0" ,"swap","OsiIndexedVector");
  
  // Swap positions i and j of the
  // indices array
  
  int isave = indices_[i];
  indices_[i] = indices_[j];
  indices_[j] = isave;
}

//#############################################################################

void
OsiIndexedVector::truncate( int n ) 
{
  reserve(n);
}

//#############################################################################

void
OsiIndexedVector::operator+=(double value) 
{
  int i,indexValue;
  for (i=0;i<nElements_;i++) {
    indexValue = indices_[i];
    elements_[indexValue] += value;
  }
}

//-----------------------------------------------------------------------------

void
OsiIndexedVector::operator-=(double value) 
{
  int i,indexValue;
  for (i=0;i<nElements_;i++) {
    indexValue = indices_[i];
    elements_[indexValue] -= value;
  }
}

//-----------------------------------------------------------------------------

void
OsiIndexedVector::operator*=(double value) 
{
  int i,indexValue;
  for (i=0;i<nElements_;i++) {
    indexValue = indices_[i];
    elements_[indexValue] *= value;
  }
}

//-----------------------------------------------------------------------------

void
OsiIndexedVector::operator/=(double value) 
{
  int i,indexValue;
  for (i=0;i<nElements_;i++) {
    indexValue = indices_[i];
    elements_[indexValue] /= value;
  }
}
//#############################################################################

void
OsiIndexedVector::reserve(int n) 
{
  int i;
  // don't make allocated space smaller but do take off values
  if ( n < capacity_ ) {
    if (n<0) 
      throw CoinError("negative capacity", "reserve", "OsiIndexedVector");
    
    int nNew=0;
    for (i=0;i<nElements_;i++) {
      int indexValue=indices_[i];
      if (indexValue<n) {
        indices_[nNew++]=indexValue;
      } else {
        elements_[indexValue]=0.0;
      }
    }
    nElements_=nNew;
  } else if (n>capacity_) {
    
    // save pointers to existing data
    int * tempIndices = indices_;
    double * tempElements = elements_;
    
    // allocate new space
    indices_ = new int [n];
    elements_ = new double [n];
    
    // copy data to new space
    // and zero out part of array
    if (nElements_ > 0) {
      CoinDisjointCopyN(tempIndices, nElements_, indices_);
      CoinDisjointCopyN(tempElements, capacity_, elements_);
      CoinFillN(elements_+capacity_,n-capacity_,0.0);
    } else {
      CoinFillN(elements_,n,0.0);
    }
    capacity_ = n;
    
    // free old data
    delete [] tempElements;
    delete [] tempIndices;
  }
  // and clear index set
  clearIndexSet();
  delete [] packedElements_;
  packedElements_=NULL;
}

//#############################################################################

OsiIndexedVector::OsiIndexedVector (bool testForDuplicateIndex) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{
  // This won't fail, the indexed vector is empty. There can't be duplicate
  // indices.
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
}

//-----------------------------------------------------------------------------

OsiIndexedVector::OsiIndexedVector(int size,
                                   const int * inds, const double * elems,
                                   bool testForDuplicateIndex) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{
  gutsOfSetVector(size, inds, elems, testForDuplicateIndex,
    "constructor for array value");
}

//-----------------------------------------------------------------------------

OsiIndexedVector::OsiIndexedVector(int size,
                                   const int * inds, double value,
                                   bool testForDuplicateIndex) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{
  gutsOfSetConstant(size, inds, value, testForDuplicateIndex,
    "constructor for constant value");
}

//-----------------------------------------------------------------------------

OsiIndexedVector::OsiIndexedVector(int size, const double * element,
                                   bool testForDuplicateIndex) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{
  setFull(size, element, testForDuplicateIndex);
}

//-----------------------------------------------------------------------------

OsiIndexedVector::OsiIndexedVector(const OsiPackedVectorBase & rhs) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{  
  gutsOfSetVector(rhs.getNumElements(), rhs.getIndices(), rhs.getElements(),
    rhs.testForDuplicateIndex(), "copy constructor from base");
}

//-----------------------------------------------------------------------------

OsiIndexedVector::OsiIndexedVector(const OsiIndexedVector & rhs) :
OsiPackedVectorBase(),
indices_(NULL),
elements_(NULL),
nElements_(0),
packedElements_(NULL),
capacity_(0)
{  
  gutsOfSetVector(rhs.capacity_,rhs.nElements_, rhs.indices_, rhs.elements_,
    rhs.testForDuplicateIndex(), "copy constructor");
}

//-----------------------------------------------------------------------------

OsiIndexedVector::~OsiIndexedVector ()
{
  delete [] indices_;
  delete [] packedElements_;
  delete [] elements_;
}
//#############################################################################

// Get element values
const double * 
OsiIndexedVector::getElements() const 
{
  if (!packedElements_)
    packedElements_ = new double[nElements_];
  int i;
  for (i=0;i<nElements_;i++) {
    int indexValue=indices_[i];
    packedElements_[i]=elements_[indexValue];
  }
  return packedElements_;
}

double * 
OsiIndexedVector::getElements()  
{
  if (!packedElements_)
    packedElements_ = new double[nElements_];
  int i;
  for (i=0;i<nElements_;i++) {
    int indexValue=indices_[i];
    packedElements_[i]=elements_[indexValue];
  }
  return packedElements_;
}
//#############################################################################

/// Return the sum of two indexed vectors
OsiIndexedVector 
OsiIndexedVector::operator+(
                            const OsiIndexedVector& op2)
{
  int i;
  int nElements=nElements_;
  int capacity = CoinMax(capacity_,op2.capacity_);
  OsiIndexedVector newOne(*this);
  newOne.reserve(capacity);
  bool needClean=false;
  // new one now can hold everything so just modify old and add new
  for (i=0;i<op2.nElements_;i++) {
    int indexValue=op2.indices_[i];
    double value=op2.elements_[indexValue];
    double oldValue=elements_[indexValue];
    if (!oldValue) {
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.elements_[indexValue]=value;
	newOne.indices_[nElements++]=indexValue;
      }
    } else {
      value += oldValue;
      newOne.elements_[indexValue]=value;
      if (fabs(value)<TINY_ELEMENT) {
	needClean=true;
      }
    }
  }
  newOne.nElements_=nElements;
  if (needClean) {
    // go through again
    nElements_=0;
    for (i=0;i<nElements;i++) {
      int indexValue=newOne.indices_[i];
      double value=newOne.elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.indices_[nElements_++]=indexValue;
      } else {
        newOne.elements_[indexValue]=0.0;
      }
    }
  }
  return newOne;
}

/// Return the difference of two indexed vectors
OsiIndexedVector 
OsiIndexedVector::operator-(
                            const OsiIndexedVector& op2)
{
  int i;
  int nElements=nElements_;
  int capacity = CoinMax(capacity_,op2.capacity_);
  OsiIndexedVector newOne(*this);
  newOne.reserve(capacity);
  bool needClean=false;
  // new one now can hold everything so just modify old and add new
  for (i=0;i<op2.nElements_;i++) {
    int indexValue=op2.indices_[i];
    double value=op2.elements_[indexValue];
    double oldValue=elements_[indexValue];
    if (!oldValue) {
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.elements_[indexValue]=-value;
	newOne.indices_[nElements++]=indexValue;
      }
    } else {
      value = oldValue-value;
      newOne.elements_[indexValue]=value;
      if (fabs(value)<TINY_ELEMENT) {
	needClean=true;
      }
    }
  }
  newOne.nElements_=nElements;
  if (needClean) {
    // go through again
    nElements_=0;
    for (i=0;i<nElements;i++) {
      int indexValue=newOne.indices_[i];
      double value=newOne.elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.indices_[nElements_++]=indexValue;
      } else {
        newOne.elements_[indexValue]=0.0;
      }
    }
  }
  return newOne;
}

/// Return the element-wise product of two indexed vectors
OsiIndexedVector 
OsiIndexedVector::operator*(
                            const OsiIndexedVector& op2)
{
  int i;
  int nElements=nElements_;
  int capacity = CoinMax(capacity_,op2.capacity_);
  OsiIndexedVector newOne(*this);
  newOne.reserve(capacity);
  bool needClean=false;
  // new one now can hold everything so just modify old and add new
  for (i=0;i<op2.nElements_;i++) {
    int indexValue=op2.indices_[i];
    double value=op2.elements_[indexValue];
    double oldValue=elements_[indexValue];
    if (oldValue) {
      value *= oldValue;
      newOne.elements_[indexValue]=value;
      if (fabs(value)<TINY_ELEMENT) {
	needClean=true;
      }
    }
  }
  newOne.nElements_=nElements;
  if (needClean) {
    // go through again
    nElements_=0;
    for (i=0;i<nElements;i++) {
      int indexValue=newOne.indices_[i];
      double value=newOne.elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.indices_[nElements_++]=indexValue;
      } else {
        newOne.elements_[indexValue]=0.0;
      }
    }
  }
  return newOne;
}

/// Return the element-wise ratio of two indexed vectors
OsiIndexedVector 
OsiIndexedVector::operator/(
                            const OsiIndexedVector& op2) 
{
  // I am treating 0.0/0.0 as 0.0
  int i;
  int nElements=nElements_;
  int capacity = CoinMax(capacity_,op2.capacity_);
  OsiIndexedVector newOne(*this);
  newOne.reserve(capacity);
  bool needClean=false;
  // new one now can hold everything so just modify old and add new
  for (i=0;i<op2.nElements_;i++) {
    int indexValue=op2.indices_[i];
    double value=op2.elements_[indexValue];
    double oldValue=elements_[indexValue];
    if (oldValue) {
      if (!value)
        throw CoinError("zero divisor", "/", "OsiIndexedVector");
      value = oldValue/value;
      newOne.elements_[indexValue]=value;
      if (fabs(value)<TINY_ELEMENT) {
	needClean=true;
      }
    }
  }
  newOne.nElements_=nElements;
  if (needClean) {
    // go through again
    nElements_=0;
    for (i=0;i<nElements;i++) {
      int indexValue=newOne.indices_[i];
      double value=newOne.elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	newOne.indices_[nElements_++]=indexValue;
      } else {
        newOne.elements_[indexValue]=0.0;
      }
    }
  }
  return newOne;
}
//#############################################################################

void 
OsiIndexedVector::sortDecrIndex()
{ 
  double * elements = getElements();
  CoinSort_2(indices_, indices_ + nElements_, elements,
    CoinFirstGreater_2<int, double>());
}

void 
OsiIndexedVector::sortIncrElement()
{ 
  double * elements = getElements();
  CoinSort_2(elements, elements + nElements_, indices_,
    CoinFirstLess_2<double, int>());
}

void 
OsiIndexedVector::sortDecrElement()
{ 
  double * elements = getElements();
  CoinSort_2(elements, elements + nElements_, indices_,
    CoinFirstGreater_2<double, int>());
}

//#############################################################################

void
OsiIndexedVector::gutsOfSetVector(int size,
                                  const int * inds, const double * elems,
                                  bool testForDuplicateIndex,
                                  const char * method) 
{
  // we are going to do a faster test for duplicates so test base class when empty
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
  // and clear index set
  clearIndexSet();
  
  if (size<0)
    throw CoinError("negative number of indices", method, "OsiIndexedVector");
  
  // find largest
  int i;
  int maxIndex=-1;
  for (i=0;i<size;i++) {
    int indexValue = inds[i];
    if (indexValue<0)
      throw CoinError("negative index", method, "OsiIndexedVector");
    if (maxIndex<indexValue)
      maxIndex = indexValue;
  }
  reserve(maxIndex+1);
  nElements_ = 0;
  // elements_ array is all zero
  bool needClean=false;
  int numberDuplicates=0;
  for (i=0;i<size;i++) {
    int indexValue=inds[i];
    double value=elems[i];
    if (elements_[indexValue]) {
      numberDuplicates++;
      value += elements_[indexValue];
      if (fabs(value)<TINY_ELEMENT) 
	needClean=true; // need to go through again
      elements_[indexValue]=value;
    } else {
      if (fabs(value)>=TINY_ELEMENT) {
	elements_[indexValue]=value;
	indices_[nElements_++]=indexValue;
      }
    }
  }
  if (needClean) {
    // go through again
    size=nElements_;
    nElements_=0;
    for (i=0;i<size;i++) {
      int indexValue=indices_[i];
      double value=elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	indices_[nElements_++]=indexValue;
      } else {
        elements_[indexValue]=0.0;
      }
    }
  }
  if (numberDuplicates &&  testForDuplicateIndex)
    throw CoinError("duplicate index", "setVector", "OsiIndexedVector");
}

//#############################################################################

void
OsiIndexedVector::gutsOfSetVector(int size, int numberIndices, 
                                  const int * inds, const double * elems,
                                  bool testForDuplicateIndex,
                                  const char * method) 
{
  // we are not going to test for duplicates so test base class when empty
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
  // and clear index set
  clearIndexSet();
  
  int i;
  reserve(size);
  if (numberIndices<0)
    throw CoinError("negative number of indices", method, "OsiIndexedVector");
  nElements_ = 0;
  // elements_ array is all zero
  bool needClean=false;
  int numberDuplicates=0;
  for (i=0;i<numberIndices;i++) {
    int indexValue=inds[i];
    double value=elems[indexValue];
    if (indexValue<0) 
      throw CoinError("negative index", method, "OsiIndexedVector");
    else if (indexValue>=size) 
      throw CoinError("too large an index", method, "OsiIndexedVector");
    if (elements_[indexValue]) {
      numberDuplicates++;
      value += elements_[indexValue];
      if (fabs(value)<TINY_ELEMENT) 
	needClean=true; // need to go through again
      elements_[indexValue]=value;
    } else {
      if (fabs(value)>=TINY_ELEMENT) {
	elements_[indexValue]=value;
	indices_[nElements_++]=indexValue;
      }
    }
  }
  if (needClean) {
    // go through again
    size=nElements_;
    nElements_=0;
    for (i=0;i<size;i++) {
      int indexValue=indices_[i];
      double value=elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	indices_[nElements_++]=indexValue;
      } else {
        elements_[indexValue]=0.0;
      }
    }
  }
  if (numberDuplicates &&  testForDuplicateIndex)
    throw CoinError("duplicate index", "setVector", "OsiIndexedVector");
}

//-----------------------------------------------------------------------------

void
OsiIndexedVector::gutsOfSetConstant(int size,
                                    const int * inds, double value,
                                    bool testForDuplicateIndex,
                                    const char * method) 
{
  // we are going to do a faster test for duplicates so test base class when empty
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
  // and clear index set
  clearIndexSet();
  
  if (size<0)
    throw CoinError("negative number of indices", method, "OsiIndexedVector");
  
  // find largest
  int i;
  int maxIndex=-1;
  for (i=0;i<size;i++) {
    int indexValue = inds[i];
    if (indexValue<0)
      throw CoinError("negative index", method, "OsiIndexedVector");
    if (maxIndex<indexValue)
      maxIndex = indexValue;
  }
  
  reserve(maxIndex+1);
  nElements_ = 0;
  int numberDuplicates=0;
  // elements_ array is all zero
  bool needClean=false;
  for (i=0;i<size;i++) {
    int indexValue=inds[i];
    if (elements_[indexValue]) {
      numberDuplicates++;
      value += elements_[indexValue];
      if (fabs(value)<TINY_ELEMENT) 
	needClean=true; // need to go through again
      elements_[indexValue]=value;
    } else {
      if (fabs(value)>=TINY_ELEMENT) {
	elements_[indexValue]=value;
	indices_[nElements_++]=indexValue;
      }
    }
  }
  if (needClean) {
    // go through again
    size=nElements_;
    nElements_=0;
    for (i=0;i<size;i++) {
      int indexValue=indices_[i];
      double value=elements_[indexValue];
      if (fabs(value)>=TINY_ELEMENT) {
	indices_[nElements_++]=indexValue;
      } else {
        elements_[indexValue]=0.0;
      }
    }
  }
  if (numberDuplicates &&  testForDuplicateIndex)
    throw CoinError("duplicate index", "setConstant", "OsiIndexedVector");
}

//#############################################################################
