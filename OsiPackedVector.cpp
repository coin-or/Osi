// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "CoinHelperFunctions.hpp"
#include "OsiPackedVector.hpp"

//#############################################################################

void
OsiPackedVector::clear()
{
   nElements_ = 0;
   clearBase();
}

//#############################################################################

OsiPackedVector &
OsiPackedVector::operator=(const OsiPackedVector & rhs)
{
   if (this != &rhs) {
      clear();
      gutsOfSetVector(rhs.getNumElements(), rhs.getIndices(), rhs.getElements(),
		      OsiPackedVectorBase::testForDuplicateIndex(),
		      "operator=");
   }
   return *this;
}

//#############################################################################

OsiPackedVector &
OsiPackedVector::operator=(const OsiPackedVectorBase & rhs)
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
#if 0
void
OsiPackedVector::assignVector(int size, int*& inds, double*& elems,
			      bool testForDuplicateIndex)
{
   clear();
   // Allocate storage
   if ( size != 0 ) {
      reserve(size);
      nElements_ = size;
      indices_ = inds;    inds = NULL;
      elements_ = elems;  elems = NULL;
      CoinIotaN(origIndices_, size, 0);
   }
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "assignVector", "OsiPackedVector");
   }
}
#else
void
OsiPackedVector::assignVector(int size, int*& inds, double*& elems,
                              bool testForDuplicateIndex)
{
  clear();
  // Allocate storage
  if ( size != 0 ) {
		  //reserve(size); //This is a BUG!!!
    nElements_ = size;
    if (indices_ != NULL) delete[] indices_;
    indices_ = inds;    inds = NULL;
    if (elements_ != NULL) delete[] elements_;
    elements_ = elems;  elems = NULL;
    if (origIndices_ != NULL) delete[] origIndices_;
    origIndices_ = new int[size];
    CoinIotaN(origIndices_, size, 0);
    capacity_ = size;
  }
  try {    
    OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
  }
  catch (CoinError e) {
    throw CoinError("duplicate index", "assignVector",
      "OsiPackedVector");
  }
}
#endif

//#############################################################################

void
OsiPackedVector::setVector(int size, const int * inds, const double * elems,
			   bool testForDuplicateIndex)
   throw(CoinError)
{
   clear();
   gutsOfSetVector(size, inds, elems, testForDuplicateIndex, "setVector");
}

//#############################################################################

void
OsiPackedVector::setConstant(int size, const int * inds, double value,
			     bool testForDuplicateIndex)
   throw(CoinError)
{
   clear();
   gutsOfSetConstant(size, inds, value, testForDuplicateIndex, "setConstant");
}

//#############################################################################

void
OsiPackedVector::setFull(int size, const double * elems,
			 bool testForDuplicateIndex) 
{
  // Clear out any values presently stored
  clear();
  
  // Allocate storage
  if ( size!=0 ) {
    reserve(size);  
    nElements_ = size;

    CoinIotaN(origIndices_, size, 0);
    CoinIotaN(indices_, size, 0);
    CoinDisjointCopyN(elems, size, elements_);
  }
  OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
}

//#############################################################################

void
OsiPackedVector::setElement(int index, double element) throw(CoinError)
{
   if ( index >= nElements_ ) 
      throw CoinError("index >= size()", "setElement", "OsiPackedVector");
   if ( index < 0 ) 
      throw CoinError("index < 0" , "setElement", "OsiPackedVector");
   elements_[index] = element;
}

//#############################################################################

void
OsiPackedVector::insert( int index, double element ) throw(CoinError)
{
   const int s = nElements_;
   if (testForDuplicateIndex()) {
      std::set<int>& is = *indexSet("insert", "OsiPackedVector");
      if (! is.insert(index).second)
	 throw CoinError("Index already exists", "insert", "OsiPackedVector");
   }

   if( capacity_ <= s ) {
      reserve( CoinMax(5, 2*capacity_) );
      assert( capacity_ > s );
   }
   indices_[s] = index;
   elements_[s] = element;
   origIndices_[s] = s;
   ++nElements_;
}

//#############################################################################

void
OsiPackedVector::append(const OsiPackedVectorBase & caboose) throw(CoinError)
{
   const int s = nElements_;
   const int cs = caboose.getNumElements();
   // Make sure there is enough room for the caboose
   if ( capacity_ < s + cs)
      reserve(CoinMax(s + cs, 2 * capacity_));

   const int * cind = caboose.getIndices();
   const double * celem = caboose.getElements();
   CoinDisjointCopyN(cind, cs, indices_ + s);
   CoinDisjointCopyN(celem, cs, elements_ + s);
   CoinIotaN(origIndices_ + s, cs, s);
   nElements_ += cs;
   if (testForDuplicateIndex()) {
      std::set<int>& is = *indexSet("append", "OsiPackedVector");
      for (int i = 0; i < cs; ++i) {
	 if (!is.insert(cind[i]).second)
	    throw CoinError("duplicate index", "append", "OsiPackedVector");
      }
   }
}

//#############################################################################

void
OsiPackedVector::swap(int i, int j) throw(CoinError)
{
   if ( i >= nElements_ ) 
      throw CoinError("index i >= size()","swap","OsiPackedVector");
   if ( i < 0 ) 
      throw CoinError("index i < 0" ,"swap","OsiPackedVector");
   if ( i >= nElements_ ) 
      throw CoinError("index j >= size()","swap","OsiPackedVector");
   if ( i < 0 ) 
      throw CoinError("index j < 0" ,"swap","OsiPackedVector");

   // Swap positions i and j of the
   // indices and elements arrays
   std::swap(indices_[i], indices_[j]);
   std::swap(elements_[i], elements_[j]);
}

//#############################################################################

void
OsiPackedVector::truncate( int n ) throw(CoinError)
{
   if ( n > nElements_ ) 
      throw CoinError("n > size()","truncate","OsiPackedVector");
   if ( n < 0 ) 
      throw CoinError("n < 0","truncate","OsiPackedVector");
   nElements_ = n;
   clearBase();
}

//#############################################################################

void
OsiPackedVector::operator+=(double value) 
{
   std::transform(elements_, elements_ + nElements_, elements_,
		  std::bind2nd(std::plus<double>(), value) );
}

//-----------------------------------------------------------------------------

void
OsiPackedVector::operator-=(double value) 
{
   std::transform(elements_, elements_ + nElements_, elements_,
		  std::bind2nd(std::minus<double>(), value) );
}

//-----------------------------------------------------------------------------

void
OsiPackedVector::operator*=(double value) 
{
   std::transform(elements_, elements_ + nElements_, elements_,
		  std::bind2nd(std::multiplies<double>(), value) );
}

//-----------------------------------------------------------------------------

void
OsiPackedVector::operator/=(double value) 
{
   std::transform(elements_, elements_ + nElements_, elements_,
		  std::bind2nd(std::divides<double>(), value) );
}

//#############################################################################

void
OsiPackedVector::sortOriginalOrder() {
  CoinSort_3(origIndices_, origIndices_ + nElements_, indices_, elements_);
}

//#############################################################################

void
OsiPackedVector::reserve(int n)
{
   // don't make allocated space smaller
   if ( n <= capacity_ )
      return;
   capacity_ = n;

   // save pointers to existing data
   int * tempIndices = indices_;
   int * tempOrigIndices = origIndices_;
   double * tempElements = elements_;

   // allocate new space
   indices_ = new int [capacity_];
   origIndices_ = new int [capacity_];
   elements_ = new double [capacity_];

   // copy data to new space
   if (nElements_ > 0) {
      CoinDisjointCopyN(tempIndices, nElements_, indices_);
      CoinDisjointCopyN(tempOrigIndices, nElements_, origIndices_);
      CoinDisjointCopyN(tempElements, nElements_, elements_);
   }

   // free old data
   delete [] tempElements;
   delete [] tempOrigIndices;
   delete [] tempIndices;
}

//#############################################################################

OsiPackedVector::OsiPackedVector (bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{
   // This won't fail, the packed vector is empty. There can't be duplicate
   // indices.
   OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
}

//-----------------------------------------------------------------------------

OsiPackedVector::OsiPackedVector(int size,
				 const int * inds, const double * elems,
				 bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{
   gutsOfSetVector(size, inds, elems, testForDuplicateIndex,
		   "constructor for array value");
}

//-----------------------------------------------------------------------------

OsiPackedVector::OsiPackedVector(int size,
				 const int * inds, double value,
				 bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{
   gutsOfSetConstant(size, inds, value, testForDuplicateIndex,
		     "constructor for constant value");
}

//-----------------------------------------------------------------------------

OsiPackedVector::OsiPackedVector(int size, const double * element,
				 bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{
   setFull(size, element, testForDuplicateIndex);
}

//-----------------------------------------------------------------------------

OsiPackedVector::OsiPackedVector(const OsiPackedVectorBase & rhs) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{  
   gutsOfSetVector(rhs.getNumElements(), rhs.getIndices(), rhs.getElements(),
		   rhs.testForDuplicateIndex(), "copy constructor from base");
}

//-----------------------------------------------------------------------------

OsiPackedVector::OsiPackedVector(const OsiPackedVector & rhs) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0),
   origIndices_(NULL),
   capacity_(0)
{  
   gutsOfSetVector(rhs.getNumElements(), rhs.getIndices(), rhs.getElements(),
		   rhs.testForDuplicateIndex(), "copy constructor");
}

//-----------------------------------------------------------------------------

OsiPackedVector::~OsiPackedVector ()
{
   delete [] indices_;
   delete [] origIndices_;
   delete [] elements_;
}

//#############################################################################

void
OsiPackedVector::gutsOfSetVector(int size,
				 const int * inds, const double * elems,
				 bool testForDuplicateIndex,
				 const char * method)
{
   if ( size != 0 ) {
      reserve(size);
      nElements_ = size;
      CoinDisjointCopyN(inds, size, indices_);
      CoinDisjointCopyN(elems, size, elements_);
      CoinIotaN(origIndices_, size, 0);
   }
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", method, "OsiPackedVector");
   }
}

//-----------------------------------------------------------------------------

void
OsiPackedVector::gutsOfSetConstant(int size,
				   const int * inds, double value,
				   bool testForDuplicateIndex,
				   const char * method)
{
   if ( size != 0 ) {
      reserve(size);
      nElements_ = size;
      CoinDisjointCopyN(inds, size, indices_);
      CoinFillN(elements_, size, value);
      CoinIotaN(origIndices_, size, 0);
   }
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", method, "OsiPackedVector");
   }
}

//#############################################################################
