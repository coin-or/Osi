// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinHelperFunctions.hpp"
#include "OsiShallowPackedVector.hpp"

//#############################################################################

void
OsiShallowPackedVector::clear()
{
   clearBase();
   indices_ = NULL;
   elements_ = NULL;
   nElements_ = 0;
}

//#############################################################################

OsiShallowPackedVector&
OsiShallowPackedVector::operator=(const OsiPackedVectorBase & x)
{
   if (&x != this) {
      indices_ = x.getIndices();
      elements_ = x.getElements();
      nElements_ = x.getNumElements();
      OsiPackedVectorBase::clearBase();
      OsiPackedVectorBase::copyMaxMinIndex(x);
      try {
	 OsiPackedVectorBase::duplicateIndex();
      }
      catch (CoinError e) {
	 throw CoinError("duplicate index", "operator= from base",
			"OsiShallowPackedVector");
      }
   }
   return *this;
}

//#############################################################################

OsiShallowPackedVector&
OsiShallowPackedVector::operator=(const OsiShallowPackedVector & x)
{
   if (&x != this) {
      indices_ = x.indices_;
      elements_ = x.elements_;
      nElements_ = x.nElements_;
      OsiPackedVectorBase::clearBase();
      OsiPackedVectorBase::copyMaxMinIndex(x);
      try {
	 OsiPackedVectorBase::duplicateIndex();
      }
      catch (CoinError e) {
	 throw CoinError("duplicate index", "operator=",
			"OsiShallowPackedVector");
      }
   }
   return *this;
}

//#############################################################################

void
OsiShallowPackedVector::setVector(int size,
				  const int * inds, const double * elems,
				  bool testForDuplicateIndex)
   throw(CoinError)
{
   indices_ = inds;
   elements_ = elems;
   nElements_ = size;
   OsiPackedVectorBase::clearBase();
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "setVector",
		     "OsiShallowPackedVector");
   }
}

//#############################################################################

//-------------------------------------------------------------------
// Default
//-------------------------------------------------------------------
OsiShallowPackedVector::OsiShallowPackedVector(bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(NULL),
   elements_(NULL),
   nElements_(0)
{
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "default constructor",
		     "OsiShallowPackedVector");
   }
}
   
//-------------------------------------------------------------------
// Explicit
//-------------------------------------------------------------------
OsiShallowPackedVector::OsiShallowPackedVector(int size, 
					       const int * inds,
					       const double * elems,
					       bool testForDuplicateIndex) :
   OsiPackedVectorBase(),
   indices_(inds),
   elements_(elems),
   nElements_(size)
{
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(testForDuplicateIndex);
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "explicit constructor",
		     "OsiShallowPackedVector");
   }
}

//-------------------------------------------------------------------
// Copy
//-------------------------------------------------------------------
OsiShallowPackedVector::OsiShallowPackedVector(const OsiPackedVectorBase& x) :
   OsiPackedVectorBase(),
   indices_(x.getIndices()),
   elements_(x.getElements()),
   nElements_(x.getNumElements())
{
   OsiPackedVectorBase::copyMaxMinIndex(x);
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(x.testForDuplicateIndex());
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "copy constructor from base",
		     "OsiShallowPackedVector");
   }
}

//-------------------------------------------------------------------
// Copy
//-------------------------------------------------------------------
OsiShallowPackedVector::OsiShallowPackedVector(const
					       OsiShallowPackedVector& x) :
   OsiPackedVectorBase(),
   indices_(x.getIndices()),
   elements_(x.getElements()),
   nElements_(x.getNumElements())
{
   OsiPackedVectorBase::copyMaxMinIndex(x);
   try {
      OsiPackedVectorBase::setTestForDuplicateIndex(x.testForDuplicateIndex());
   }
   catch (CoinError e) {
      throw CoinError("duplicate index", "copy constructor",
		     "OsiShallowPackedVector");
   }
}
