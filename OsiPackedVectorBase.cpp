// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <numeric>

#include "CoinHelperFunctions.hpp"
#include "OsiPackedVectorBase.hpp"

//#############################################################################

double *
OsiPackedVectorBase::denseVector(int denseSize) const throw(CoinError)
{
   if (getMaxIndex() >= denseSize)
      throw CoinError("Dense vector size is less than max index",
		     "denseVector", "OsiPackedVectorBase");

   double * dv = new double[denseSize];
   CoinFillN(dv, denseSize, 0.0);
   const int s = getNumElements();
   const int * inds = getIndices();
   const double * elems = getElements();
   for (int i = 0; i < s; ++i)
      dv[inds[i]] = elems[i];
   return dv;
}

//-----------------------------------------------------------------------------

double
OsiPackedVectorBase::operator[](int i) const throw (CoinError)
{
   if (! testedDuplicateIndex_)
      duplicateIndex("operator[]", "OsiPackedVectorBase");

   // Get a reference to a map of full storage indices to 
   // packed storage location.
   const std::set<int> & sv = *indexSet("operator[]", "OsiPackedVectorBase");
#if 1
   if (sv.find(i) == sv.end())
      return 0.0;
   return getElements()[findIndex(i)];
#else
   // LL: suggested change, somthing is wrong with this
   const size_t ind = std::distance(sv.begin(), sv.find(i));
   return (ind == sv.size()) ? 0.0 : getElements()[ind];
#endif
      
}

//#############################################################################

void
OsiPackedVectorBase::setTestForDuplicateIndex(bool test) const
{
   if (test == true) {
      testForDuplicateIndex_ = true;
      duplicateIndex("setTestForDuplicateIndex", "OsiPackedVectorBase");
   } else {
      testForDuplicateIndex_ = false;
      testedDuplicateIndex_ = false;
   }
}

//#############################################################################

int
OsiPackedVectorBase::getMaxIndex() const
{
   findMaxMinIndices();
   return maxIndex_;
}

//-----------------------------------------------------------------------------

int
OsiPackedVectorBase::getMinIndex() const
{
   findMaxMinIndices();
   return minIndex_;
}

//-----------------------------------------------------------------------------

void
OsiPackedVectorBase::duplicateIndex(const char* methodName,
				    const char * className)
   const throw(CoinError)
{
   if (testForDuplicateIndex())
      indexSet(methodName, className);
   testedDuplicateIndex_ = true;
}

//-----------------------------------------------------------------------------

bool
OsiPackedVectorBase::isExistingIndex(int i) const
{
   if (! testedDuplicateIndex_)
      duplicateIndex("indexExists", "OsiPackedVectorBase");

   const std::set<int> & sv = *indexSet("indexExists", "OsiPackedVectorBase");
   return sv.find(i) != sv.end();
}


int
OsiPackedVectorBase::findIndex(int i) const
{   
   const int * inds = getIndices();
   int retVal = std::find(inds, inds + getNumElements(), i) - inds;
   if (retVal == getNumElements() ) retVal = -1;
   return retVal;
}

//#############################################################################

bool
OsiPackedVectorBase::operator==(const OsiPackedVectorBase& rhs) const
{
   return (getNumElements()==rhs.getNumElements() &&
	   std::equal(getIndices(), getIndices() + getNumElements(),
		      rhs.getIndices()) &&
	   std::equal(getElements(), getElements() + getNumElements(),
		      rhs.getElements()));
}

//-----------------------------------------------------------------------------

bool
OsiPackedVectorBase::operator!=(const OsiPackedVectorBase& rhs) const
{
   return !( (*this)==rhs );
}

//-----------------------------------------------------------------------------

double
OsiPackedVectorBase::dotProduct(const double* dense) const
{
   const double * elems = getElements();
   const int * inds = getIndices();
   double dp = 0.0;
   for (int i = getNumElements() - 1; i >= 0; --i)
      dp += elems[i] * dense[inds[i]];
   return dp;
}

//-----------------------------------------------------------------------------

double
OsiPackedVectorBase::sum() const
{
   return std::accumulate(getElements(), getElements() + getNumElements(), 0.0);
}

//#############################################################################

OsiPackedVectorBase::OsiPackedVectorBase() :
   maxIndex_(/*std::numeric_limits<int>::max()*/INT_MIN/*0*/),
   minIndex_(/*std::numeric_limits<int>::min()*/INT_MAX/*0*/),
   indexSetPtr_(NULL),
   testForDuplicateIndex_(true),
   testedDuplicateIndex_(false) {}

//-----------------------------------------------------------------------------

OsiPackedVectorBase::~OsiPackedVectorBase()
{
   delete indexSetPtr_;
}

//#############################################################################
//#############################################################################

void
OsiPackedVectorBase::findMaxMinIndices() const
{
   if ( getNumElements()==0 ) 
      return;
   // if indexSet exists then grab begin and rend to get min & max indices
   else if ( indexSetPtr_ != NULL ) {
      maxIndex_ = *indexSetPtr_->rbegin();
      minIndex_ = *indexSetPtr_-> begin();
   } else {
      // Have to scan through vector to find min and max.
      maxIndex_ = *(std::max_element(getIndices(),
				     getIndices() + getNumElements()));
      minIndex_ = *(std::min_element(getIndices(),
				     getIndices() + getNumElements()));
   }
}

//-------------------------------------------------------------------

std::set<int> *
OsiPackedVectorBase::indexSet(const char* methodName,
			      const char * className) const throw(CoinError)
{
   testedDuplicateIndex_ = true;
   if ( indexSetPtr_ == NULL ) {
      // create a set of the indices
      indexSetPtr_ = new std::set<int>;
      const int s = getNumElements();
      const int * inds = getIndices();
      for (int j=0; j < s; ++j) {
	 if (!indexSetPtr_->insert(inds[j]).second) {
	    testedDuplicateIndex_ = false;
	    delete indexSetPtr_;
	    indexSetPtr_ = NULL;
	    if (methodName != NULL) {
	       throw CoinError("Duplicate index found", methodName, className);
	    } else {
	       throw CoinError("Duplicate index found",
			      "indexSet", "OsiPackedVectorBase");
	    }
	 }
      }
   }
   return indexSetPtr_;
}

//-----------------------------------------------------------------------------

void
OsiPackedVectorBase::clearIndexSet() const
{
   delete indexSetPtr_;
   indexSetPtr_ = NULL;
}

//-----------------------------------------------------------------------------

void
OsiPackedVectorBase::clearBase() const
{
   clearIndexSet();
   maxIndex_ = /*std::numeric_limits<int>::max()*/INT_MIN/*0*/;
   minIndex_ = /*std::numeric_limits<int>::min()*/INT_MAX/*0*/;
   testedDuplicateIndex_ = false;
}

//#############################################################################
