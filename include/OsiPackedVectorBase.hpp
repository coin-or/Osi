// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiPackedVectorBase_H
#define OsiPackedVectorBase_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <set>
#include <map>
#include "CoinError.hpp"
#include "OsiFloatEqual.hpp"

class OsiPackedVector;

/** Abstract base class for various sparse vectors.

    Since this class is abstract, no object of this type can be created. The
    sole purpose of this class is to provide access to a <em>constant</em>
    packed vector. All members of this class are const methods, they can't
    change the object. */

class OsiPackedVectorBase  {
  
public:
   /**@name Virtual methods that the derived classes must provide */
   //@{
   /// Get length of indices and elements vectors
   virtual int getNumElements() const = 0;
   /// Get indices of elements
   virtual const int * getIndices() const = 0;
   /// Get element values
   virtual const double * getElements() const = 0;
   //@}

   /**@name Methods related to whether a check for duplicate 
       indices will be performed.

       If the checking for duplicate indices is turned off, then
       some OsiPackedVector methods may not work correctly if there
       are duplicate indices.
       Turning off the checking for duplicate indices may result in
       better run time performance.
      */
   //@{
   /** Set to the argument value whether to test for duplicate indices in the
       vector whenever they can occur. */
   void setTestForDuplicateIndex(bool test) const;
   /** Returns true if the vector should be tested for duplicate indices when
       they can occur. */
   bool testForDuplicateIndex() const { return testForDuplicateIndex_; }
   //@}

   /**@name Methods for getting info on the packed vector as a full vector */
   //@{
   /** Get the vector as a dense vector. The argument specifies how long this
       dense vector is. <br>
       <strong>NOTE</code>: The user needs to <code>delete[]</code> this
       pointer after it's not needed anymore.
   */
   double * denseVector(int denseSize) const throw(CoinError);
   /** Access the i'th element of the full storage vector.
       If the i'th is not stored, then zero is returned. The initial use of
       this method has some computational and storage overhead associated with
       it.<br>
       <strong>NOTE</strong>: This is <em>very</em> expensive. It is probably
       much better to use <code>denseVector()</code>.
   */
   double operator[](int i) const throw(CoinError); 
   //@}

   /**@name Index methods */
   //@{
   /// Get value of maximum index
   int getMaxIndex() const;
   /// Get value of minimum index
   int getMinIndex() const;

   /// Throw an exception if there are duplicate indices
   void duplicateIndex(const char* methodName = NULL,
		       const char * className = NULL)  const throw(CoinError);

   /** Return true if the i'th element of the full storage vector exists in
       the packed storage vector.*/
   bool isExistingIndex(int i) const;
   
   /** Return the position of the i'th element of the full storage vector.
       If index does not exist then -1 is returned  */
   int findIndex(int i) const;
 
   //@}
  
   /**@name Comparison operators on two packed vectors */
   //@{
   /** Equal. Returns true if vectors have same length and corresponding
       element of each vector is equal. */
   bool operator==(const OsiPackedVectorBase & rhs) const;
   /// Not equal
   bool operator!=(const OsiPackedVectorBase & rhs) const;

   /** equivalent - If shallow packed vector A & B are equivalent, then they
       are still equivalent no matter how they are sorted.
       In this method the FloatEqual function operator can be specified. The
       default equivalence test is that the entries are relatively equal.<br> 
       <strong>NOTE</strong>: This is a relatively expensive method as it
       sorts the two shallow packed vectors.
   */
   template <class FloatEqual> bool
   isEquivalent(const OsiPackedVectorBase& rhs, const FloatEqual& eq) const
      throw(CoinError)
   {
      if (getNumElements() != rhs.getNumElements())
	 return false;

      duplicateIndex("equivalent", "OsiPackedVector");
      rhs.duplicateIndex("equivalent", "OsiPackedVector");

      std::map<int,double> mv;
      const int * inds = getIndices();
      const double * elems = getElements();
      int i;
      for ( i = getNumElements() - 1; i >= 0; --i) {
	 mv.insert(std::make_pair(inds[i], elems[i]));
      }

      std::map<int,double> mvRhs;
      inds = rhs.getIndices();
      elems = rhs.getElements();
      for ( i = getNumElements() - 1; i >= 0; --i) {
	 mvRhs.insert(std::make_pair(inds[i], elems[i]));
      }

      std::map<int,double>::const_iterator mvI = mv.begin();
      std::map<int,double>::const_iterator mvIlast = mv.end();
      std::map<int,double>::const_iterator mvIrhs = mvRhs.begin();
      while (mvI != mvIlast) {
	 if (mvI->first != mvIrhs->first || ! eq(mvI->second, mvIrhs->second))
	    return false;
	 ++mvI;
	 ++mvIrhs;
      }
      return true;
   };
   bool isEquivalent(const OsiPackedVectorBase& rhs) const
   {
      return isEquivalent(rhs,  OsiRelFltEq());
   }
   //@}


   /**@name Arithmetic operators. */
   //@{
   /// Create the dot product with a full vector
   double dotProduct(const double* dense) const;
   /// Sum elements of vector.
   double sum() const;

   template <class BinaryFunction> OsiPackedVector
   binaryOp(double value, BinaryFunction bf) const
   {
      OsiPackedVector retVal;
      retVal.setTestForDuplicateIndex(true);

      const int s = getNumElements();
      if (s > 0) {
	 retVal.reserve(s);
	 const int * inds = getIndices();
	 const double * elems = getElements();
	 for (int i=0; i<s; ++i ) {
	    retVal.insert(inds[i], bf(value, elems[i]));
	 }
      }
      return retVal;
   }

   template <class BinaryFunction> OsiPackedVector
   binaryOp(BinaryFunction bf, double value) const
   {
      OsiPackedVector retVal;
      retVal.setTestForDuplicateIndex(true);

      const int s = getNumElements();
      if (s > 0) {
	 retVal.reserve(s);
	 const int * inds = getIndices();
	 const double * elems = getElements();
	 for (int i=0; i<s; ++i ) {
	    retVal.insert(inds[i], bf(elems[i], value));
	 }
      }
      return retVal;
   }

   template <class BinaryFunction> OsiPackedVector
   binaryOp(const OsiPackedVectorBase& op2, BinaryFunction bf) const
   {
      OsiPackedVector retVal;
      retVal.setTestForDuplicateIndex(true);
      if (op2.getNumElements() == 0)
	 return retVal;

      // loop once for each element in *this
      const int s1 = getNumElements();
      const int * inds1 = getIndices();
      const double * elems1 = getElements();
      const int maxind2 = op2.getMaxIndex();
      double * full2p = op2.denseVector(maxind2 + 1);
      const double *full2 = full2p;
      int i;
      for ( i=0; i<s1; ++i ) {
	 const int index = inds1[i];
	 const double val = bf(elems1[i], index>maxind2 ? 0.0 : full2[index]);
	 //	 if (val != 0.0) // *THINK* : should we put in only nonzeros?
	    retVal.insert(index, val);
      }
      delete[] full2p;

      const int s2 = op2.getNumElements();
      const int * inds2 = op2.getIndices();
      const double * elems2 = op2.getElements();

      // loop once for each element in operand2  
      for ( i=0; i<s2; ++i ) {
	 const int index = inds2[i];
	 // if index exists in *this, then element was processed in prior loop
	 if ( isExistingIndex(index) )
	    continue;
	 // Index does not exist in *this, so the element value must be zero
	 const double val = bf(0.0, elems2[i]);
	 //	 if (val != 0.0) // *THINK* : should we put in only nonzeros?
	    retVal.insert(index, val);
      }

      return retVal;
   }
   //@}
  
protected:

   /**@name Constructors, destructor<br>
      <strong>NOTE</strong>: All constructors are protected. There's no need
      to expose them, after all, this is an abstract class. */
   //@{
   /** Default constructor. */
   OsiPackedVectorBase();
   /** Destructor */
   virtual ~OsiPackedVectorBase();
   //@}

private:
   /**@name Disabled methods */
   //@{
   /** The copy constructor. <br>
       This must be at least protected, but we make it private. The reason is
       that when, say, a shallow packed vector is created, first the
       underlying class, it this one is constructed. However, at that point we
       don't know how much of the data members of this class we need to copy
       over. Therefore the copy constructor is not used. */
   OsiPackedVectorBase(const OsiPackedVectorBase&);
   /** This class provides <em>const</em> access to packed vectors, so there's
       no need to provide an assignment operator. */
   OsiPackedVectorBase& operator=(const OsiPackedVectorBase&);
   //@}
   
protected:
    
   /**@name Protected methods */
   //@{      
   /// Find Maximum and Minimum Indices
   void findMaxMinIndices() const;

   /// Return indexSetPtr_ (create it if necessary).
   std::set<int> * indexSet(const char* methodName = NULL,
			    const char * className = NULL) const
      throw(CoinError);

   /// Delete the indexSet
   void clearIndexSet() const;
   void clearBase() const;
   void copyMaxMinIndex(const OsiPackedVectorBase & x) const {
      maxIndex_ = x.maxIndex_;
      minIndex_ = x.minIndex_;
   }
   //@}
    
private:
   /**@name Protected member data */
   //@{
   /// Contains max index value or -infinity
   mutable int maxIndex_;
   /// Contains minimum index value or infinity
   mutable int minIndex_;
   /** Store the indices in a set. This set is only created if it is needed.
       Its primary use is testing for duplicate indices.
    */
   mutable std::set<int> * indexSetPtr_;
   /** True if the vector should be tested for duplicate indices when they can
       occur. */
   mutable bool testForDuplicateIndex_;
   /** True if the vector has already been tested for duplicate indices. Most
       of the operations in OsiPackedVector preserves this flag. */
   mutable bool testedDuplicateIndex_;
   //@}
};

#endif
