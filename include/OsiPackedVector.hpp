// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiPackedVector_H
#define OsiPackedVector_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <map>

#include "OsiPackedVectorBase.hpp"
#include "CoinSort.hpp"

/** Sparse Vector

Stores vector of indices and associated element values.
Supports sorting of vector while maintaining the original indices.

Here is a sample usage:
<pre>
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };

    // Create vector and set its value
    OsiPackedVector r(ne,inx,el);

    // access each index and element
    assert( r.indices ()[0]== 1  );
    assert( r.elements()[0]==10. );
    assert( r.indices ()[1]== 4  );
    assert( r.elements()[1]==40. );
    assert( r.indices ()[2]== 0  );
    assert( r.elements()[2]== 1. );
    assert( r.indices ()[3]== 2  );
    assert( r.elements()[3]==50. );

    // access original position of index
    assert( r.originalPosition()[0]==0 );
    assert( r.originalPosition()[1]==1 );
    assert( r.originalPosition()[2]==2 );
    assert( r.originalPosition()[3]==3 );

    // access as a full storage vector
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);

    // sort Elements in increasing order
    r.sortIncrElement();

    // access each index and element
    assert( r.indices ()[0]== 0  );
    assert( r.elements()[0]== 1. );
    assert( r.indices ()[1]== 1  );
    assert( r.elements()[1]==10. );
    assert( r.indices ()[2]== 4  );
    assert( r.elements()[2]==40. );
    assert( r.indices ()[3]== 2  );
    assert( r.elements()[3]==50. );    

    // access original position of index    
    assert( r.originalPosition()[0]==2 );
    assert( r.originalPosition()[1]==0 );
    assert( r.originalPosition()[2]==1 );
    assert( r.originalPosition()[3]==3 );

    // access as a full storage vector
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);

    // Restore orignal sort order
    r.sortOriginalOrder();
    
    assert( r.indices ()[0]== 1  );
    assert( r.elements()[0]==10. );
    assert( r.indices ()[1]== 4  );
    assert( r.elements()[1]==40. );
    assert( r.indices ()[2]== 0  );
    assert( r.elements()[2]== 1. );
    assert( r.indices ()[3]== 2  );
    assert( r.elements()[3]==50. );

    // Tests for equality and equivalence
    OsiPackedVector r1;
    r1=r;
    assert( r==r1 );
    assert( r.equivalent(r1) );
    r.sortIncrElement();
    assert( r!=r1 );
    assert( r.equivalent(r1) );

    // Add packed vectors.
    // Similarly for subtraction, multiplication,
    // and division.
    OsiPackedVector add = r + r1;
    assert( add[0] ==  1.+ 1. );
    assert( add[1] == 10.+10. );
    assert( add[2] == 50.+50. );
    assert( add[3] ==  0.+ 0. );
    assert( add[4] == 40.+40. );

    assert( r.sum() == 10.+40.+1.+50. );
</pre>
*/
class OsiPackedVector : public OsiPackedVectorBase {
   friend void OsiPackedVectorUnitTest();
  
public:
   /**@name Get methods. */
   //@{
   /// Get the size
   virtual int getNumElements() const { return nElements_; }
   /// Get indices of elements
   virtual const int * getIndices() const { return indices_; }
   /// Get element values
   virtual const double * getElements() const { return elements_; }
   /// Get indices of elements
   int * getIndices() { return indices_; }
   /// Get element values
   double * getElements() { return elements_; }
   /** Get pointer to int * vector of original postions.
       If the packed vector has not been sorted then this
       function returns the vector: 0, 1, 2, ..., size()-1. */
   const int * getOriginalPosition() const { return origIndices_; }
   //@}
 
   //-------------------------------------------------------------------
   // Set indices and elements
   //------------------------------------------------------------------- 
   /**@name Set methods */
   //@{
   /// Reset the vector (as if were just created an empty vector)
   void clear();
   /** Assignment operator. <br>
       <strong>NOTE</strong>: This operator keeps the current
       <code>testForDuplicateIndex</code> setting, and affter copying the data
       it acts accordingly. */
   OsiPackedVector & operator=(const OsiPackedVector &);
   /** Assignment operator <em>for a PackedVectorBase</em>. <br>
       <strong>NOTE</strong>: This operator keeps the current
       <code>testForDuplicateIndex</code> setting, and affter copying the data
       it acts accordingly. */
   OsiPackedVector & operator=(const OsiPackedVectorBase & rhs);

   /** Assign the ownership of the arguments to this vector.
       Size is the length of both the indices and elements vectors.
       The indices and elements vectors are copied into this class instance's
       member data. The last argument indicates whether this vector will have
       to be tested for duplicate indices.
   */
   void assignVector(int size, int*& inds, double*& elems,
		     bool testForDuplicateIndex = true);

   /** Set vector size, indices, and elements.
       Size is the length of both the indices and elements vectors.
       The indices and elements vectors are copied into this class instance's
       member data. The last argument specifies whether this vector will have
       to be checked for duplicate indices whenever that can happen. */
   void setVector(int size, const int * inds, const double * elems,
		  bool testForDuplicateIndex = true) throw(CoinError);
  
   /** Elements set to have the same scalar value */
   void setConstant(int size, const int * inds, double elems,
		    bool testForDuplicateIndex = true) throw(CoinError);
  
   /** Indices are not specified and are taken to be 0,1,...,size-1 */
   void setFull(int size, const double * elems,
		bool testForDuplicateIndex = true);

   /** Set an existing element in the packed vector
       The first argument is the "index" into the elements() array
   */
   void setElement(int index, double element) throw(CoinError);

   /// Insert an element into the vector
   void insert(int index, double element) throw(CoinError);
   /// Append a OsiPackedVector to the end
   void append(const OsiPackedVectorBase & caboose) throw(CoinError);

   /// Swap values in positions i and j of indices and elements
   void swap(int i, int j) throw(CoinError); 

   /** Resize the packed vector to be the first newSize elements.
       Problem with truncate: what happens with origIndices_ ??? */
   void truncate(int newSize) throw(CoinError); 
   //@}

   /**@name Arithmetic operators. */
   //@{
   /// add <code>value</code> to every entry
   void operator+=(double value);
   /// subtract <code>value</code> from every entry
   void operator-=(double value);
   /// multiply every entry by <code>value</code>
   void operator*=(double value);
   /// divide every entry by <code>value</code>
   void operator/=(double value);
   //@}

   /**@name Sorting */
   //@{ 
   /** Sort the packed storage vector.
       Typcical usages:
       <pre> 
       packedVector.sort(OsiIncrIndexOrdered());   //increasing indices
       packedVector.sort(OsiIncrElementOrdered()); // increasing elements
       </pre>
   */ 
   template <class OsiCompare3>
   void sort(const OsiCompare3 & tc)
   { CoinSort_3(indices_, indices_ + nElements_, origIndices_, elements_,
		tc); }

   void sortIncrIndex()
   { CoinSort_3(indices_, indices_ + nElements_, origIndices_, elements_,
		CoinFirstLess_3<int, int, double>()); }

   void sortDecrIndex()
   { CoinSort_3(indices_, indices_ + nElements_, origIndices_, elements_,
		CoinFirstGreater_3<int, int, double>()); }
  
   void sortIncrElement()
   { CoinSort_3(elements_, elements_ + nElements_, origIndices_, indices_,
		CoinFirstLess_3<double, int, int>()); }

   void sortDecrElement()
   { CoinSort_3(elements_, elements_ + nElements_, origIndices_, indices_,
		CoinFirstGreater_3<double, int, int>()); }
  

   /** Sort in original order.
       If the vector has been sorted, then this method restores
       to its orignal sort order.
   */
   void sortOriginalOrder();
   //@}

   /**@name Memory usage */
   //@{
   /** Reserve space.
       If one knows the eventual size of the packed vector,
       then it may be more efficient to reserve the space.
   */
   void reserve(int n);
   /** capacity returns the size which could be accomodated without
       having to reallocate storage.
   */
   int capacity() const { return capacity_; }
   //@}

   /**@name Constructors and destructors */
   //@{
   /** Default constructor */
   OsiPackedVector(bool testForDuplicateIndex = true);
   /** Alternate Constructors - set elements to vector of doubles */
   OsiPackedVector(int size, const int * inds, const double * elems,
		   bool testForDuplicateIndex = true);
   /** Alternate Constructors - set elements to same scalar value */
   OsiPackedVector(int size, const int * inds, double element,
		   bool testForDuplicateIndex = true);
   /** Alternate Constructors - construct full storage with indices 0 through
       size-1. */
   OsiPackedVector(int size, const double * elements,
		   bool testForDuplicateIndex = true);
   /** Copy constructor. */
   OsiPackedVector(const OsiPackedVector &);
   /** Copy constructor <em>from a PackedVectorBase</em>. */
   OsiPackedVector(const OsiPackedVectorBase & rhs);
   /** Destructor */
   virtual ~OsiPackedVector ();
   //@}
    
private:
   /**@name Private methods */
   //@{  
   /// Copy internal date
   void gutsOfSetVector(int size,
			const int * inds, const double * elems,
			bool testForDuplicateIndex,
			const char * method);
   ///
   void gutsOfSetConstant(int size,
			  const int * inds, double value,
			  bool testForDuplicateIndex,
			  const char * method);
   //@}

private:
   /**@name Private member data */
   //@{
   /// Vector indices
   int * indices_;
   ///Vector elements
   double * elements_;
   /// Size of indices and elements vectors
   int nElements_;
   /// original unsorted indices
   int * origIndices_;
   /// Amount of memory allocated for indices_, origIndices_, and elements_.
   int capacity_;
   //@}
};

//#############################################################################

/**@name Arithmetic operators on packed vectors.

   <strong>NOTE</strong>: These methods operate on those positions where at
   least one of the arguments has a value listed. At those positions the
   appropriate operation is executed, Otherwise the result of the operation is
   considered 0.<br>
   <strong>NOTE 2</strong>: Because these methods return an object (they can't
   return a reference, though they could return a pointer...) they are
   <em>very</em> inefficient...
 */
//@{
/// Return the sum of two packed vectors
inline OsiPackedVector operator+(const OsiPackedVectorBase& op1,
				 const OsiPackedVectorBase& op2)
{ return op1.binaryOp(op2, std::plus<double>()); }

/// Return the difference of two packed vectors
inline OsiPackedVector operator-(const OsiPackedVectorBase& op1,
				 const OsiPackedVectorBase& op2)
{ return op1.binaryOp(op2, std::minus<double>()); }

/// Return the element-wise product of two packed vectors
inline OsiPackedVector operator*(const OsiPackedVectorBase& op1,
				 const OsiPackedVectorBase& op2)
{ return op1.binaryOp(op2, std::multiplies<double>()); }

/// Return the element-wise ratio of two packed vectors
inline OsiPackedVector operator/(const OsiPackedVectorBase& op1,
				 const OsiPackedVectorBase& op2)
{ return op1.binaryOp(op2, std::divides<double>()); }
//@}


/**@name Arithmetic operators on packed vector and a constant. <br>
   These functions create a packed vector as a result. That packed vector will
   have the same indices as <code>op1</code> and the specified operation is
   done entry-wise with the given value. */
//@{
/// Return the sum of a packed vector and a constant
inline OsiPackedVector operator+(const OsiPackedVectorBase& op1, double value)
{ return op1.binaryOp(std::plus<double>(), value); }

/// Return the difference of a packed vector and a constant
inline OsiPackedVector operator-(const OsiPackedVectorBase& op1, double value)
{ return op1.binaryOp(std::minus<double>(), value); }

/// Return the element-wise product of a packed vector and a constant
inline OsiPackedVector operator*(const OsiPackedVectorBase& op1, double value)
{ return op1.binaryOp(std::multiplies<double>(), value); }

/// Return the element-wise ratio of a packed vector and a constant
inline OsiPackedVector operator/(const OsiPackedVectorBase& op1, double value)
{ return op1.binaryOp(std::divides<double>(), value); }

/// Return the sum of a constant and a packed vector
inline OsiPackedVector operator+(double value, const OsiPackedVectorBase& op1)
{ return op1.binaryOp(value, std::plus<double>()); }

/// Return the difference of a constant and a packed vector
inline OsiPackedVector operator-(double value, const OsiPackedVectorBase& op1)
{ return op1.binaryOp(value, std::minus<double>()); }

/// Return the element-wise product of a constant and a packed vector
inline OsiPackedVector operator*(double value, const OsiPackedVectorBase& op1)
{ return op1.binaryOp(value, std::multiplies<double>()); }

/// Return the element-wise ratio of a a constant and packed vector
inline OsiPackedVector operator/(double value, const OsiPackedVectorBase& op1)
{ return op1.binaryOp(value, std::divides<double>()); }
//@}

//#############################################################################
/** A function that tests the methods in the OsiPackedVector class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
OsiPackedVectorUnitTest();

#endif
