// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiPackedMatrix_H
#define OsiPackedMatrix_H

#include "CoinError.hpp"
#include "OsiPackedVectorBase.hpp"
#include "OsiShallowPackedVector.hpp"

/** Sparse Matrix Base Class

  This class is used for storing a matrix by rows or columns.

  The sparse represention can be completely compact or it can 
  have "extra" space. The extra space can be added at the end 
  of rows and/or columns.  Incorporating extra space into the 
  sparse matrix representation can improve performance in 
  cases where new data needs to be inserted into the packed 
  matrix.

  For example if the matrix:
  <pre>
     3  1  0   -2   -1  0  0   -1                 
     0  2  1.1  0    0  0  0    0                       
     0  0  1    0    0  1  0    0         
     0  0  0    2.8  0  0 -1.2  0   
   5.6  0  0    0    1  0  0    1.9

  was stored by rows (with no extra space) in 
  OsiPackedMatrix r then: 
    r.elements() returns a vector containing: 
      3 1 -2 -1 -1 2 1.1 1 1 2.8 -1.2 5.6 1 1.9 
    r.indices() returns a vector containing: 
      0 1  3  4  7 1 2   2 5 3    6   0   4 7 
    r.vectorStarts() returns a vector containing: 
      0 5 7 9 11 14 
    r.size() returns 14. 
    r.majorDim() returns 5. 
    r.vectorSize(0) returns 5. 
    r.vectorSize(1) returns 2. 
    r.vectorSize(2) returns 2. 
    r.vectorSize(3) returns 2. 
    r.vectorSize(4) returns 3. 
 
  If stored by columns (with no extra space) then: 
    c.elements() returns a vector containing: 
      3 5.6 1 2 1.1 1 -2 2.8 -1 1 1 -1.2 -1 1.9 
    c.indices() returns a vector containing: 
      0  4  0 1 1   2  0 3    0 4 2  3    0 4 
    c.vectorStarts() returns a vector containing: 
      0 2 4 6 8 10 11 12 14 
    c.size() returns 14. 
    c.majorDim() returns 8. 
  </pre>
*/
class OsiPackedMatrix  {
   friend void OsiPackedMatrixUnitTest();

public:

  //---------------------------------------------------------------------------
  /**@name Query members */
  //@{
    /** Return the current setting of the extra gap. */
    double getExtraGap() const { return extraGap_; }
    /** Return the current setting of the extra major. */
    double getExtraMajor() const { return extraMajor_; }

    /** Reserve sufficient space for appending majorordered vectors. */
    void reserve(const int newMaxMajorDim, const int newMaxSize);
    /** Clear the data, but do not free any arrays */
    void clear();

    /** Whether the packed matrix is column major ordered or not. */
    bool isColOrdered() const { return colOrdered_; }
    /** Number of entries in the packed matrix. */
    int getNumElements() const { return size_; }
    /** Number of columns. */
    int getNumCols() const { return colOrdered_ ? majorDim_ : minorDim_; }
    /** Number of rows. */
    int getNumRows() const { return colOrdered_ ? minorDim_ : majorDim_; }

    /** A vector containing the elements in the packed matrix. Note that there
	might be gaps in this list, entries that do not belong to any
	major-dimension vector. To get the actual elements one should look at
	this vector together with vectorStarts and vectorLengths. */
    const double * getElements() const { return element_; }
    /** A vector containing the minor indices of the elements in the packed
        matrix. Note that there might be gaps in this list, entries that do not
        belong to any major-dimension vector. To get the actual elements one
        should look at this vector together with vectorStarts and
        vectorLengths. */
    const int * getIndices() const { return index_; }

    /** The size of the <code>vectorStarts</code> array */
    int getSizeVectorStarts()const { return majorDim_ > 0 ? majorDim_+1 : 0;}
    /** The size of the <code>vectorLengths</code> array */
    int getSizeVectorLengths() const { return majorDim_; }
    /** The positions where the major-dimension vectors start in elements and
        indices. */
    const int * getVectorStarts() const { return start_; }
    /** The lengths of the major-dimension vectors. */
    const int * getVectorLengths() const { return length_; }


    /** The position of the first element in the i'th major-dimension vector.
     */
    int getVectorFirst(const int i) const throw(CoinError) {
      if (i < 0 || i >= majorDim_)
	throw CoinError("bad index", "vectorFirst", "OsiPackedMatrix");
      return start_[i];
    }
    /** The position of the last element (well, one entry <em>past</em> the
        last) in the i'th major-dimension vector. */
    int getVectorLast(const int i) const throw(CoinError) {
      if (i < 0 || i >= majorDim_)
	throw CoinError("bad index", "vectorLast", "OsiPackedMatrix");
      return start_[i] + length_[i];
    }
    /** The length of i'th vector. */
    int getVectorSize(const int i) const throw(CoinError) {
      if (i < 0 || i >= majorDim_)
	throw CoinError("bad index", "vectorSize", "OsiPackedMatrix");
      return length_[i];
    }
  
    /** Return the i'th vector in matrix. */
    const OsiShallowPackedVector getVector(int i) const throw(CoinError) {
      if (i < 0 || i >= majorDim_)
	throw CoinError("bad index", "vector", "OsiPackedMatrix");
      return OsiShallowPackedVector(length_[i],
  				    index_ + start_[i],
  				    element_ + start_[i],
  				    false);
    }
  //@}

  //---------------------------------------------------------------------------
  /**@name Modifying members. */
  //@{
    /** Set the dimansions of the matrix. In effect, append new empty
	columns/rows to the matrix. A negative number for either dimension
	means that that dimension doesn't change. Otherwise the new dimensions
	MUST be at least as large as the current ones otherwise an exception
	is thrown. */
    void setDimensions(int numrows, int numcols) throw(CoinError);
   
    /** Set the extra gap to be allocated to the specified value. */
    void setExtraGap(const double newGap) throw(CoinError);
    /** Set the extra major to be allocated to the specified value. */
    void setExtraMajor(const double newMajor) throw(CoinError);

    /** Append a column to the end of the matrix. When libosi is compiled with
	a OSI_DEBUG defined then this method throws an exception if the new
	column contains an index that's larger than the number of rows (-1).
	Otherwise the method assumes that every index fits into the matrix. */
    void appendCol(const OsiPackedVectorBase& vec) throw(CoinError);
    /** Append a column to the end of the matrix. When libosi is compiled with
	a OSI_DEBUG defined then this method throws an exception if the new
	column contains an index that's larger than the number of rows (-1).
	Otherwise the method assumes that every index fits into the matrix. */
    void appendCol(const int vecsize,
  		  const int *vecind, const double *vecelem) throw(CoinError);
    /** Append a set of columns to the end of the matrix. When libosi is
	compiled with a OSI_DEBUG defined then this method throws an exception
	if any of the new columns contain an index that's larger than the
	number of rows (-1). Otherwise the method assumes that every index
	fits into the matrix. */
    void appendCols(const int numcols,
		    const OsiPackedVectorBase * const * cols) throw(CoinError);

  /** Append a row to the end of the matrix. When libosi is compiled with
	a OSI_DEBUG defined then this method throws an exception if the new
	row contains an index that's larger than the number of columns (-1).
	Otherwise the method assumes that every index fits into the matrix. */
    void appendRow(const OsiPackedVectorBase& vec) throw(CoinError);
    /** Append a row to the end of the matrix. When libosi is compiled with
	a OSI_DEBUG defined then this method throws an exception if the new
	row contains an index that's larger than the number of columns (-1).
	Otherwise the method assumes that every index fits into the matrix. */
    void appendRow(const int vecsize,
  		  const int *vecind, const double *vecelem) throw(CoinError);
    /** Append a set of rows to the end of the matrix. When libosi is
	compiled with a OSI_DEBUG defined then this method throws an exception
	if any of the new rows contain an index that's larger than the
	number of columns (-1). Otherwise the method assumes that every index
	fits into the matrix. */
    void appendRows(const int numrows,
		    const OsiPackedVectorBase * const * rows) throw(CoinError);
  
    /** Append the argument to the "right" of the current matrix. Imagine this
        as adding new columns (don't worry about how the matrices are ordered,
        that is taken care of). An exception is thrown if the number of rows
        is different in the matrices. */
    void rightAppendPackedMatrix(const OsiPackedMatrix& matrix)
      throw(CoinError);
    /** Append the argument to the "bottom" of the current matrix. Imagine this
        as adding new rows (don't worry about how the matrices are ordered,
        that is taken care of). An exception is thrown if the number of columns
        is different in the matrices. */
    void bottomAppendPackedMatrix(const OsiPackedMatrix& matrix)
      throw(CoinError);
  
    /** Delete the columns whose indices are listed in <code>indDel</code>. */
    void deleteCols(const int numDel, const int * indDel);
    /** Delete the rows whose indices are listed in <code>indDel</code>. */
    void deleteRows(const int numDel, const int * indDel);

    /** Replace the elements of a vector.  The indices remain the same.
	At most the number specified will be replaced.
        The index is between 0 and major dimension of matrix */
    void replaceVector(const int index,
		       const int numReplace, const double * newElements);
  //@}

  //---------------------------------------------------------------------------
  /**@name Methods that reorganize the whole matrix */
  //@{
    /** Remove the gaps from the matrix if there were any */
    void removeGaps();
 
    /** Extract a submatrix from matrix. Those major-dimension vectors of
	the matrix comprise the submatrix whose indices are given in the
	arguments. */
    void submatrixOf(const OsiPackedMatrix& matrix,
		     const int numMajor, const int * indMajor)
       throw(CoinError);
#if 0
    /** Extract a submatrix from matrix. Those major/minor-dimension vectors of
	the matrix comprise the submatrix whose indices are given in the
	arguments. */
    void submatrixOf(const OsiPackedMatrix& matrix,
		     const int numMajor, const int * indMajor,
		     const int numMinor, const int * indMinor) throw(CoinError);
#endif

    /** Copy method. This method makes an exact replica of the argument,
        including the extra space parameters. */
    void copyOf(const OsiPackedMatrix& rhs);
    /** Copy the arguments to the matrix. If <code>len<code> is a NULL pointer
        then the matrix is assumed to have no gaps in it and <code>len</code>
        will be created accordingly. */
    void copyOf(const bool colordered,
 	       const int minor, const int major, const int numels,
 	       const double * elem, const int * ind,
 	       const int * start, const int * len,
 	       const double extraMajor=0.0, const double extraGap=0.0);
    /** Reverse copy method. This method makes an exact replica of the
        argument, but the major ordering reversed. The extra space parameters
        are copied and reversed, too. */
    void reverseOrderedCopyOf(const OsiPackedMatrix& rhs);
    /** Assign the arguments to the matrix. If <code>len<code> is a NULL pointer
        then the matrix is assumed to have no gaps in it and <code>len</code>
        will be created accordingly. <br>
        <strong>NOTE 1</strong>: After this method returns the pointers
        passed to the method will be NULL pointers! <br>
        <strong>NOTE 2</strong>: When the matrix is eventually destructed the
        arrays will be deleted by <code>delete[]</code>. Hence one should use
        this method ONLY if all array swere allocated by <code>new[]</code>! */
    void assignMatrix(const bool colordered,
 		     const int minor, const int major, const int numels,
 		     double *& elem, int *& ind,
 		     int *& start, int *& len,
 		     const int maxmajor = -1, const int maxsize = -1);
 
 
 
    /** Assignment operator. This copies out the data, but uses the current
        matrix's extra space parameters. */
    OsiPackedMatrix & operator=(const OsiPackedMatrix& rhs);
 
    /** Reverse the ordering of the packed matrix. */
    void reverseOrdering();
    /** Transpose the matrix. <br>
        NOTE: All this routine does is to flip the ordering! Of course, then
        the matrix describes the transposed matrix. To get the matrix
        physically transposed (e.g., for a column ordered matrix to get the
        transpose in column ordered format) one has to invoke this method AND
        <code>reverseOrdering()</code>. */
    void transpose();
 
    /** Swap the content of the two packed matrix. */
    void swap(OsiPackedMatrix& matrix);
   
  //@}

  //---------------------------------------------------------------------------
  /**@name Matrix times vector methods */
  //@{
    /** Return <code>A * x</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numColumns()</code>
        @precond <code>y<code> must be of size <code>numRows()</code> */
    void times(const double * x, double * y) const;
    /** Return <code>A * x</code> in <code>y</code>. Same as the previous
        method, just <code>x</code> is given in the form of a packed vector. */
    void times(const OsiPackedVectorBase& x, double * y) const;
    /** Return <code>x * A</code> in <code>y</code>.
        @precond <code>x<code> must be of size <code>numRows()</code>
        @precond <code>y<code> must be of size <code>numColumns()</code> */
    void transposeTimes(const double * x, double * y) const;
    /** Return <code>x * A</code> in <code>y</code>. Same as the previous
        method, just <code>x</code> is given in the form of a packed vector. */
    void transposeTimes(const OsiPackedVectorBase& x, double * y) const;
  //@}

  //---------------------------------------------------------------------------
  /**@name Helper functions used internally, but maybe useful externally.

     These methods do not worry about testing whether the packed matrix is
     row or column major ordered; they operate under the assumption that the
     correct version is invoked. In fact, a number of other methods simply
     just call one of these after testing the ordering of the matrix. */
  //@{

    //-------------------------------------------------------------------------
    /**@name Queries */
    //@{
      /** Count the number of entries in every minor-dimension vector and
	  return an array containing these lengths. The returned array is
	  allocated with <code>new int[]</code>, free it with
	  <code>delete[]</code>. */
      int * countOrthoLength() const;
      /** Major dimension. For row ordered matrix this would be the number of
          rows. */
      int getMajorDim() const { return majorDim_; }
      /** Minor dimension. For row ordered matrix this would be the number of
	  columns. */
      int getMinorDim() const { return minorDim_; }
      /** Current maximum for major dimension. For row ordered matrix this many
          rows can be added without reallocating the vector related to the
	  major dimension (<code>start_</code> and <code>length_</code>). */
      int getMaxMajorDim() const { return maxMajorDim_; }

      /** Dump the matrix on stdout. When in dire straits this method can
	  help. */
      void dumpMatrix(const char* fname = NULL) const;
    //@}

    //-------------------------------------------------------------------------
    /**@name Append vectors. <br>
       When libosi is compiled with a OSI_DEBUG defined then these methods
       throws an exception if the major (minor) vector contains an index
       that's larger than the minor (major) dimension (-1). Otherwise the
       methods assume that every index fits into the matrix. */
    //@{
      /** Append a major-dimension vector to the end of the matrix. */
      void appendMajorVector(const OsiPackedVectorBase& vec) throw(CoinError);
      /** Append a major-dimension vector to the end of the matrix. */
      void appendMajorVector(const int vecsize, const int *vecind,
			     const double *vecelem) throw(CoinError);
      /** Append several major-dimensonvectors to the end of the matrix */
      void appendMajorVectors(const int numvecs,
			      const OsiPackedVectorBase * const * vecs)
	throw(CoinError);

      /** Append a minor-dimension vector to the end of the matrix. */
      void appendMinorVector(const OsiPackedVectorBase& vec) throw(CoinError);
      /** Append a minor-dimension vector to the end of the matrix. */
      void appendMinorVector(const int vecsize, const int *vecind,
			     const double *vecelem) throw(CoinError);
      /** Append several minor-dimensonvectors to the end of the matrix */
      void appendMinorVectors(const int numvecs,
			      const OsiPackedVectorBase * const * vecs)
	throw(CoinError);
    //@}

    //-------------------------------------------------------------------------
    /**@name Append matrices.

       We'll document these methods assuming that the current matrix is
       column major ordered (Hence in the <code>...SameOrdered()</code>
       methods the argument is column ordered, in the
       <code>OrthoOrdered()</code> methods the argument is row ordered.)
    */
    //@{
      /** Append the columns of the argument to the right end of this matrix.
	  @precondition <code>minorDim_ == matrix.minorDim_</code> <br>
	  This method throws an exception if the minor dimensions are not the
	  same. */
      void majorAppendSameOrdered(const OsiPackedMatrix& matrix)
	 throw(CoinError);
      /** Append the columns of the argument to the bottom end of this matrix.
	  @precondition <code>majorDim_ == matrix.majorDim_</code> <br>
	  This method throws an exception if the major dimensions are not the
	  same. */
      void minorAppendSameOrdered(const OsiPackedMatrix& matrix)
	 throw(CoinError);
      /** Append the rows of the argument to the right end of this matrix.
	  @precondition <code>minorDim_ == matrix.majorDim_</code> <br>
	  This method throws an exception if the minor dimension of the
	  current matrix is not the same as the major dimension of the
	  argument matrix. */
      void majorAppendOrthoOrdered(const OsiPackedMatrix& matrix)
	 throw(CoinError);
      /** Append the rows of the argument to the bottom end of this matrix.
	  @precondition <code>majorDim_ == matrix.minorDim_</code> <br>
	  This method throws an exception if the major dimension of the
	  current matrix is not the same as the minor dimension of the
	  argument matrix. */
      void minorAppendOrthoOrdered(const OsiPackedMatrix& matrix)
	 throw(CoinError);
      //@}

      //-----------------------------------------------------------------------
      /**@name Delete vectors */
      //@{
      /** Delete the major-dimension vectors whose indices are listed in
	  <code>indDel</code>. */
      void deleteMajorVectors(const int numDel,
			      const int * indDel) throw(CoinError);
      /** Delete the minor-dimension vectors whose indices are listed in
	  <code>indDel</code>. */
      void deleteMinorVectors(const int numDel,
			      const int * indDel) throw(CoinError);
      //@}

      //-----------------------------------------------------------------------
      /**@name Various dot products. */
      //@{
      /** Return <code>A * x</code> (multiplied from the "right" direction) in
	  <code>y</code>.
	  @precond <code>x<code> must be of size <code>majorDim()</code>
	  @precond <code>y<code> must be of size <code>minorDim()</code> */
      void timesMajor(const double * x, double * y) const;
      /** Return <code>A * x</code> (multiplied from the "right" direction) in
	  <code>y</code>. Same as the previous method, just <code>x</code> is
	  given in the form of a packed vector. */
      void timesMajor(const OsiPackedVectorBase& x, double * y) const;
      /** Return <code>A * x</code> (multiplied from the "right" direction) in
	  <code>y</code>.
	  @precond <code>x<code> must be of size <code>minorDim()</code>
	  @precond <code>y<code> must be of size <code>majorDim()</code> */
      void timesMinor(const double * x, double * y) const;
      /** Return <code>A * x</code> (multiplied from the "right" direction) in
	  <code>y</code>. Same as the previous method, just <code>x</code> is
	  given in the form of a packed vector. */
      void timesMinor(const OsiPackedVectorBase& x, double * y) const;
      //@}
   //@}

   //--------------------------------------------------------------------------
   /**@name Logical Operations. */
   //@{
   /** Equivalence.
       Two matrices are equivalent if they are both by rows or both by columns,
       they have the same dimensions, and each vector is equivalent. 
       In this method the FloatEqual function operator can be specified. 
   */
   template <class FloatEqual> bool 
   isEquivalent(const OsiPackedMatrix& rhs, const FloatEqual& eq) const
   {
      // Both must be column order or both row ordered and must be of same size
      if ((isColOrdered() ^ rhs.isColOrdered()) ||
	  (getNumCols() != rhs.getNumCols()) ||
	  (getNumRows() != rhs.getNumRows()) ||
	  (getNumElements() != rhs.getNumElements()))
	 return false;
     
      for (int i=getMajorDim()-1; i >= 0; --i) {
        OsiShallowPackedVector pv = getVector(i);
        OsiShallowPackedVector rhsPv = rhs.getVector(i);
        if ( !pv.isEquivalent(rhsPv,eq) )
          return false;
      }
      return true;
   }
   
   ///The default equivalence test is that the entries are relatively equal.
   bool isEquivalent(const OsiPackedMatrix& rhs) const
   {
      return isEquivalent(rhs,  OsiRelFltEq());
   }
  bool isEquivalent2(const OsiPackedMatrix& rhs) const;
   //@}

   //--------------------------------------------------------------------------
   /**@name Constructors, destructors and major modifying methods*/
   //@{
   /// Default Constructor creates an empty column ordered packed matrix
   OsiPackedMatrix();

   /// A constructor where the ordering and the gaps are specified
   OsiPackedMatrix(const bool colordered,
		   const double extraMajor, const double extraGap);

   OsiPackedMatrix(const bool colordered,
		   const int minor, const int major, const int numels,
		   const double * elem, const int * ind,
		   const int * start, const int * len,
		   const double extraMajor, const double extraGap);

   OsiPackedMatrix(const bool colordered,
		   const int minor, const int major, const int numels,
		   const double * elem, const int * ind,
		   const int * start, const int * len);

   /** Create packed matrix from ordered triples.
       If colordered is true then the created matrix will be column ordered.
       Duplicate matrix elements are allowed. The created matrix will have 
       the sum of the duplicates. <br>
       For example if: <br>
         rowIndices[0]=2; colIndices[0]=5; elements[0]=2.0 <br>
         rowIndices[1]=2; colIndices[1]=5; elements[1]=0.5 <br>
       then the created matrix will contain a value of 2.5 in row 2 and column 5.<br>
       The matrix is created without gaps.
   */
   OsiPackedMatrix(const bool colordered,
     const int * rowIndices, 
     const int * colIndices, 
     const double * elements, 
     int numels ); 

   /// Copy constructor 
   OsiPackedMatrix(const OsiPackedMatrix& m);

   /// Destructor 
   virtual ~OsiPackedMatrix();    
   //@}

   //--------------------------------------------------------------------------
protected:
   void gutsOfDestructor();
   void gutsOfCopyOf(const bool colordered,
		     const int minor, const int major, const int numels,
		     const double * elem, const int * ind,
		     const int * start, const int * len,
		     const double extraMajor=0.0, const double extraGap=0.0);
   void gutsOfOpEqual(const bool colordered,
		      const int minor, const int major, const int numels,
		      const double * elem, const int * ind,
		      const int * start, const int * len);
   void resizeForAddingMajorVectors(const int numVec, const int * lengthVec);
   void resizeForAddingMinorVectors(const int * addedEntries);
private:
   inline int getLastStart() const {
      return majorDim_ == 0 ? 0 : start_[majorDim_];
   }

   //--------------------------------------------------------------------------
protected:
   /**@name Data members
      The data members are protected to allow access for derived classes. */
   //@{
   /** A flag indicating whether the matrix is column or row major ordered. */
   bool     colOrdered_;
   /** This much times more space should be allocated for each major-dimension
       vector (with respect to the number of entries in the vector) when the
       matrix is resized. The purpose of these gaps is to allow fast insertion
       of new minor-dimension vectors. */
   double   extraGap_;
   /** his much times more space should be allocated for major-dimension
       vectors when the matrix is resized. The purpose of these gaps is to
       allow fast addition of new major-dimension vectors. */
   double   extraMajor_;

   /** List of nonzero element values. The entries in the gaps between
       major-dimension vectors are undefined. */
   double  *element_;
   /** List of nonzero element minor-dimension indices. The entries in the gaps
       between major-dimension vectors are undefined. */
   int     *index_;
   /** Starting positions of major-dimension vectors. */
   int     *start_;
   /** Lengths of major-dimension vectors. */
   int     *length_;

   /// number of vectors in matrix
   int majorDim_;
   /// size of other dimension
   int minorDim_;
   /// the number of nonzero entries
   int size_;

   /// max space allocated for major-dimension
   int maxMajorDim_;
   /// max space allocated for entries
   int maxSize_;
   //@}
};

//#############################################################################
/** A function that tests the methods in the OsiPackedMatrix class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
OsiPackedMatrixUnitTest();

#endif
