// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiMpsReader_H
#define OsiMpsReader_H

#include "OsiPackedMatrix.hpp"

//#############################################################################

/** MPS Reader Interface

    This can be used to read in mps files without a solver.
    After reading the file this contains all relevant data which
    may be more than the OsiSolverInterface allows for.  Items may be
    deleted to allow for flexibility of data storage.

    This design makes it look very like a dummy solver as the same
    conventions are used.
*/

class OsiMpsReader {
   friend void OsiMpsReaderUnitTest(const std::string & mpsDir);

public:
  /**@name Problem information methods 

   These methods return
   information about the problem referred to by the current object.
   Querying a problem that has no data associated with it result in
   zeros for the number of rows and columns, and NULL pointers from
   the methods that return vectors.

   Const pointers returned from any data-query method are always valid
  */
  //@{
    /// Get number of columns
    int getNumCols() const;

    /// Get number of rows
    int getNumRows() const;

    /// Get number of nonzero elements
    int getNumElements() const;

    /// Get pointer to array[getNumCols()] of column lower bounds
    const double * getColLower() const;

    /// Get pointer to array[getNumCols()] of column upper bounds
    const double * getColUpper() const;

    /** Get pointer to array[getNumRows()] of row constraint senses.
	<ul>
	<li>'L': <= constraint
	<li>'E': =  constraint
	<li>'G': >= constraint
	<li>'R': ranged constraint
	<li>'N': free constraint
	</ul>
    */
    const char * getRowSense() const;

    /** Get pointer to array[getNumRows()] of rows right-hand sides
	<ul>
	  <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
	  <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
	  <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
	  <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
	</ul>
    */
    const double * getRightHandSide() const;

    /** Get pointer to array[getNumRows()] of row ranges.
	<ul>
          <li> if rowsense()[i] == 'R' then
                  rowrange()[i] == rowupper()[i] - rowlower()[i]
          <li> if rowsense()[i] != 'R' then
                  rowrange()[i] is 0.0
        </ul>
    */
    const double * getRowRange() const;

    /// Get pointer to array[getNumRows()] of row lower bounds
    const double * getRowLower() const;

    /// Get pointer to array[getNumRows()] of row upper bounds
    const double * getRowUpper() const;

    /// Get pointer to array[getNumCols()] of objective function coefficients
    const double * getObjCoefficients() const;

#if 0
    /// Get objective function sense (1 for min (default), -1 for max)
    double getObjSense() const;

    /// Get pointer to array[getNumCols()] of primal solution vector
    const double * getColSolution() const;

    /// Get pointer to array[getNumRows()] of dual prices
    const double * getRowPrice() const;

    /// Get objective function value
    double getObjValue() const;

#endif

  /// Sets infinity!
  void setInfinity(double value);
   /// Gets infinity
   double getInfinity() const;
    /// Return true if column is continuous
    bool isContinuous(int colNumber) const;

    /** Return true if column is integer.
        Note: This function returns true if the the column
        is binary or a general integer.
    */
    bool isInteger(int columnNumber) const;


    /// Get pointer to row-wise copy of matrix
    const OsiPackedMatrix * getMatrixByRow() const;

    /// Get pointer to column-wise copy of matrix
    const OsiPackedMatrix * getMatrixByCol() const;
  
  /**@name Methods to input a problem (returns number of errors) 
    -1 if file not opened */
  //@{
    /// Read an mps file from the given filename
    int readMps(const char *filename,
			 const char *extension = "mps");
    /// Read an mps file from previously given filename (or stdin)
    int readMps();
    /// Set file name
    void setFileName(const char * name);
    /// Get file name
    const char * getFileName() const;
    /// Test if current file exists and readable
    const bool fileReadable() const;
    // Could later allow for file pointers and positioning
    /// Sets default upper bound for integer variables
    void setDefaultBound(int value);
    /// gets default upper bound for integer variables
    int getDefaultBound() const;
  //@}


  /**@name Constructors and destructors */
  //@{
    /// Default Constructor
    OsiMpsReader(); 
      
    /// Copy constructor 
    OsiMpsReader (const OsiMpsReader &);
  
    /// Assignment operator 
    OsiMpsReader &    
    operator=(const OsiMpsReader& rhs);
  
    /// Destructor 
    ~OsiMpsReader ();
  //@}

  /**@name Things not in OsiSolverInterface */
  //@{
    /// array saying if each variable integer
    const char * integerColumns() const;
    /// names - returns NULL if out of range
    const char * rowName(int index) const;
    const char * columnName(int index) const;
  /** names - returns -1 if name not found
      at present rowIndex will return numberRows for objective
      and > numberRows for dropped free rows */
    int rowIndex(const char * name) const;
    int columnIndex(const char * name) const;
  /** objective offset - this is RHS entry for objective row */
  double objectiveOffset() const;
  //@}


  /**@name release storage */
  //@{
    /** Release all information which can be re-calculated e.g. rowsense
	also any row copies OR hash tables for names */
    void releaseRedundantInformation();
    /// Release all row information (lower, upper)
    void releaseRowInformation();
    /// Release all column information (lower, upper, objective)
    void releaseColumnInformation();
    /// Release integer information
    void releaseIntegerInformation();
    /// Release row names
    void releaseRowNames();
    /// Release column names
    void releaseColumnNames();
    /// Release matrix information
    void releaseMatrixInformation();
  
  //@}

private:
  
  
  /**@name Private methods */
  //@{
  
    /// The real work of a copy constructor (used by copy and assignment)
    void gutsOfDestructor();
    void gutsOfCopy(const OsiMpsReader &);
  
    /// The real work of a destructor (used by copy and assignment)
    void freeAll();

  /** A quick inlined function to convert from lb/ub stryle constraint
    definition to sense/rhs/range style */
  inline void
  convertBoundToSense(const double lower, const double upper,
				    char& sense, double& right,
				    double& range) const;
  /** A quick inlined function to convert from sense/rhs/range stryle constraint
    definition to lb/ub style */
  inline void
  convertSenseToBound(const char sense, const double right,
				    const double range,
				    double& lower, double& upper) const;
  //@}

  
  // for hashing
  typedef struct {
    int index, next;
  } OsiHashLink;
  /**@name hash stuff */
  //@{
  /// Creates hash list for names (0 rows, 1 columns)
  void startHash ( char **names, const int number , int section );
  /// This one does it when names are already in
  void startHash ( int section ) const;
  /// Deletes hash storage
  void stopHash ( int section );
  /// Finds match using hash,  -1 not found
  int findHash ( const char *name , int section ) const;
  //@}

  /**@name Protected member data */
  //@{
    /**@name Cached information */
    //@{
      /// Pointer to dense vector of row sense indicators
      mutable char    *rowsense_;
  
      /// Pointer to dense vector of row right-hand side values
      mutable double  *rhs_;
  
     /** Pointer to dense vector of slack upper bounds for range 
         constraints (undefined for non-range rows)
     */
     mutable double  *rowrange_;
  
     /// Pointer to row-wise copy of problem matrix coefficients.
     mutable OsiPackedMatrix *matrixByRow_;  
     /// Pointer to column-wise copy of problem matrix coefficients.
     OsiPackedMatrix *matrixByColumn_;  
     double * rowlower_;
     double * rowupper_;
     double * collower_;
     double * colupper_;
     double * objective_;
     char * integerType_;
     char * fileName_;
     /// Number in hash table
     int numberHash_[2];
     /// Hash tables
     mutable OsiHashLink *hash_[2];
     /// Names linked to hash (0 - row names, 1 column names)
     char **names_[2];
     int numberRows_;
     int numberColumns_;
     int numberElements_;
     /// Upper bound when no bounds for integers
     int defaultBound_; 
     double infinity_;
     /// offset for objective function (i.e. rhs of OBJ row)
     double objectiveOffset_;
    //@}
  //@}

};

//#############################################################################
/** A function that tests the methods in the OsiMpsReader class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. Also, if this method is compiled with
    optimization, the compilation takes 10-15 minutes and the machine pages
    (has 256M core memory!)... */
void
OsiMpsReaderUnitTest(const std::string & mpsDir);

#endif
