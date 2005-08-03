//-----------------------------------------------------------------------------
// name:     OSI Interface for GLPK
//-----------------------------------------------------------------------------
// Copyright (C) 2001, Vivian De Smedt, Braden Hunsaker
// Copyright (C) 2003  University of Pittsburgh
//   University of Pittsburgh coding done by Brady Hunsaker
// All Rights Reserved.

#ifndef OsiGlpkSolverInterface_H
#define OsiGlpkSolverInterface_H

#include <string>
extern "C" {
#include "glpk.h"
}
#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

/** GPLK Solver Interface

    Instantiation of OsiGlpkSolverInterface for GPLK
*/

class OsiGlpkSolverInterface : virtual public OsiSolverInterface {
  friend void OsiGlpkSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);
  
public:
  
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial LP relaxation 
  virtual void initialSolve();
  
  /// Resolve an LP relaxation after problem modification
  virtual void resolve();
  
  /// Invoke solver's built-in enumeration algorithm
  virtual void branchAndBound();
  //@}
  
  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the volume algorithm), the parameter is not yet implemented
     for the solver or simply the value of the parameter is out of the range
     the solver accepts. If a parameter setting call returns false check the
     details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.
  */
  //@{
    // Set an integer parameter
    bool setIntParam(OsiIntParam key, int value);
    // Set an double parameter
    bool setDblParam(OsiDblParam key, double value);
    // Set a string parameter
    bool setStrParam(OsiStrParam key, const std::string & value);
    // Get an integer parameter
    bool getIntParam(OsiIntParam key, int& value) const;
    // Get an double parameter
    bool getDblParam(OsiDblParam key, double& value) const;
    // Get a string parameter
    bool getStrParam(OsiStrParam key, std::string& value) const;
  //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there a numerical difficulties?
    virtual bool isAbandoned() const;
    /// Is optimality proven?
    virtual bool isProvenOptimal() const;
    /// Is primal infeasiblity proven?
    virtual bool isProvenPrimalInfeasible() const;
    /// Is dual infeasiblity proven?
    virtual bool isProvenDualInfeasible() const;
    /// Is the given primal objective limit reached?
    virtual bool isPrimalObjectiveLimitReached() const;
    /// Is the given dual objective limit reached?
    virtual bool isDualObjectiveLimitReached() const;
    /// Iteration limit reached?
    virtual bool isIterationLimitReached() const;
  //@}

  //---------------------------------------------------------------------------
  /**@name WarmStart related methods */
  //@{
  /*! \brief Get an empty warm start object
    
    This routine returns an empty CoinWarmStartBasis object. Its purpose is
    to provide a way to give a client a warm start basis object of the
    appropriate type, which can resized and modified as desired.
  */
  inline CoinWarmStart *getEmptyWarmStart () const
  { return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartBasis())) ; }
    /// Get warmstarting information
    virtual CoinWarmStart* getWarmStart() const;
    /** Set warmstarting information. Return true/false depending on whether
	the warmstart information was accepted or not. */
    virtual bool setWarmStart(const CoinWarmStart* warmstart);
  //@}

  //---------------------------------------------------------------------------
  /**@name Hotstart related methods (primarily used in strong branching). <br>
     The user can create a hotstart (a snapshot) of the optimization process
     then reoptimize over and over again always starting from there.<br>
     <strong>NOTE</strong>: between hotstarted optimizations only
     bound changes are allowed. */
  //@{
    /// Create a hotstart point of the optimization process
    virtual void markHotStart();
    /// Optimize starting from the hotstart
    virtual void solveFromHotStart();
    /// Delete the snapshot
    virtual void unmarkHotStart();
  //@}

  //---------------------------------------------------------------------------
  /**@name Problem information methods 
     
     These methods call the solver's query routines to return
     information about the problem referred to by the current object.
     Querying a problem that has no data associated with it result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.
     
     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
  */
  //@{
    /**@name Methods related to querying the input data */
    //@{
      /// Get number of columns
      virtual int getNumCols() const;
  
      /// Get number of rows
      virtual int getNumRows() const;
  
      /// Get number of nonzero elements
      virtual int getNumElements() const;
	  
      /// Get pointer to array[getNumCols()] of column lower bounds
      virtual const double * getColLower() const;
  
      /// Get pointer to array[getNumCols()] of column upper bounds
      virtual const double * getColUpper() const;
  
      /** Get pointer to array[getNumRows()] of row constraint senses.
  	<ul>
  	<li>'L': <= constraint
  	<li>'E': =  constraint
  	<li>'G': >= constraint
  	<li>'R': ranged constraint
  	<li>'N': free constraint
  	</ul>
      */
      virtual const char * getRowSense() const;
  
      /** Get pointer to array[getNumRows()] of rows right-hand sides
  	<ul>
  	  <li> if rowsense()[i] == 'L' then rhs()[i] == rowupper()[i]
  	  <li> if rowsense()[i] == 'G' then rhs()[i] == rowlower()[i]
  	  <li> if rowsense()[i] == 'R' then rhs()[i] == rowupper()[i]
  	  <li> if rowsense()[i] == 'N' then rhs()[i] == 0.0
  	</ul>
      */
      virtual const double * getRightHandSide() const;
  
      /** Get pointer to array[getNumRows()] of row ranges.
  	<ul>
            <li> if rowsense()[i] == 'R' then
                    rowrange()[i] == rowupper()[i] - rowlower()[i]
            <li> if rowsense()[i] != 'R' then
                    rowrange()[i] is 0.0
          </ul>
      */
      virtual const double * getRowRange() const;
  
      /// Get pointer to array[getNumRows()] of row lower bounds
      virtual const double * getRowLower() const;
  
      /// Get pointer to array[getNumRows()] of row upper bounds
      virtual const double * getRowUpper() const;
  
      /// Get pointer to array[getNumCols()] of objective function coefficients
      virtual const double * getObjCoefficients() const;
  
      /// Get objective function sense (1 for min (default), -1 for max)
      virtual double getObjSense() const;

      /// Return true if column is continuous
      virtual bool isContinuous(int colNumber) const;

#if 0
      /// Return true if column is binary
      virtual bool isBinary(int columnNumber) const;
  
      /** Return true if column is integer.
          Note: This function returns true if the the column
          is binary or a general integer.
      */
      virtual bool isInteger(int columnNumber) const;
  
      /// Return true if column is general integer
      virtual bool isIntegerNonBinary(int columnNumber) const;
  
      /// Return true if column is binary and not fixed at either bound
      virtual bool isFreeBinary(int columnNumber) const;
#endif
  
      /// Get pointer to row-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByRow() const;
  
      /// Get pointer to column-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByCol() const;
  
      /// Get solver's value for infinity
      virtual double getInfinity() const;
    //@}
    
    /**@name Methods related to querying the solution */
    //@{
      /// Get pointer to array[getNumCols()] of primal solution vector
      virtual const double * getColSolution() const;
  
      /// Get pointer to array[getNumRows()] of dual prices
      virtual const double * getRowPrice() const;
  
      /// Get a pointer to array[getNumCols()] of reduced costs
      virtual const double * getReducedCost() const;
  
      /** Get pointer to array[getNumRows()] of row activity levels (constraint
  	matrix times the solution vector */
      virtual const double * getRowActivity() const;
  
      /// Get objective function value
      virtual double getObjValue() const;
  
      /** Get how many iterations it took to solve the problem (whatever
	  "iteration" mean to the solver. */
      virtual int getIterationCount() const;
  
      /** Get as many dual rays as the solver can provide. (In case of proven
          primal infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumRows() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getDualRays(int maxNumRays) const;
      /** Get as many primal rays as the solver can provide. (In case of proven
          dual infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumCols() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getPrimalRays(int maxNumRays) const;
  
#if 0
      /** Get vector of indices of solution which are integer variables 
  	presently at fractional values */
      virtual OsiVectorInt getFractionalIndices(const double etol=1.e-05)
	const;
#endif
    //@}
  //@}
  
  //---------------------------------------------------------------------------

  /**@name Problem modifying methods */
  //@{
    //-------------------------------------------------------------------------
    /**@name Changing bounds on variables and constraints */
    //@{
      /** Set an objective function coefficient */
      virtual void setObjCoeff( int elementIndex, double elementValue );

      /** Set a single column lower bound<br>
    	  Use -DBL_MAX for -infinity. */
      virtual void setColLower( int elementIndex, double elementValue );
      
      /** Set a single column upper bound<br>
    	  Use DBL_MAX for infinity. */
      virtual void setColUpper( int elementIndex, double elementValue );
      
      /** Set a single column lower and upper bound<br>
    	  The default implementation just invokes setColLower() and
    	  setColUpper() */
      virtual void setColBounds( int elementIndex,
    				 double lower, double upper );
    
      /** Set the bounds on a number of columns simultaneously<br>
    	  The default implementation just invokes setColLower() and
    	  setColUpper() over and over again.
    	  @param indexFirst,indexLast pointers to the beginning and after the
	         end of the array of the indices of the variables whose
		 <em>either</em> bound changes
    	  @param boundList the new lower/upper bound pairs for the variables
      */
      virtual void setColSetBounds(const int* indexFirst,
				   const int* indexLast,
				   const double* boundList);
      
      /** Set a single row lower bound<br>
    	  Use -DBL_MAX for -infinity. */
      virtual void setRowLower( int elementIndex, double elementValue );
      
      /** Set a single row upper bound<br>
    	  Use DBL_MAX for infinity. */
      virtual void setRowUpper( int elementIndex, double elementValue );
    
      /** Set a single row lower and upper bound<br>
    	  The default implementation just invokes setRowLower() and
    	  setRowUpper() */
      virtual void setRowBounds( int elementIndex,
    				 double lower, double upper );
    
      /** Set the type of a single row<br> */
      virtual void setRowType(int index, char sense, double rightHandSide,
    			      double range);
    
      /** Set the bounds on a number of rows simultaneously<br>
    	  The default implementation just invokes setRowLower() and
    	  setRowUpper() over and over again.
    	  @param indexFirst,indexLast pointers to the beginning and after the
	         end of the array of the indices of the constraints whose
		 <em>either</em> bound changes
    	  @param boundList the new lower/upper bound pairs for the constraints
      */
      virtual void setRowSetBounds(const int* indexFirst,
    				   const int* indexLast,
    				   const double* boundList);
    
      /** Set the type of a number of rows simultaneously<br>
    	  The default implementation just invokes setRowType()
    	  over and over again.
    	  @param indexFirst,indexLast pointers to the beginning and after the
	         end of the array of the indices of the constraints whose
		 <em>any</em> characteristics changes
    	  @param senseList the new senses
    	  @param rhsList   the new right hand sides
    	  @param rangeList the new ranges
      */
      virtual void setRowSetTypes(const int* indexFirst,
				  const int* indexLast,
				  const char* senseList,
				  const double* rhsList,
				  const double* rangeList);
    //@}
    
    //-------------------------------------------------------------------------
    /**@name Integrality related changing methods */
    //@{
      /** Set the index-th variable to be a continuous variable */
      virtual void setContinuous(int index);
      /** Set the index-th variable to be an integer variable */
      virtual void setInteger(int index);
      /** Set the variables listed in indices (which is of length len) to be
	  continuous variables */
      virtual void setContinuous(const int* indices, int len);
      /** Set the variables listed in indices (which is of length len) to be
	  integer variables */
      virtual void setInteger(const int* indices, int len);
    //@}
    
    //-------------------------------------------------------------------------
    /// Set objective function sense (1 for min (default), -1 for max,)
    virtual void setObjSense(double s);
    
    /** Set the primal solution column values
    
    	colsol[numcols()] is an array of values of the problem column
    	variables. These values are copied to memory owned by the
    	solver object or the solver.  They will be returned as the
    	result of colsol() until changed by another call to
    	setColsol() or by a call to any solver routine.  Whether the
    	solver makes use of the solution in any way is
    	solver-dependent. 
    */
    virtual void setColSolution(const double * colsol);
    
    /** Set dual solution vector
    
    	rowprice[numrows()] is an array of values of the problem row
    	dual variables. These values are copied to memory owned by the
    	solver object or the solver.  They will be returned as the
    	result of rowprice() until changed by another call to
    	setRowprice() or by a call to any solver routine.  Whether the
    	solver makes use of the solution in any way is
    	solver-dependent. 
    */
    virtual void setRowPrice(const double * rowprice);
    
    //-------------------------------------------------------------------------
    /**@name Methods to expand a problem.<br>
       Note that if a column is added then by default it will correspond to a
       continuous variable. */
    //@{
      /** */
      virtual void addCol(const CoinPackedVectorBase& vec,
			  const double collb, const double colub,   
			  const double obj);
      /** */
      virtual void addCols(const int numcols,
			   const CoinPackedVectorBase * const * cols,
			   const double* collb, const double* colub,   
			   const double* obj);
      /** */
      virtual void deleteCols(const int num, const int * colIndices);
    
      /** */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const double rowlb, const double rowub);
      /** */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const char rowsen, const double rowrhs,   
    			  const double rowrng);
      /** */
      virtual void addRows(const int numrows,
			   const CoinPackedVectorBase * const * rows,
			   const double* rowlb, const double* rowub);
      /** */
      virtual void addRows(const int numrows,
			   const CoinPackedVectorBase * const * rows,
    			   const char* rowsen, const double* rowrhs,   
    			   const double* rowrng);
      /** */
      virtual void deleteRows(const int num, const int * rowIndices);
    
#if 0
  // ??? implemented in OsiSolverInterface
      //-----------------------------------------------------------------------
      /** Apply a collection of cuts.<br>
    	  Only cuts which have an <code>effectiveness >= effectivenessLb</code>
    	  are applied.
    	  <ul>
    	    <li> ReturnCode.numberIneffective() -- number of cuts which were
                 not applied because they had an
    	         <code>effectiveness < effectivenessLb</code>
    	    <li> ReturnCode.numberInconsistent() -- number of invalid cuts
    	    <li> ReturnCode.numberInconsistentWrtIntegerModel() -- number of
                 cuts that are invalid with respect to this integer model
            <li> ReturnCode.numberInfeasible() -- number of cuts that would
    	         make this integer model infeasible
            <li> ReturnCode.numberApplied() -- number of integer cuts which
    	         were applied to the integer model
            <li> cs.size() == numberIneffective() +
                              numberInconsistent() +
    			      numberInconsistentWrtIntegerModel() +
    			      numberInfeasible() +
    			      nubmerApplied()
          </ul>
      */
      virtual ApplyCutsReturnCode applyCuts(const OsiCuts & cs,
    					    double effectivenessLb = 0.0);
#endif
    //@}
  //@}

  //---------------------------------------------------------------------------

  /**@name Methods to input a problem */
  //@{
    /** Load in an problem by copying the arguments (the constraints on the
        rows are given by lower and upper bounds). If a pointer is 0 then the
        following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub);
			    
    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by lower and upper bounds). For
        default values see the previous method. <br>
	<strong>WARNING</strong>: The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       double*& rowlb, double*& rowub);

    /** Load in an problem by copying the arguments (the constraints on the
	rows are given by sense/rhs/range triplets). If a pointer is 0 then the
	following values are the default:
	<ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
	  <li> <code>obj</code>: all variables have 0 objective coefficient
          <li> <code>rowsen</code>: all rows are >=
          <li> <code>rowrhs</code>: all right hand sides are 0
          <li> <code>rowrng</code>: 0 for the ranged rows
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng);

    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by sense/rhs/range triplets). For
        default values see the previous method. <br>
	<strong>WARNING</strong>: The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       char*& rowsen, double*& rowrhs,
			       double*& rowrng);

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const int* start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub);

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const int* start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng);

    /** Read an mps file from the given filename */
    virtual int readMps(const char *filename,
			 const char *extension = "mps");

    /** Write the problem into an mps file of the given filename.
     If objSense is non zero then -1.0 forces the code to write a
    maximization objective and +1.0 to write a minimization one.
    If 0.0 then solver can do what it wants */
    virtual void writeMps(const char *filename,
			  const char *extension = "mps",
			  double objSense=0.0) const;
  //@}

  //---------------------------------------------------------------------------

  /**@name GLPK specific public interfaces */
  //@{
  enum keepCachedFlag
  {
    /// discard all cached data (default)
    KEEPCACHED_NONE    = 0,
    /// column information: objective values, lower and upper bounds, variable types
    KEEPCACHED_COLUMN  = 1,
    /// row information: right hand sides, ranges and senses, lower and upper bounds for row
    KEEPCACHED_ROW     = 2,
    /// problem matrix: matrix ordered by column and by row
    KEEPCACHED_MATRIX  = 4,
    /// LP solution: primal and dual solution, reduced costs, row activities
    KEEPCACHED_RESULTS = 8,
    /// only discard cached LP solution
    KEEPCACHED_PROBLEM = KEEPCACHED_COLUMN | KEEPCACHED_ROW | KEEPCACHED_MATRIX,
    /// keep all cached data (similar to getMutableLpPtr())
    KEEPCACHED_ALL     = KEEPCACHED_PROBLEM | KEEPCACHED_RESULTS,
    /// free only cached column and LP solution information
    FREECACHED_COLUMN  = KEEPCACHED_PROBLEM & ~KEEPCACHED_COLUMN,
    /// free only cached row and LP solution information
    FREECACHED_ROW     = KEEPCACHED_PROBLEM & ~KEEPCACHED_ROW,
    /// free only cached matrix and LP solution information
    FREECACHED_MATRIX  = KEEPCACHED_PROBLEM & ~KEEPCACHED_MATRIX,
    /// free only cached LP solution information
    FREECACHED_RESULTS = KEEPCACHED_ALL & ~KEEPCACHED_RESULTS
  };

  /// Get pointer to GLPK model
  LPX * getModelPtr();

  //@}

  /**@name Constructors and destructor */
  //@{
  /// Default Constructor
  OsiGlpkSolverInterface(); 
  
  /// Clone
  virtual OsiSolverInterface * clone(bool copyData = true) const;
  
  /// Copy constructor 
  OsiGlpkSolverInterface( const OsiGlpkSolverInterface& );
  
  /// Assignment operator 
  OsiGlpkSolverInterface& operator=( const OsiGlpkSolverInterface& rhs );
  
  /// Destructor 
  virtual ~OsiGlpkSolverInterface();

  /// Resets as if default constructor
  virtual void reset();
  //@}
  
protected:
  
  /**@name Protected methods */
  //@{
  /// Apply a row cut. Return true if cut was applied.
  virtual void applyRowCut( const OsiRowCut & rc );
  
  /** Apply a column cut (bound adjustment). 
      Return true if cut was applied.
  */
  virtual void applyColCut( const OsiColCut & cc );

  /// Pointer to the model
  LPX * getMutableModelPtr() const;

  //@}
 
private: 
  /**@name Private methods */
  //@{
    
  /// The real work of a copy constructor (used by copy and assignment)
  void gutsOfCopy( const OsiGlpkSolverInterface & source );
  
  /// The real work of the constructor
  void gutsOfConstructor();
  
  /// The real work of the destructor
  void gutsOfDestructor();

  /// free cached column rim vectors
  void freeCachedColRim();

  /// free cached row rim vectors
  void freeCachedRowRim();

  /// free cached result vectors
  void freeCachedResults();
  
  /// free cached matrices
  void freeCachedMatrix();

  /// free all cached data (except specified entries, see getLpPtr())
  void freeCachedData( int keepCached = KEEPCACHED_NONE );

  /// free all allocated memory
  void freeAllMemory();

  /// Just for testing purposes
  void printBounds(); 

  /// Fill cached collumn bounds
  void fillColBounds() const;
  //@}
  
  
  /**@name Private member data */
  //@{
  /// GPLK model represented by this class instance
  mutable LPX* lp_;

  // Remember whether simplex or b&b was most recently done
  // 0 = simplex;  1 = b&b
  int bbWasLast_; 

  // Int parameters.
  int maxIteration_;
  int hotStartMaxIteration_;

  // Double parameters.
  double dualObjectiveLimit_;
  double primalObjectiveLimit_;
  double dualTolerance_;
  double primalTolerance_;

  /// Hotstart information

  int hotStartCStatSize_;
  int *hotStartCStat_;
  double *hotStartCVal_;
  double *hotStartCDualVal_;

  int hotStartRStatSize_;
  int *hotStartRStat_;
  double *hotStartRVal_;
  double *hotStartRDualVal_;

  bool isIterationLimitReached_;
  bool isAbandoned_;
  bool isObjLowerLimitReached_;
  bool isObjUpperLimitReached_;
  bool isPrimInfeasible_;
  bool isDualInfeasible_;

  /**@name Cached information derived from the GLPK model */
  //@{

  /// Number of iterations
  mutable int iter_used_;

  /// Pointer to objective vector
  mutable double  *obj_;
  
  /// Pointer to dense vector of variable lower bounds
  mutable double  *collower_;
  
  /// Pointer to dense vector of variable lower bounds
  mutable double  *colupper_;
  
  /// Pointer to dense vector of variable types (continous, binary, integer)
  mutable char    *ctype_;
  
  /// Pointer to dense vector of row sense indicators
  mutable char    *rowsense_;
  
  /// Pointer to dense vector of row right-hand side values
  mutable double  *rhs_;
  
  /// Pointer to dense vector of slack upper bounds for range constraints (undefined for non-range rows)
  mutable double  *rowrange_;
  
  /// Pointer to dense vector of row lower bounds
  mutable double  *rowlower_;
  
  /// Pointer to dense vector of row upper bounds
  mutable double  *rowupper_;
  
  /// Pointer to primal solution vector
  mutable double  *colsol_;
  
  /// Pointer to dual solution vector
  mutable double  *rowsol_;

  /// Pointer to reduced cost vector
  mutable double  *redcost_;

  /// Pointer to row activity (slack) vector
  mutable double  *rowact_;

  /// Pointer to row-wise copy of problem matrix coefficients.
  mutable CoinPackedMatrix *matrixByRow_;  
  
  /// Pointer to row-wise copy of problem matrix coefficients.
  mutable CoinPackedMatrix *matrixByCol_;  
  //@}
  //@}
};

//#############################################################################
/** A function that tests the methods in the OsiGlpkSolverInterface class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void OsiGlpkSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);

#endif // OsiGlpkSolverInterface_H
