// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiClpSolverInterface_H
#define OsiClpSolverInterface_H

#include <string>
#include <cfloat>
#include <map>

#include "ClpSimplex.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"

class OsiRowCut;

static const double OsiClpInfinity = DBL_MAX;

//#############################################################################

/** Clp Solver Interface

    Instantiation of OsiClpSolverInterface for the Model Algorithm.

    It has a trivial branch and bound code for completeness.  It would be
    an interesting project to expand it and move it to base class.

*/

class OsiClpSolverInterface :
   public OsiSolverInterface {
   friend void OsiClpSolverInterfaceUnitTest(const std::string & mpsDir);

public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve();

    /// Resolve an LP relaxation after problem modification
    virtual void resolve();

    /// Invoke solver's built-in enumeration algorithm
    virtual void branchAndBound() {
      throw CoinError("Sorry, Clp doesn't implement B&B",
		     "branchAndBound", "OsiClpSolverInterface");
    };
  //@}

  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the clp algorithm), the parameter is not yet implemented
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
      virtual int getNumCols() const {
        return modelPtr_->numberColumns(); }
  
      /// Get number of rows
      virtual int getNumRows() const {
        return modelPtr_->numberRows(); }
  
      /// Get number of nonzero elements
      virtual int getNumElements() const {
        return modelPtr_->matrix()->getNumElements(); }

      /// Get pointer to array[getNumCols()] of column lower bounds
      virtual const double * getColLower() const { return modelPtr_->columnLower(); }
  
      /// Get pointer to array[getNumCols()] of column upper bounds
      virtual const double * getColUpper() const { return modelPtr_->columnUpper(); }
  
      /** Get pointer to array[getNumRows()] of row constraint senses.
  	<ul>
  	<li>'L' <= constraint
  	<li>'E' =  constraint
  	<li>'G' >= constraint
  	<li>'R' ranged constraint
  	<li>'N' free constraint
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
      virtual const double * getRightHandSide() const ;
  
      /** Get pointer to array[getNumRows()] of row ranges.
  	<ul>
            <li> if rowsense()[i] == 'R' then
                    rowrange()[i] == rowupper()[i] - rowlower()[i]
            <li> if rowsense()[i] != 'R' then
                    rowrange()[i] is undefined
          </ul>
      */
      virtual const double * getRowRange() const ;
  
      /// Get pointer to array[getNumRows()] of row lower bounds
      virtual const double * getRowLower() const { return modelPtr_->rowLower(); }
  
      /// Get pointer to array[getNumRows()] of row upper bounds
      virtual const double * getRowUpper() const { return modelPtr_->rowUpper(); }
  
      /// Get pointer to array[getNumCols()] of objective function coefficients
      virtual const double * getObjCoefficients() const 
       { return modelPtr_->objective(); }
  
      /// Get objective function sense (1 for min (default), -1 for max)
      virtual double getObjSense() const 
      { return modelPtr_->optimizationDirection(); }
  
       /// Return true if column is continuous
      virtual bool isContinuous(int colNumber) const;
  
  
      /// Get pointer to row-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByRow() const;
  
      /// Get pointer to column-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByCol() const;
  
      /// Get solver's value for infinity
      virtual double getInfinity() const { return OsiClpInfinity; }
    //@}
    
    /**@name Methods related to querying the solution */
    //@{
      /// Get pointer to array[getNumCols()] of primal solution vector
      virtual const double * getColSolution() const 
      { return modelPtr_->primalColumnSolution(); }
  
      /// Get pointer to array[getNumRows()] of dual prices
      virtual const double * getRowPrice() const
      { return modelPtr_->dualRowSolution(); }
  
      /// Get a pointer to array[getNumCols()] of reduced costs
      virtual const double * getReducedCost() const 
      { return modelPtr_->dualColumnSolution(); }
  
      /** Get pointer to array[getNumRows()] of row activity levels (constraint
  	matrix times the solution vector */
      virtual const double * getRowActivity() const 
      { return modelPtr_->primalRowSolution(); }
  
      /// Get objective function value
      virtual double getObjValue() const { return modelPtr_->objectiveValue(); }
  
      /** Get how many iterations it took to solve the problem (whatever
	  "iteration" mean to the solver. */
      virtual int getIterationCount() const 
       { return modelPtr_->numberIterations(); }
  
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
  
    //@}
  //@}

  //---------------------------------------------------------------------------

  /**@name Problem modifying methods */
  //@{
    //-------------------------------------------------------------------------
    /**@name Changing bounds on variables and constraints */
    //@{
      /** Set an objective function coefficient */
      virtual void setObjCoeff( int elementIndex, double elementValue ) {
        modelPtr_->objective()[elementIndex] = elementValue;
      }

      /** Set a single column lower bound<br>
    	  Use -DBL_MAX for -infinity. */
      virtual void setColLower( int elementIndex, double elementValue ) {
	modelPtr_->columnLower()[elementIndex] = elementValue;
      }
      
      /** Set a single column upper bound<br>
    	  Use DBL_MAX for infinity. */
      virtual void setColUpper( int elementIndex, double elementValue ) {
	modelPtr_->columnUpper()[elementIndex] = elementValue;
      }

      /** Set a single column lower and upper bound */
      virtual void setColBounds( int elementIndex,
    				 double lower, double upper ) {
	modelPtr_->columnLower()[elementIndex] = lower;
	modelPtr_->columnUpper()[elementIndex] = upper;
      }

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
      virtual void setRowLower( int elementIndex, double elementValue ) {
	modelPtr_->rowLower()[elementIndex] = elementValue;
      }
      
      /** Set a single row upper bound<br>
    	  Use DBL_MAX for infinity. */
      virtual void setRowUpper( int elementIndex, double elementValue ) {
	modelPtr_->rowUpper()[elementIndex] = elementValue;
      }
    
      /** Set a single row lower and upper bound */
      virtual void setRowBounds( int elementIndex,
    				 double lower, double upper ) {
	modelPtr_->rowLower()[elementIndex] = lower;
	modelPtr_->rowUpper()[elementIndex] = upper;
      }
    
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
    virtual void setObjSense(double s ) 
     { modelPtr_->setOptimizationDirection( s < 0 ? -1 : 1); }
    
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
    
      //-----------------------------------------------------------------------
      /** Apply a collection of row cuts which are all effective.
	  applyCuts seems to do one at a time which seems inefficient.
      */
      virtual void applyRowCuts(int numberCuts, const OsiRowCut * cuts);
    //@}
  //@}

  //---------------------------------------------------------------------------

public:
   
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
        default values see the previous method.  <br>
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
#if 0    
    /** Read an mps file from the given filename */
    virtual int readMps(const char *filename,
			 const char *extension = "mps");
#endif

    /** Write the problem into an mps file of the given filename */
    virtual void writeMps(const char *filename,
			  const char *extension = "mps") const;
    /** Write the problem into an mps file of the given filename,
	names may be null.  formatType is
	0 - normal
	1 - extra accuracy 
	2 - IEEE hex (later)

	Returns non-zero on I/O error
    */
    virtual int writeMps(const char *filename, 
		  const char ** rowNames, const char ** columnNames,
		  int formatType=0,int numberAcross=2) const ;
  //@}

  /**@name Message handling (extra for Clp messages).
   Normally I presume you would want the same language.
   If not then you could use underlying model pointer */
  //@{
  /// Set language
  void newLanguage(CoinMessages::Language language);
  void setLanguage(CoinMessages::Language language)
  {newLanguage(language);};
  //@}
  //---------------------------------------------------------------------------

  /**@name Clp specific public interfaces */
  //@{
    /// Get pointer to Clp model
  ClpSimplex * getModelPtr() const ;
  //@}

  //---------------------------------------------------------------------------

  /**@name Constructors and destructors */
  //@{
    /// Default Constructor
    OsiClpSolverInterface ();
    
    /// Clone
    virtual OsiSolverInterface * clone(bool copyData = true) const;
    
    /// Copy constructor 
    OsiClpSolverInterface (const OsiClpSolverInterface &);
    
    /// Borrow constructor - only delete one copy
    OsiClpSolverInterface (ClpSimplex *);

    /// Releases so won't error
    void releaseClp();
    
    /// Assignment operator 
    OsiClpSolverInterface & operator=(const OsiClpSolverInterface& rhs);
    
    /// Destructor 
    virtual ~OsiClpSolverInterface ();
  //@}

  //---------------------------------------------------------------------------

protected:
  ///@name Protected methods
  //@{
    /** Apply a row cut (append to constraint matrix). */
    virtual void applyRowCut(const OsiRowCut& rc);

    /** Apply a column cut (adjust one or more bounds). */
    virtual void applyColCut(const OsiColCut& cc);
  //@}

  //---------------------------------------------------------------------------

private:
  /**@name Private methods */
  //@{
    /// The real work of a copy constructor (used by copy and assignment)
    void gutsOfDestructor();
  
    /// Deletes all mutable stuff
    void freeCachedResults() const;

    /// A method that fills up the rowsense_, rhs_ and rowrange_ arrays
    void extractSenseRhsRange() const;

    ///
    void fillParamMaps();
    /// Warm start
    CoinWarmStartBasis getBasis(ClpSimplex * model) const;
    /// Sets up working basis as a copy of input
    void setBasis( const CoinWarmStartBasis & basis, ClpSimplex * model);
  //@}
  
  /**@name Private member data */
  //@{
     /// Clp model represented by this class instance
     mutable ClpSimplex * modelPtr_;

    /**@name Cached information derived from the OSL model */
    //@{
      /// Pointer to dense vector of row sense indicators
      mutable char    *rowsense_;
  
      /// Pointer to dense vector of row right-hand side values
      mutable double  *rhs_;
  
     /** Pointer to dense vector of slack upper bounds for range 
         constraints (undefined for non-range rows)
     */
     mutable double  *rowrange_;

     /** A pointer to the warmstart information to be used in the hotstarts.
	 This is NOT efficient and more thought should be given to it... */
     mutable CoinWarmStartBasis* ws_;
     /** Warmstart information to be used in resolves. */
     CoinWarmStartBasis basis_;
     /** The original iteration limit before hotstarts started. */
     int itlimOrig_;

     /// Last algorithm used
     int lastAlgorithm_;

     /// To say if destructor should delete underlying model
     bool notOwned_;
  
     /// Pointer to row-wise copy of problem matrix coefficients.
     mutable CoinPackedMatrix *matrixByRow_;  

     /// Pointer to integer information
     char * integerInformation_;

     std::map<OsiIntParam, ClpIntParam> intParamMap_;
     std::map<OsiDblParam, ClpDblParam> dblParamMap_;
     std::map<OsiStrParam, ClpStrParam> strParamMap_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the OsiClpSolverInterface class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. Also, if this method is compiled with
    optimization, the compilation takes 10-15 minutes and the machine pages
    (has 256M core memory!)... */
void
OsiClpSolverInterfaceUnitTest(const std::string & mpsDir);

#endif
