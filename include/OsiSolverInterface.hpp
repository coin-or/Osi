// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiSolverInterface_H
#define OsiSolverInterface_H

#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "CoinPackedVectorBase.hpp"

#include "OsiCollections.hpp"
#include "OsiSolverParameters.hpp"

class CoinPackedMatrix;
class CoinWarmStart;

class OsiCuts;
class OsiRowCut;
class OsiRowCutDebugger;
class CoinSet;
#include "CoinFinite.hpp"

//#############################################################################

/** Solver Interface Abstract Base Class

  Abstract Base Class for describing an interface to a solver.

  Many OsiSolverInterface query methods return a const pointer to the
  requested read-only data. If the model data is changed or the solver
  is called, these pointers may no longer be valid and should be 
  refreshed by invoking the member function to obtain an updated copy
  of the pointer.
  For example:
  \code
      OsiSolverInterface solverInterfacePtr ;
      const double * ruBnds = solverInterfacePtr->getRowUpper();
      solverInterfacePtr->applyCuts(someSetOfCuts);
      // ruBnds is no longer a valid pointer and must be refreshed
      ruBnds = solverInterfacePtr->getRowUpper();
  \endcode

  Querying a problem that has no data associated with it will result in
  zeros for the number of rows and columns, and NULL pointers from
  the methods that return vectors.
*/

class OsiSolverInterface  {
   friend void OsiSolverInterfaceCommonUnitTest(
      const OsiSolverInterface* emptySi,
      const std::string & mpsDir,
      const std::string & netlibDir);
   friend void OsiSolverInterfaceMpsUnitTest(
      const std::vector<OsiSolverInterface*> & vecSiP,
      const std::string & mpsDir);

public:
  /// Internal class for obtaining status from the applyCuts method 
  class ApplyCutsReturnCode {
    friend class OsiSolverInterface;
    friend class OsiOslSolverInterface;

  public:
    ///@name Constructors and desctructors
    //@{
      /// Default constructor
      ApplyCutsReturnCode():
	 intInconsistent_(0),
	 extInconsistent_(0),
	 infeasible_(0),
	 ineffective_(0),
	 applied_(0) {} 
      /// Copy constructor
      ApplyCutsReturnCode(const ApplyCutsReturnCode & rhs):
	 intInconsistent_(rhs.intInconsistent_),
	 extInconsistent_(rhs.extInconsistent_),
	 infeasible_(rhs.infeasible_),
	 ineffective_(rhs.ineffective_),
	 applied_(rhs.applied_) {} 
      /// Assignment operator
      ApplyCutsReturnCode & operator=(const ApplyCutsReturnCode& rhs)
      { 
	if (this != &rhs) { 
	  intInconsistent_ = rhs.intInconsistent_;
	  extInconsistent_ = rhs.extInconsistent_;
	  infeasible_      = rhs.infeasible_;
	  ineffective_     = rhs.ineffective_;
	  applied_         = rhs.applied_;
	}
	return *this;
      }
      /// Destructor
      ~ApplyCutsReturnCode(){}
    //@}

    /**@name Accessing return code attributes */
    //@{
      /// Number of logically inconsistent cuts
      inline int getNumInconsistent(){return intInconsistent_;}
      /// Number of cuts inconsistent with the current model
      inline int getNumInconsistentWrtIntegerModel(){return extInconsistent_;}
      /// Number of cuts that cause obvious infeasibility
      inline int getNumInfeasible(){return infeasible_;}
      /// Number of redundant or ineffective cuts
      inline int getNumIneffective(){return ineffective_;}
      /// Number of cuts applied
      inline int getNumApplied(){return applied_;}
    //@}

  private: 
    /**@name Private methods */
    //@{
      /// Increment logically inconsistent cut counter 
      inline void incrementInternallyInconsistent(){intInconsistent_++;}
      /// Increment model-inconsistent counter
      inline void incrementExternallyInconsistent(){extInconsistent_++;}
      /// Increment infeasible cut counter
      inline void incrementInfeasible(){infeasible_++;}
      /// Increment ineffective cut counter
      inline void incrementIneffective(){ineffective_++;}
      /// Increment applied cut counter
      inline void incrementApplied(){applied_++;}
    //@}

    ///@name Private member data
    //@{
      /// Counter for logically inconsistent cuts
      int intInconsistent_;
      /// Counter for model-inconsistent cuts
      int extInconsistent_;
      /// Counter for infeasible cuts
      int infeasible_;
      /// Counter for ineffective cuts
      int ineffective_;
      /// Counter for applied cuts
      int applied_;
    //@}
  };

  //---------------------------------------------------------------------------

public:
  ///@name Solve methods 
  //@{
    /// Solve initial LP relaxation 
    virtual void initialSolve() = 0; 

    /// Resolve an LP relaxation after problem modification
    virtual void resolve() = 0;

    /// Invoke solver's built-in enumeration algorithm
    virtual void branchAndBound() = 0;
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

     \note
     There is a default implementation of the set/get
     methods, namely to store/retrieve the given value using an array in the
     base class. A specific solver implementation can use this feature, for
     example, to store parameters that should be used later on. Implementors
     of a solver interface should overload these functions to provide the
     proper interface to and accurately reflect the capabilities of a
     specific solver.

     The format for hints is slightly different in that the value is 
     boolean and there is an enum to show strength of hint.
     There is also an optional void pointer to allow for any eventuality.
     Hints should be initialised when a solver is instantiated.
     (See OsiSolverParameters.hpp for defined hint parameters and strength.)
     A value of true means to work with the hint, false to work against it.
     For example,
     <ul>
       <li> \code setHintParam(OsiDoScale,true,OsiHintTry) \endcode
	    is a mild suggestion to the solver to scale the constraint
	    system.
       <li> \code setHintParam(OsiDoScale,false,OsiForceDo) \endcode
	    tells the solver to disable scaling, or throw an exception if
	    it cannot comply.
     </ul>
     As another example, a solver interface could use the value and strength
     of the \c OsiDoReducePrint hint to adjust the amount of information
     printed by the interface and/or solver.  The extent to which a solver
     obeys hints is left to the solver.  The value and strength returned by
     \c getHintParam will match the most recent call to \c setHintParam,
     and will not necessarily reflect the solver's ability to comply with the
     hint.  If the hint strength is \c OsiForceDo, the solver is required to
     throw an exception if it cannot perform the specified action.

     \note
     As with the other set/get methods, there is a default implementation
     which maintains arrays in the base class for hint value and strength.
     The default implementation does not store the void pointer, and always
     throws an exception for strength \c OsiForceDo. Implementors of a solver
     interface should overload these functions to provide the proper interface
     to and accurately reflect the capabilities of a specific solver.
  */
  //@{
    // Set an integer parameter
    virtual bool setIntParam(OsiIntParam key, int value) {
      if (key == OsiLastIntParam) return (false) ;
      intParam_[key] = value;
      return true;
    }
    // Set an double parameter
    virtual bool setDblParam(OsiDblParam key, double value) {
      if (key == OsiLastDblParam) return (false) ;
      dblParam_[key] = value;
      return true;
    }
    // Set an string parameter
    virtual bool setStrParam(OsiStrParam key, const std::string & value) {
      if (key == OsiLastStrParam) return (false) ;
      strParam_[key] = value;
      return true;
    }
    // Set a hint parameter
    virtual bool setHintParam(OsiHintParam key, bool yesNo=true,
			      OsiHintStrength strength=OsiHintTry,
			      void * otherInformation=NULL) {
      if (key==OsiLastHintParam)
	return false; 
      hintParam_[key] = yesNo;
      hintStrength_[key] = strength;
      if (strength == OsiForceDo)
	throw CoinError("OsiForceDo illegal",
			"setHintParam", "OsiSolverInterface");
      return true;
    }
    // Get an integer parameter
    virtual bool getIntParam(OsiIntParam key, int& value) const {
      if (key == OsiLastIntParam) return (false) ;
      value = intParam_[key];
      return true;
    }
    // Get an double parameter
    virtual bool getDblParam(OsiDblParam key, double& value) const {
      if (key == OsiLastDblParam) return (false) ;
      value = dblParam_[key];
      return true;
    }
    // Get a string parameter
    virtual bool getStrParam(OsiStrParam key, std::string& value) const {
      if (key == OsiLastStrParam) return (false) ;
      value = strParam_[key];
      return true;
    }
    // get a hint parameter
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength,
			      void *& otherInformation) const {
      if (key==OsiLastHintParam)
	return false; 
      yesNo = hintParam_[key];
      strength = hintStrength_[key];
      otherInformation=NULL;
      return true;
    }
    // get a hint parameter (less information)
    virtual bool getHintParam(OsiHintParam key, bool& yesNo,
			      OsiHintStrength& strength) const {
      if (key==OsiLastHintParam)
	return false; 
      yesNo = hintParam_[key];
      strength = hintStrength_[key];
      return true;
    }
    // get a hint parameter (even less information)
    virtual bool getHintParam(OsiHintParam key, bool& yesNo) const {
      if (key==OsiLastHintParam)
	return false; 
      yesNo = hintParam_[key];
      return true;
    }
    // copy all parameters in this section from one solver to another
    void copyParameters(OsiSolverInterface & rhs);
  //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there numerical difficulties?
    virtual bool isAbandoned() const = 0;
    /// Is optimality proven?
    virtual bool isProvenOptimal() const = 0;
    /// Is primal infeasiblity proven?
    virtual bool isProvenPrimalInfeasible() const = 0;
    /// Is dual infeasiblity proven?
    virtual bool isProvenDualInfeasible() const = 0;
    /// Is the given primal objective limit reached?
    virtual bool isPrimalObjectiveLimitReached() const = 0;
    /// Is the given dual objective limit reached?
    virtual bool isDualObjectiveLimitReached() const = 0;
    /// Iteration limit reached?
    virtual bool isIterationLimitReached() const = 0;
  //@}

  //---------------------------------------------------------------------------
  /**@name Warm start methods */
  //@{
    /*! \brief Get an empty warm start object
      
      This routine returns an empty warm start object. Its purpose is
      to provide a way to give a client a warm start object of the
      appropriate type, which can resized and modified as desired.
    */

    virtual CoinWarmStart *getEmptyWarmStart () const = 0 ;

    /** Get warm start information.

      If there is no valid solution, an empty warm start object (0 rows, 0
      columns) wil be returned.
    */
    virtual CoinWarmStart* getWarmStart() const = 0;

    /** Set warm start information.
    
      Return true/false depending on whether the warm start information was
      accepted or not. */
    virtual bool setWarmStart(const CoinWarmStart* warmstart) = 0;
  //@}

  //---------------------------------------------------------------------------
  /**@name Hot start methods
  
     Primarily used in strong branching. The user can create a hot start
     object --- a snapshot of the optimization process --- then reoptimize
     over and over again, starting from the same point.

     \note
     <ul>
     <li> Between hot started optimizations only bound changes are allowed.
     <li> The copy constructor and assignment operator should NOT copy any
          hot start information.
     <li> The default implementation simply extracts a warm start object in
          \c markHotStart, resets to the warm start object in
	  \c solveFromHotStart, and deletes the warm start object in
	  \c unmarkHotStart.
	  <em>Actual solver implementations are encouraged to do better.</em>
     </ul>

  */
  //@{
    /// Create a hot start snapshot of the optimization process.
    virtual void markHotStart();
    /// Optimize starting from the hot start snapshot.
    virtual void solveFromHotStart();
    /// Delete the hot start snapshot.
    virtual void unmarkHotStart();
  //@}

  //---------------------------------------------------------------------------
    /**@name Problem query methods

     Querying a problem that has no data associated with it will result in
     zeros for the number of rows and columns, and NULL pointers from
     the methods that return vectors.
     
     Const pointers returned from any data-query method are valid as
     long as the data is unchanged and the solver is not called.
    */
    //@{
      /// Get number of columns
      virtual int getNumCols() const = 0;
  
      /// Get number of rows
      virtual int getNumRows() const = 0;
  
      /// Get number of nonzero elements
      virtual int getNumElements() const = 0;
  
      /// Get pointer to array[getNumCols()] of column lower bounds
      virtual const double * getColLower() const = 0;
  
      /// Get pointer to array[getNumCols()] of column upper bounds
      virtual const double * getColUpper() const = 0;
  
      /** Get pointer to array[getNumRows()] of row constraint senses.
  	<ul>
  	<li>'L': <= constraint
  	<li>'E': =  constraint
  	<li>'G': >= constraint
  	<li>'R': ranged constraint
  	<li>'N': free constraint
  	</ul>
      */
      virtual const char * getRowSense() const = 0;
  
      /** Get pointer to array[getNumRows()] of row right-hand sides
  	<ul>
  	  <li> if getRowSense()[i] == 'L' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'G' then
	       getRightHandSide()[i] == getRowLower()[i]
  	  <li> if getRowSense()[i] == 'R' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'N' then
	       getRightHandSide()[i] == 0.0
  	</ul>
      */
      virtual const double * getRightHandSide() const = 0;
  
      /** Get pointer to array[getNumRows()] of row ranges.
  	<ul>
            <li> if getRowSense()[i] == 'R' then
                    getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
            <li> if getRowSense()[i] != 'R' then
                    getRowRange()[i] is 0.0
          </ul>
      */
      virtual const double * getRowRange() const = 0;
  
      /// Get pointer to array[getNumRows()] of row lower bounds
      virtual const double * getRowLower() const = 0;
  
      /// Get pointer to array[getNumRows()] of row upper bounds
      virtual const double * getRowUpper() const = 0;
  
      /// Get pointer to array[getNumCols()] of objective function coefficients
      virtual const double * getObjCoefficients() const = 0;
  
      /// Get objective function sense (1 for min (default), -1 for max)
      virtual double getObjSense() const = 0;
  
      /// Return true if variable is continuous
      virtual bool isContinuous(int colIndex) const = 0;
  
      /// Return true if variable is binary
      virtual bool isBinary(int colIndex) const;
  
      /** Return true if column is integer.
          Note: This function returns true if the the column
          is binary or a general integer.
      */
      virtual bool isInteger(int colIndex) const;
  
      /// Return true if variable is general integer
      virtual bool isIntegerNonBinary(int colIndex) const;
  
      /// Return true if variable is binary and not fixed at either bound
      virtual bool isFreeBinary(int colIndex) const; 
    
      /// Get pointer to row-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByRow() const = 0;
  
      /// Get pointer to column-wise copy of matrix
      virtual const CoinPackedMatrix * getMatrixByCol() const = 0;
  
      /// Get solver's value for infinity
      virtual double getInfinity() const = 0;
    //@}
    
    /**@name Solution query methods */
    //@{
      /// Get pointer to array[getNumCols()] of primal variable values
      virtual const double * getColSolution() const = 0;
  
      /// Get pointer to array[getNumRows()] of dual variable values
      virtual const double * getRowPrice() const = 0;
  
      /// Get a pointer to array[getNumCols()] of reduced costs
      virtual const double * getReducedCost() const = 0;
  
      /** Get pointer to array[getNumRows()] of row activity levels (constraint
  	matrix times the solution vector). */
      virtual const double * getRowActivity() const = 0;
  
      /// Get objective function value
      virtual double getObjValue() const = 0;

      /** Get the number of iterations it took to solve the problem (whatever
	  ``iteration'' means to the solver). */
      virtual int getIterationCount() const = 0;
  
      /** Get as many dual rays as the solver can provide. In case of proven
          primal infeasibility there should be at least one.
     
          \note
	  Implementors of solver interfaces note that
          the double pointers in the vector should point to arrays of length
          getNumRows() and they should be allocated via new[].
     
          \note
	  Clients of solver interfaces note that
          it is the client's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getDualRays(int maxNumRays) const = 0;
      /** Get as many primal rays as the solver can provide. (In case of proven
          dual infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumCols() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getPrimalRays(int maxNumRays) const = 0;
  
      /** Get vector of indices of primal variables which are integer variables 
  	  but have fractional values in the current solution. */
      virtual OsiVectorInt getFractionalIndices(const double etol=1.e-05)
	const;
    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to modify the objective, bounds, and solution

       For functions which take a set of indices as parameters
       (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
       \c setRowSetTypes()), the parameters follow the C++ STL iterator
       convention: \c indexFirst points to the first index in the
       set, and \c indexLast points to a position one past the last index
       in the set.
    
    */
    //@{
      /** Set an objective function coefficient */
      virtual void setObjCoeff( int elementIndex, double elementValue ) = 0;

      /** Set a set of objective function coefficients */
      virtual void setObjCoeffSet(const int* indexFirst,
				  const int* indexLast,
				  const double* coeffList);

      /** Set a single column lower bound.
    	  Use -getInfinity() for -infinity. */
      virtual void setColLower( int elementIndex, double elementValue ) = 0;
      
      /** Set a single column upper bound.
    	  Use getInfinity() for infinity. */
      virtual void setColUpper( int elementIndex, double elementValue ) = 0;
      
      /** Set a single column lower and upper bound.
    	  The default implementation just invokes setColLower() and
    	  setColUpper() */
      virtual void setColBounds( int elementIndex,
    				 double lower, double upper ) {
    	 setColLower(elementIndex, lower);
    	 setColUpper(elementIndex, upper);
      }
    
      /** Set the upper and lower bounds of a set of columns.
    	  The default implementation just invokes setColBounds() over
	  and over again.
	  For each column, boundList must contain both a lower and
	  upper bound, in that order.
      */
      virtual void setColSetBounds(const int* indexFirst,
				   const int* indexLast,
				   const double* boundList);
      
      /** Set a single row lower bound.
    	  Use -getInfinity() for -infinity. */
      virtual void setRowLower( int elementIndex, double elementValue ) = 0;
      
      /** Set a single row upper bound.
    	  Use getInfinity() for infinity. */
      virtual void setRowUpper( int elementIndex, double elementValue ) = 0;
    
      /** Set a single row lower and upper bound.
    	  The default implementation just invokes setRowLower() and
    	  setRowUpper() */
      virtual void setRowBounds( int elementIndex,
    				 double lower, double upper ) {
    	 setRowLower(elementIndex, lower);
    	 setRowUpper(elementIndex, upper);
      }
    
      /** Set the type of a single row */
      virtual void setRowType(int index, char sense, double rightHandSide,
    			      double range) = 0;
    
      /** Set the bounds on a set of rows.
    	  The default implementation just invokes setRowBounds()
    	  over and over again.
      */
      virtual void setRowSetBounds(const int* indexFirst,
    				   const int* indexLast,
    				   const double* boundList);
    
      /** Set the type of a set of rows.
    	  The default implementation just invokes setRowType()
    	  over and over again.
      */
      virtual void setRowSetTypes(const int* indexFirst,
				  const int* indexLast,
				  const char* senseList,
				  const double* rhsList,
				  const double* rangeList);
    
    /// Set the objective function sense.
    /// (1 for min (default), -1 for max)
    virtual void setObjSense(double s) = 0;
    
    /** Set the primal solution variable values
    
	colsol[getNumCols()] is an array of values for the primal variables.
	These values are copied to memory owned by the solver interface object
	or the solver.  They will be returned as the result of getColSolution()
	until changed by another call to setColSolution() or by a call to any
	solver routine.  Whether the solver makes use of the solution in any
	way is solver-dependent.
    */
    virtual void setColSolution(const double *colsol) = 0;

    /** Set dual solution variable values

	rowprice[getNumRows()] is an array of values for the dual
	variables. These values are copied to memory owned by the solver
	interface object or the solver.  They will be returned as the result of
	getRowPrice() until changed by another call to setRowPrice() or by a
	call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */

   virtual void setRowPrice(const double * rowprice) = 0;

    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to set variable type */
    //@{
      /** Set the index-th variable to be a continuous variable */
      virtual void setContinuous(int index) = 0;
      /** Set the index-th variable to be an integer variable */
      virtual void setInteger(int index) = 0;
      /** Set the variables listed in indices (which is of length len) to be
	  continuous variables */
      virtual void setContinuous(const int* indices, int len);
      /** Set the variables listed in indices (which is of length len) to be
	  integer variables */
      virtual void setInteger(const int* indices, int len);
    //@}
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    /**@name Methods to expand a problem.

       Note that new columns are added as continuous variables.

    */
    //@{
      /** Add a column (primal variable) to the problem. */
      virtual void addCol(const CoinPackedVectorBase& vec,
			  const double collb, const double colub,   
			  const double obj) = 0;
      /** Add a set of columns (primal variables) to the problem.
      
	The default implementation simply makes repeated calls to
	addCol().
      */
      virtual void addCols(const int numcols,
			   const CoinPackedVectorBase * const * cols,
			   const double* collb, const double* colub,   
			   const double* obj);
      /** Add a column (primal variable) to the problem. */
      virtual void addCol(int numberElements, const int * rows, const double * elements,
			  const double collb, const double colub,   
			  const double obj) ;
      /** Add a set of columns (primal variables) to the problem.
      
	The default implementation simply makes repeated calls to
	addCol().
      */
      virtual void addCols(const int numcols,
			   const int * columnStarts, const int * rows, const double * elements,
			   const double* collb, const double* colub,   
			   const double* obj);
#if 0
      /** */
      virtual void addCols(const CoinPackedMatrix& matrix,
			   const double* collb, const double* colub,   
			   const double* obj);
#endif
      /** Remove a set of columns (primal variables) from the problem.  */
      virtual void deleteCols(const int num, const int * colIndices) = 0;
    
      /** Add a row (constraint) to the problem. */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const double rowlb, const double rowub) = 0;
      /** */
      virtual void addRow(const CoinPackedVectorBase& vec,
    			  const char rowsen, const double rowrhs,   
    			  const double rowrng) = 0;
      /** Add a set of rows (constraints) to the problem.
      
	The default implementation simply makes repeated calls to
	addRow().
      */
      virtual void addRows(const int numrows,
			   const CoinPackedVectorBase * const * rows,
			   const double* rowlb, const double* rowub);
      /** Add a set of rows (constraints) to the problem.
      
	The default implementation simply makes repeated calls to
	addRow().
      */
      virtual void addRows(const int numrows,
			   const CoinPackedVectorBase * const * rows,
    			   const char* rowsen, const double* rowrhs,   
    			   const double* rowrng);
#if 0
      /** */
      virtual void addRows(const CoinPackedMatrix& matrix,
    			   const double* rowlb, const double* rowub);
      /** */
      virtual void addRows(const CoinPackedMatrix& matrix,
    			   const char* rowsen, const double* rowrhs,   
    			   const double* rowrng);
#endif
      /** Delete a set of rows (constraints) from the problem. */
      virtual void deleteRows(const int num, const int * rowIndices) = 0;
    
      //-----------------------------------------------------------------------
      /** Apply a collection of cuts.

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
      /** Apply a collection of row cuts which are all effective.
	  applyCuts seems to do one at a time which seems inefficient.
	  Would be even more efficient to pass an array of pointers.
      */
      virtual void applyRowCuts(int numberCuts, const OsiRowCut * cuts);
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
			     const double* rowlb, const double* rowub) = 0;
			    
    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by lower and upper bounds).
	For default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       double*& rowlb, double*& rowub) = 0;

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
			     const double* rowrng) = 0;

    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by sense/rhs/range triplets). For
        default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       char*& rowsen, double*& rowrhs,
			       double*& rowrng) = 0;

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub) = 0;

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng) = 0;

    /** Read a problem in MPS format from the given filename.
    
	The default implementation uses CoinMpsIO::readMps() to read
	the MPS file and returns the number of errors encountered.
   */
    virtual int readMps(const char *filename,
			 const char *extension = "mps") ;

    /** Read a problem in MPS format from the given full filename.
    
	This uses CoinMpsIO::readMps() to read
	the MPS file and returns the number of errors encountered.
	It also may return an array of set information
   */
  virtual int readMps(const char *filename, const char*extension,
			int & numberSets, CoinSet ** & sets);

    /** Write the problem in MPS format to the specified file.

      If objSense is non-zero, a value of -1.0 causes the problem to be
      written with a maximization objective; +1.0 forces a minimization
      objective. If objSense is zero, the choice is left to implementation.
    */
    virtual void writeMps(const char *filename,
			  const char *extension = "mps",
			  double objSense=0.0) const = 0;

    /** Write the problem in MPS format to the specified file.

	Row and column names may be null.
	formatType is
	<ul>
	  <li> 0 - normal
	  <li> 1 - extra accuracy 
	  <li> 2 - IEEE hex (later)
	</ul>

	Returns non-zero on I/O error
    */
    int writeMpsNative(const char *filename, 
		  const char ** rowNames, const char ** columnNames,
		  int formatType=0,int numberAcross=2,
		 double objSense=0.0) const ;
  //@}

  //---------------------------------------------------------------------------

  /**@name Setting/Accessing application data */
  //@{
    /** Set application data.

	This is a pointer that the application can store into and
	retrieve from the solver interface.
	This field is available for the application to optionally
	define and use.
    */
    void setApplicationData (void * appData);

    /// Get application data
    void * getApplicationData() const;
  //@}
  //---------------------------------------------------------------------------

  /**@name Message handling
  
    See the COIN library documentation for additional information about
    COIN message facilities.
  
  */
  //@{
  /** Pass in a message handler
  
    It is the client's responsibility to destroy a message handler installed
    by this routine; it will not be destroyed when the solver interface is
    destroyed. 
  */
  void passInMessageHandler(CoinMessageHandler * handler);
  /// Set language
  void newLanguage(CoinMessages::Language language);
  void setLanguage(CoinMessages::Language language)
  {newLanguage(language);};
  /// Return a pointer to the current message handler
  CoinMessageHandler * messageHandler() const
  {return handler_;};
  /// Return the current set of messages
  CoinMessages messages() 
  {return messages_;};
  /// Return a pointer to the current set of messages
  CoinMessages * messagesPointer() 
  {return &messages_;};
  //@}
  //---------------------------------------------------------------------------

  /**@name Methods related to testing generated cuts */
  //@{
    /** Activate the row cut debugger.

        If the model name passed is on list of known models
	then all cuts are checked to see that they do NOT cut
	off the known optimal solution.  
    */
    virtual void activateRowCutDebugger (const char * modelName);

    /** Get the row cut debugger.

	If there is a row cut debugger object associated with
	model AND if the known optimal solution is within the
	current feasible region then a pointer to the object is
	returned which may be used to test validity of cuts.

	Otherwise NULL is returned
    */
    const OsiRowCutDebugger * getRowCutDebugger() const;
  //@} 

  //---------------------------------------------------------------------------

  ///@name Constructors and destructors
  //@{
    /// Default Constructor
    OsiSolverInterface(); 
    
    /** Clone

      The result of calling clone(false) is defined to be equivalent to
      calling the default constructor OsiSolverInterface().
    */
    virtual OsiSolverInterface * clone(bool copyData = true) const = 0;
  
    /// Copy constructor 
    OsiSolverInterface(const OsiSolverInterface &);
  
    /// Assignment operator 
    OsiSolverInterface & operator=(const OsiSolverInterface& rhs);
  
    /// Destructor 
    virtual ~OsiSolverInterface ();

    /** Reset the solver interface.

    A call to reset() returns the solver interface to the same state as
    it would have if it had just been constructed by calling the default
    constructor OsiSolverInterface().
    */
    virtual void reset();
  //@}

  //---------------------------------------------------------------------------

protected:
  ///@name Protected methods
  //@{
    /** Apply a row cut (append to the constraint matrix). */
    virtual void applyRowCut( const OsiRowCut & rc ) = 0;

    /** Apply a column cut (adjust the bounds of one or more variables). */
    virtual void applyColCut( const OsiColCut & cc ) = 0;

    /** A quick inlined function to convert from the lb/ub style of
	constraint definition to the sense/rhs/range style */
    inline void
    convertBoundToSense(const double lower, const double upper,
			char& sense, double& right, double& range) const;
    /** A quick inlined function to convert from the sense/rhs/range style
	of constraint definition to the lb/ub style */
    inline void
    convertSenseToBound(const char sense, const double right,
			const double range,
			double& lower, double& upper) const;
    /** A quick inlined function to force a value to be between a minimum and
	a maximum value */
    template <class T> inline T
    forceIntoRange(const T value, const T lower, const T upper) const {
      return value < lower ? lower : (value > upper ? upper : value);
    }
    /** Set OsiSolverInterface object state for default constructor

      This routine establishes the initial values of data fields in the
      OsiSolverInterface object when the object is created using the
      default constructor.
    */
    void setInitialData();
  //@}

  ///@name Protected member data
  //@{
    /// Pointer to row cut debugger object
    OsiRowCutDebugger * rowCutDebugger_;
  //@}
  
  //---------------------------------------------------------------------------

private:
  ///@name Private member data 
  //@{
    /// Pointer to user-defined data structure
    void * appData_;
    /// Array of integer parameters
    int intParam_[OsiLastIntParam];
    /// Array of double parameters
    double dblParam_[OsiLastDblParam];
    /// Array of string parameters
    std::string strParam_[OsiLastStrParam];
    /// Array of hint parameters
    bool hintParam_[OsiLastHintParam];
    /// Array of hint strengths
    OsiHintStrength hintStrength_[OsiLastHintParam];
    /** Warm start information used for hot starts when the default
       hot start implementation is used. */
    CoinWarmStart* ws_;
  // Why not just make useful stuff protected?
protected:
   /// Message handler
  CoinMessageHandler * handler_;
  /** Flag to say if the currrent handler is the default handler.
      Indicates if the solver interface object is responsible
      for destruction of the handler (true) or if the client is
      responsible (false).
  */
  bool defaultHandler_;
  /// Messages
  CoinMessages messages_;
 //@}
};

//#############################################################################
/** A function that tests the methods in the OsiSolverInterface class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. Also, if this method is compiled with
    optimization, the compilation takes 10-15 minutes and the machine pages
    (has 256M core memory!)... */
void
OsiSolverInterfaceCommonUnitTest(
   const OsiSolverInterface* emptySi,
   const std::string & mpsDir,
   const std::string & netlibDir);

//#############################################################################
/** A function that tests that a lot of problems given in MPS files (mostly
    the NETLIB problems) solve properly with all the specified solvers. */
void
OsiSolverInterfaceMpsUnitTest(
   const std::vector<OsiSolverInterface*> & vecSiP,
   const std::string & mpsDir);

//#############################################################################
/** A quick inlined function to convert from the lb/ub style of constraint
    definition to the sense/rhs/range style */
inline void
OsiSolverInterface::convertBoundToSense(const double lower, const double upper,
					char& sense, double& right,
					double& range) const
{
  double inf = getInfinity();
  range = 0.0;
  if (lower > -inf) {
    if (upper < inf) {
      right = upper;
      if (upper==lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < inf) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
/** A quick inlined function to convert from the sense/rhs/range style of
    constraint definition to the lb/ub style */
inline void
OsiSolverInterface::convertSenseToBound(const char sense, const double right,
					const double range,
					double& lower, double& upper) const
{
  double inf=getInfinity();
  switch (sense) {
  case 'E':
    lower = upper = right;
    break;
  case 'L':
    lower = -inf;
    upper = right;
    break;
  case 'G':
    lower = right;
    upper = inf;
    break;
  case 'R':
    lower = right - range;
    upper = right;
    break;
  case 'N':
    lower = -inf;
    upper = inf;
    break;
  }
}

#endif
