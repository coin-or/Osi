#ifndef OsiSolution_H
#define OsiSolution_H

#include <iostream>

class OsiPackedMatrix;
class OsiWarmStartBasis;
class OsiFactorization;


/** This deals with results of solving for a given basis.

    Returns solutions (can also check against given)
    If not feasible gives more information

    It is really just for debug.  It does not have a seperate
    test suite as this is done in OsiFactorizationTest.

    Oh dear, oh dear - the cleanest way to do this is to have a model
    class to contain stuff - this could get dangerous!

    All coding trusts users to behave and not modify model 

*/
/** This struct says which bits of model are in what format and who they
    belong to */
typedef struct {
  unsigned int rhsSense:1; // If rhs in sense format
  unsigned int rhsMine:1; // If I own rhs (and it will be in lo,up format!) 
  unsigned int matrixType:2; //0 - packed, 1 no gaps, 2 gaps 
  unsigned int matrixMine:1; // If I own (and it will be type 1) 
  unsigned int columnMine:1; // If I own column stuff
  unsigned int basisMine:1;  // If I own basis information
  unsigned int factorizationMine:1;  // If I own factorization
  unsigned int extraInfo:24;
} OSIDescription;

class OsiSimplexModel {

public:

  /**@name Constructors and destructor 
     NOTE - At present there are no copy operators to remove 
     ambiguity on who owns what.  So you have to reconstruct model.
   */
  //@{
  /// Default constructor
    OsiSimplexModel (  );

  /// Destructor
   ~OsiSimplexModel (  );
  //@}

  /**@name Borrow model - borrows some stuff and initializes others */
  //@{
    /** Borrows a problem (the constraints on the
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
  void borrowModel (  const OsiPackedMatrix& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub);

    /** Borrows a problem (the constraints on the
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
  void borrowModel (  const OsiPackedMatrix& matrix,
		     const double* collb, const double* colub,
		     const double* obj,
		     const char* rowsen, const double* rowrhs,   
		     const double* rowrng);
  /** Just like the other borrowModel() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
  void borrowModel (  const int numcols, const int numrows,
		     const int* start, const int* index,
		     const double* value,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub);

  void borrowModel (  const int numcols, const int numrows,
		     const int* start, const int* index,
		     const double* value,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const char* rowsen, const double* rowrhs,   
		     const double* rowrng);
  /// Factorizes a problem 
  int factorize ( OsiFactorization & factorization, 
		  const OsiWarmStartBasis & basis);

  /** Factorizes a problem and sets non-basic solution
      Returns +1 if any super-basic variables in warmstart (and okay),
      returns +2 if any super-basic variables not free
      otherwise factorization status */

  int factorizeAndSet ( OsiFactorization & factorization, 
		  const OsiWarmStartBasis & basis);

  //@}

  /**@name gets and sets */
  //@{ 
  /// Primal tolerance to use
  inline double primalTolerance() const 
          { return primalTolerance_;} ;
  void setPrimalTolerance( double value) ;
  /// Dual tolerance to use
  inline double dualTolerance() const 
          { return dualTolerance_;} ;
  void setDualTolerance( double value) ;
  /// Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore
  inline int optimizationDirection() const
          { return optimizationDirection_;};
  void setOptimizationDirection(int value);
  /// Primal row solution
  inline const double * primalRowSolution() const 
          { return primalSolution_[0];} ;
  /// Primal column solution
  inline const double * primalColumnSolution() const 
          { return primalSolution_[1];} ;
  /// Dual row solution
  inline const double * dualRowSolution() const 
          { return dualSolution_[0];} ;
  /// Reduced costs
  inline const double * dualColumnSolution() const 
          { return dualSolution_[1];} ;
  /// How many iterative refinements to do
  inline int numberRefinements() const 
          { return numberRefinements_;} ;
  void setnumberRefinements( int value) ;
  /// factorization 
  inline OsiFactorization * factorization() const 
          { return factorization_;};
  //@}

  /**@private or protected methods */
  //@{
private:
  /// Does most of deletion
  void gutsOfDelete();
protected:
  /// gets lower and upper bounds on rows
  void getRowBound(int iRow, double& lower, double& upper) const;
  //@}


////////////////// data //////////////////
protected:

  /**@name data */
  //@{
  /// Primal tolerance to use
  double primalTolerance_;
  /// Dual tolerance to use
  double dualTolerance_;
  /// Direction of optimization (1 - minimize, -1 - maximize, 0 - ignore
  int optimizationDirection_;
  /// How many iterative refinements to do
  int numberRefinements_;
  /// Number of rows
  int numberRows_;
  /// Number of columns
  int numberColumns_;
  /// Primal arrays - 0 = rows, 1 = columns
  double * primalSolution_[2];
  /// Dual arrays - 0 = rows, 1 = columns
  double * dualSolution_[2];
  /// Basic variables pivoting on which rows
  int * pivotVariable_;
  /// factorization 
  OsiFactorization * factorization_;
  /// Gives information on storage methods
  OSIDescription description_;
  /// Rowsense (if used)
  const char* rowSense_;
  /// Row lower or rhs
  const double* rowLower_;
  /// Row upper or range
  const double* rowUpper_;
  /// Objective
  const double * objective_;
  /// Column Lower
  const double * colLower_;
  /// Column Upper
  const double * colUpper_;
  /// Packed matrix if given
  const OsiPackedMatrix * matrix_;
  /// Column starts if given
  const int * columnStart_;
  /// Column lengths if given
  const int * columnLength_;
  /// Row indices if given
  const int * row_;
  /// Elements if given
  const double * element_;
  /// Warm start
  const OsiWarmStartBasis * basis_;
  //@}
};


class OsiSolution : public OsiSimplexModel {

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    OsiSolution (  );

  /// Destructor
   ~OsiSolution (  );
  //@}

  /**@name Create Solution (returns status of factorization) */
  //@{
    /** Solves a problem 
    */
  int getSolution (  const double * rowActivities,
		     const double * columnActivities);
  int getSolution ();
  //@}

  /**@name gets and sets */
  //@{ 
  /// Worst column primal infeasibility
  inline double columnPrimalInfeasibility() const 
          { return columnPrimalInfeasibility_;} ;
  /// Sequence of worst (-1 if feasible)
  inline int columnPrimalSequence() const 
          { return columnPrimalSequence_;} ;
  /// Worst row primal infeasibility
  inline double rowPrimalInfeasibility() const 
          { return rowPrimalInfeasibility_;} ;
  /// Sequence of worst (-1 if feasible)
  inline int rowPrimalSequence() const 
          { return rowPrimalSequence_;} ;
  /// If problem is primal feasible
  inline bool primalFeasible() const
          { return (rowPrimalSequence_==-1&&columnPrimalSequence_==-1);};
  /** Worst column dual infeasibility (note - these may not be as meaningful
      if the problem is primal infeasible */
  inline double columnDualInfeasibility() const 
          { return columnDualInfeasibility_;} ;
  /// Sequence of worst (-1 if feasible)
  inline int columnDualSequence() const 
          { return columnDualSequence_;} ;
  /// Worst row dual infeasibility
  inline double rowDualInfeasibility() const 
          { return rowDualInfeasibility_;} ;
  /// Sequence of worst (-1 if feasible)
  inline int rowDualSequence() const 
          { return rowDualSequence_;} ;
  /// If problem is dual feasible
  inline bool dualFeasible() const
          { return (rowDualSequence_==-1&&columnDualSequence_==-1);};
  /// Primal tolerance needed to make dual feasible (<largeTolerance)
  inline double primalToleranceToGetOptimal() const 
          { return primalToleranceToGetOptimal_;} ;
  /// Remaining largest dual infeasibility
  inline double dualInfeasibilityWithoutFree() const 
          { return dualInfeasibilityWithoutFree_;} ;
  /// Large bound value (for complementarity etc)
  inline double largeValue() const 
          { return largeValue_;} ;
  void setLargeValue( double value) ;
  /// Largest error on Ax-b
  inline double largestPrimalError() const
          { return largestPrimalError_;} ;
  /// Largest error on basic duals
  inline double largestDualError() const
          { return largestDualError_;} ;
  /// Largest difference between input primal solution and computed
  inline double largestSolutionError() const
          { return largestSolutionError_;} ;
  /// Status of solution i.e. there or not there
  inline bool gotSolution() const 
          { return gotSolution_;} ;
  //@}
  private:
  /**@private methods */
  //@{
  void gutsOfSolution ( const double * rowActivities,
			const double * columnActivities,
			const int numberRows,
			const int numberColumns, 
			const int* columnStart, 
			const int * columnLength,
			const int* row,
			const double * element,
			const double* collb, const double* colub,   
			const double* obj,
			const double* rowlb, const double* rowub);
  //@}

////////////////// data //////////////////
private:

  /**@name data */
  //@{
  /// Worst column primal infeasibility
  double columnPrimalInfeasibility_;
  /// Sequence of worst (-1 if feasible)
  int columnPrimalSequence_;
  /// Worst row primal infeasibility
  double rowPrimalInfeasibility_;
  /// Sequence of worst (-1 if feasible)
  int rowPrimalSequence_;
  /// Worst column dual infeasibility
  double columnDualInfeasibility_;
  /// Sequence of worst (-1 if feasible)
  int columnDualSequence_;
  /// Worst row dual infeasibility
  double rowDualInfeasibility_;
  /// Sequence of worst (-1 if feasible)
  int rowDualSequence_;
  /// Primal tolerance needed to make dual feasible (<largeTolerance)
  double primalToleranceToGetOptimal_;
  /// Remaining largest dual infeasibility
  double dualInfeasibilityWithoutFree_;
  /// Large bound value (for complementarity etc)
  double largeValue_;
  /// Largest error on Ax-b
  double largestPrimalError_;
  /// Largest error on basic duals
  double largestDualError_;
  /// Largest difference between input primal solution and computed
  double largestSolutionError_;
  /// Status of solution i.e. there or not there
  bool gotSolution_;
  //@}
};
#endif
