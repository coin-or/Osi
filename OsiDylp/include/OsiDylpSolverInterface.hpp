//-----------------------------------------------------------------------------
// Copyright (C) 2002, Lou Hafer, Stephen Tse, International Business Machines
// Corporation and others.  All Rights Reserved.
//-----------------------------------------------------------------------------

#ifdef COIN_USE_DYLP
#ifndef OsiDylpSolverInterface_H
#define OsiDylpSolverInterface_H


/*! \file OsiDylpSolverInterface.hpp
    \brief Declaration of COIN/OSI interface for dylp.

  This file contains the declaration of the class OsiDylpSolverInterface, an
  implementation of the COIN-OSI interface for dylp, the lp solver for the
  bonsaiG MILP code.
*/

/*
  %W%	%G%
*/

#include <OsiPackedMatrix.hpp>
#include <OsiSolverInterface.hpp>
#include <OsiWarmStart.hpp>

#define DYLP_INTERNAL
extern "C" {
#include "dylp.h"
}

/*! \class OsiDylpSolverInterface
    \brief COIN-OSI interface implementation for dylp

  Very roughly: the class OsiDylpSolverInterface (ODSI) is the top level
  object. It implements the public functions defined for an OsiSolverInterface
  object.

  The OSI model (at least as expressed in the test suite) expects that the
  constraint system read back from the solver will be the same as the
  constraint system given to the solver.  Dylp expects that any grooming of
  the constraint system will be done before it's called, and furthermore
  expects that this grooming will include conversion of >= constraints to <=
  constraints. This creates a bit of a conflict. The solution here is to turn
  off constraint system grooming normally done in mpsio.c:mpsin, and move it
  to a wrapper, dylp.c:dylp_dolp. dylp_dolp only implements the conversion
  from >= to <=, which is reversible. dylp can tolerate empty constraints.
  dylp_dolp also embodies a rudimentary strategy that attempts to recover
  from numerical inaccuracy by refactoring the basis more often (on the
  premise that this will reduce numerical inaccuracy in calculations
  involving the basis inverse).

  Within dylp, a constraint system is held in a row- and column-linked
  structure called a consys_struct. The lp problem (constraint system plus
  additional problem information --- primal and dual variables, basis, status,
  etc.) is held in a structure called an lpprob_struct.
  
  Down in the private section of the ODSI class, ODSI.consys is the
  constraint system. ODSI.lpprob is the lpprob_struct that's used to pass the
  problem to the lp solver and return the results.  There are three other
  structures, options, tolerance, and statistics, that respectively hold
  control parameters, tolerances, and collect statistics.
*/

class OsiDylpSolverInterface: public OsiSolverInterface

{ friend void OsiDylpSolverInterfaceUnitTest(const std::string & mpsDir) ;

/*
  Consult the COIN OSI documentation or relevant source code for details
  not covered here. Supported functions are listed first, followed by
  unsupported functions.
*/

public:

/*! \name Constructors and Destructors */
//@{

  /*! \brief Default constructor */

  OsiDylpSolverInterface() ;

  /*! \brief Copy constructor */

  OsiDylpSolverInterface(const OsiDylpSolverInterface& src) ;

  /*! \brief Clone */

  OsiSolverInterface* clone(bool copyData = true) const ;

  /*! \brief Default destructor */

  ~OsiDylpSolverInterface() ;

//@}

/*! \name Methods to load a problem */
//@{

  /*! \brief Read a problem description in MPS format from a file.  */

  int readMps(const char *, const char* = "mps") ;

  /*! \brief Load a problem description (OSI packed matrix, row sense,
	  parameters unaffected).
  */
  void loadProblem(const OsiPackedMatrix&, const double *, const double *, 
    const double *, const char *, const double *, const double*) ;

  /*! \brief Load a problem description (OSI packed matrix, row bounds,
	  parameters unaffected).
  */
  void loadProblem(const OsiPackedMatrix&, const double *, const double *,
    const double *, const double *, const double*) ;
  
  /*! \brief Load a problem description (standard column-major packed
	  matrix, row sense, parameters unaffected)
  */
  void loadProblem (const int colcnt, const int rowcnt,
    const int *start, const int *index, const double *value,
    const double* col_lower, const double* col_upper, const double* obj,
    const char* sense, const double* rhsin, const double* range) ;

  /*! \brief Load a problem description (standard column-major packed
	   matrix, row bounds, parameters unaffected)
  */
  void loadProblem(const int colcnt, const int rowcnt,
    const int *start, const int *index, const double *value,
    const double* col_lower, const double* col_upper, const double* obj,
    const double* row_lower, const double* row_upper) ;

  /*! \brief Load a problem description (OSI packed matrix, row sense,
	  parameters destroyed).
  */
  void assignProblem(OsiPackedMatrix*&, double*&, double*&, double*&, 
    char*&, double*&, double*&) ;

  /*! \brief Load a problem description (OSI packed matrix, row bounds,
	  parameters destroyed).
  */
  void assignProblem(OsiPackedMatrix*&, double*&, double*&, 
    double*&, double*&, double*&) ;

//@}

/*! \name Methods to obtain problem information */

//@{

  /*! \brief Get the number of columns (variables) */

  int getNumCols() const ;

  /*! \brief Get the number of rows (constraints) */

  int getNumRows() const ;

  /*! \brief Get the number of non-zero coefficients */

  int getNumElements() const ;

  /*! \brief Get the column (variable) lower bound vector */

  const double* getColLower() const ;

  /*! \brief Get the column (variable) upper bound vector */

  const double* getColUpper() const ;

  /*! \brief Return true if the variable is continuous */

  bool isContinuous(int) const ;

  /*! \brief Return true if the variable is binary */

  bool isBinary(int i) const ;

  /*! \brief Return true if the variable is general integer */

  bool isIntegerNonBinary (int i) const ;

  /*! \brief Return true if the variable is integer (general or binary) */

  bool isInteger (int i) const ;

  /*! \brief Get the row sense (constraint type) vector */

  const char* getRowSense() const ;

  /*! \brief Get the row (constraint) right-hand-side vector */

  const double* getRightHandSide() const ;

  /*! \brief Get the row (constraint) range vector */

  const double* getRowRange() const ;

  /*! \brief Get the row (constraint) lower bound vector */

  const double* getRowLower() const ;

  /*! \brief Get the row (constraint) upper bound vector */

  const double* getRowUpper() const ;

  /*! \brief Get the objective function coefficient vector */

  const double* getObjCoefficients() const ;

  /*! \brief Get the objective function sense (min/max) */

  double getObjSense() const ;

  /*! \brief Get a pointer to a row-major copy of the constraint matrix */

  const OsiPackedMatrix* getMatrixByRow() const ;

  /*! \brief Get a pointer to a column-major copy of the constraint matrix */

  const OsiPackedMatrix* getMatrixByCol() const ;
//@}

/*! \name Methods to modify the problem */
//@{

  /*! \brief Set a single variable to be continuous. */

  void setContinuous(int) ;

  /*! \brief Set a list of variables to be continuous. */

  void setContinuous(const int*, int) ;

  /*! \brief Set a single variable to be integer. */

  void setInteger(int) ;

  /*! \brief Set a list of variables to be integer. */

  void setInteger(const int*, int) ;

  /*! \brief Set the lower bound on a column (variable) */

  void setColLower(int, double) ;

  /*! \brief Set the upper bound on a column (variable) */

  void setColUpper(int, double) ;

  /*! \brief Set the lower bound on a row (constraint) */

  void setRowLower(int, double) ;

  /*! \brief Set the upper bound on a row (constraint) */

  void setRowUpper(int, double) ;

  /*! \brief Set the type of a row (constraint) */

  void setRowType(int, char, double, double) ;

  /*! \brief Set an objective function coefficient */

  void setObjCoeff (int, double) ;

  /*! \brief Set the sense (min/max) of the objective */

  void setObjSense(double) ;

  /*! \brief Add a column (variable) to the problem */

  void addCol(const OsiPackedVectorBase&,
	      const double, const double, const double) ;

  /*! \brief Remove column(s) (variable(s)) from the problem */

  void deleteCols(const int, const int*) ;

  /*! \brief Add a row (constraint) to the problem */

  void addRow(const OsiPackedVectorBase&, const double, const double) ;

  /*! \brief Add a row (constraint) to the problem */

  void addRow(const OsiPackedVectorBase&,
	      const char, const double, const double) ;

  /*! \brief Delete row(s) (constraint(s)) from the problem */

  void deleteRows(const int, const int*) ;

  /*! \brief Apply a row (constraint) cut (add one constraint) */

  void applyRowCut(const OsiRowCut&) ;

  /*! \brief Apply a column (variable) cut (adjust one or more bounds) */

  void applyColCut(const OsiColCut&) ;
//@}

/*! \name Solve methods */
//@{

  /*! \brief Solve an lp from scratch */

  void initialSolve() ;

  /*! \brief Build a warm start object for the current lp solution. */

  OsiWarmStart* getWarmStart() const ;

  /*! \brief Apply a warm start object. */

  bool setWarmStart(const OsiWarmStart*) ;

  /*! \brief Call dylp to reoptimize (warm or hot start). */

  void resolve() ;

  /*! \brief Set options to allow a hot start. */

  void markHotStart() ;

  /*! \brief Call dylp to reoptimize (hot start). */

  void solveFromHotStart() ;

  /*! \brief For dylp, a noop */

  void unmarkHotStart() ;

//@}

/*! \name Methods returning solver termination status */
//@{

  /*! \brief True if dylp abandoned the problem */

  bool isAbandoned() const ;

  /*! \brief True if dylp reported an optimal solution */

  bool isProvenOptimal() const ;

  /*! \brief True if dylp reported the problem to be primal infeasible */

  bool isProvenPrimalInfeasible() const ;

  /*! \brief True if dylp reported the problem to be dual infeasible (primal
	     unbounded)
  */
  bool isProvenDualInfeasible() const ;

  /*! \brief True if dylp reached the iteration limit */

  bool isIterationLimitReached() const ;

  /*! \brief Get the number of iterations for the last lp */

  int getIterationCount() const ;

//@}


/*! \name Methods to set/get solver parameters */
//@{

  /*! \brief Get dylp's value for infinity */

  double getInfinity() const ;

  /*! \brief Set an OSI integer parameter */

  bool setIntParam(OsiIntParam key, int value) ;

  /*! \brief Set an OSI double parameter */

  bool setDblParam(OsiDblParam key, double value) ;

  /*! \brief Set an OSI string parameter */

  bool setStrParam(OsiStrParam key, const std::string& value) ;

  /*! \brief Get an OSI integer parameter */

  bool getIntParam(OsiIntParam key, int& value) const ;

  /*! \brief Get an OSI double parameter */

  bool getDblParam(OsiDblParam key, double& value) const ;

  /*! \brief Get an OSI string parameter */

  bool getStrParam(OsiStrParam key, std::string& value) const ;

//@}

/*! \name Methods to obtain solution information */
//@{

  /*! \brief Return the vector of primal variables for the solution */

  const double* getColSolution() const ;

  /*! \brief Return the vector of dual variables for the solution */

  const double* getRowPrice() const ;

  /*! \brief Return the vector of reduced costs for the solution */

  const double* getReducedCost() const ;

  /*! \brief Return the vector of row activity for the solution */

  const double* getRowActivity() const ;

  /*! \brief Get the objective function value for the solution */

  double getObjValue() const ;

//@}

/*! \name Dylp-specific methods */
//@{

  /*! \brief Process an options (.spc) file */

  void dylp_controlfile(const char* name, const bool silent,
			const bool mustexist = true) ;

  /*! \brief Establish a log file */

  void dylp_logfile(const char* name, bool echo) ;

//@}

/*! \name Unsupported functions */
//@{


  /*! \brief Invoke the solver's built-in branch-and-bound algorithm. */

  void branchAndBound() ;

  /*! \brief Set the value of the primal variables in the problem solution */

  void setColSolution(const double*) ;

  /*! \brief Set the value of the dual variables in the problem solution */

  void setRowPrice(const double*) ;

  /*! \brief Is the primal objective limit reached? */

  bool isPrimalObjectiveLimitReached() const ;

  /*! \brief Is the dual objective limit reached? */

  bool isDualObjectiveLimitReached() const ;

  /*! \brief Get as many dual rays as the solver can provide */

  std::vector<double *> getDualRays(int) const ;

  /*! \brief Get as many primal rays as the solver can provide */

  std::vector<double *> getPrimalRays(int) const ;

  /*! \brief Write the problem into the specified file in MPS format. */

  void writeMps(const char *, const char* = "mps") const ;
//@}

/*
  These two structures contain dylp control parameters and tolerances. Leave
  them visible to the public for the nonce, until a better programmatic
  interface is available.
*/

  lpopts_struct* options ;
  lptols_struct* tolerances ;

private:

/*
  Private implementation state and helper functions. If you're contemplating
  using these, you should have a look at the code.
  See OsiDylpSolverInterface.ccp for descriptions.
*/ 
  consys_struct* consys ;
  lpprob_struct* lpprob ;
  lpstats_struct* statistics ;

/*! \name Dylp residual control variables */
//@{

  static int reference_count ;
  static bool yla05_ready ;
  static OsiDylpSolverInterface *dylp_owner ;

//@}


/*! \name Dylp instance control variables

  These variables track the state of individual ODSI instances.
*/
//@{

  /*! \brief Log file for this ODSI instance */

  ioid local_logchn ;

  /*! \brief Log echo control for this ODSI instance */

  bool local_gtxecho ;

  /*! \brief Result of last call to solver for this ODSI instance */

  lpret_enum lp_retval ;

  /*! \brief Objective function sense for this ODSI instance

    Coded 1.0 to minimize (default), -1.0 to maximize.
  */

  double obj_sense ;

//@}


/*! \name OSI parameters

  These variables support the implementation of the parameters defined in the
  OSI solver interface.
*/
//@{
/*!
  \var double osi_dual_tolerance
  \brief Value of the OsiDualTolerance parameter
*/
/*!
  \var double osi_primal_tolerance
  \brief Value of the OsiPrimalTolerance parameter
*/
/*!
  \var double osi_obj_offset
  \brief Value of the OsiObjOffset parameter
*/
/*!
  \var int osi_iterlim
  \brief Limit for total simplex iterations.
*/
/*!
  \var int osi_hotiterlim
  \brief Limit for total simplex iterations for an lp initiated from
	 solveFromHotStart.
*/
/*!
  \var string osi_probname
  \brief Problem name (from user or mps file).
*/

  double osi_dual_tolerance ;
  double osi_primal_tolerance ;
  double osi_obj_offset ;
  int osi_iterlim ;
  int osi_hotiterlim ;
  std::string osi_probname ;

//@}


/*! \name Cached problem information

  Problem information is cached for efficiency, to avoid repeated
  reconstruction of OSI structures from dylp structures.
*/
//@{

  mutable double* _col_x ;
  mutable double* _col_obj ;
  mutable double* _col_cbar ;
  mutable double* _row_lhs ;
  mutable double* _row_lower ;
  mutable double* _row_price ;
  mutable double* _row_range ;
  mutable double* _row_rhs ;
  mutable char* _row_sense ;
  mutable double* _row_upper ;
  mutable OsiPackedMatrix* _matrix_by_row ;
  mutable OsiPackedMatrix* _matrix_by_col ;

//@}

/*! \name Helper functions for problem construction */

//@{
  void construct_lpprob() ;
  void construct_options() ;
  void dylp_ioinit() ;
  void gen_rowparms(int rowcnt,
		    double *rhs, double *rhslow, contyp_enum *ctyp,
		    const double *rowlb, const double *rowub) ;
  void gen_rowparms(int rowcnt,
	       double *rhs, double *rhslow, contyp_enum *ctyp,
	       const char *sense, const double *rhsin, const double *range) ;
  void load_problem(const OsiPackedMatrix& matrix,
	 const double* col_lower, const double* col_upper, const double* obj,
	 const contyp_enum *ctyp, const double* rhs, const double* rhslow) ;
  void load_problem (const int colcnt, const int rowcnt,
	 const int *start, const int *index, const double *value,
	 const double* col_lower, const double* col_upper, const double* obj,
	 const contyp_enum *ctyp, const double* rhs, const double* rhslow) ;
//@}

/*! \name Destructor helpers */
//@{
  void destruct_col_cache() ;
  void destruct_row_cache() ;
  void destruct_cache() ;
  void destruct_problem() ;
  void detach_dylp() ;
//@}


/*! \name Helper functions for problem modification */
/*
  There are separate groups for member and static methods so that doxygen
  won't promote the group to the top level.
*/
//@{
  
  void add_col(const OsiPackedVectorBase& osi_coli,
    vartyp_enum vtypi,double vlbi, double vubi, double obji) ;
  void add_row(const OsiPackedVectorBase& osi_rowi, 
    char clazzi, contyp_enum ctypi, double rhsi, double rhslowi) ;

//@}

/*! \name Helper functions for problem modification */
//@{
  static contyp_enum bound_to_type(double lower, double upper) ;
  static contyp_enum sense_to_type(char type) ;
  static char type_to_sense(contyp_enum type) ;
  static void gen_rowiparms(contyp_enum* ctypi, double* rhsi, double* rhslowi, 
			    char sensei, double rhsini, double rangei) ;
  static void gen_rowiparms(contyp_enum* ctypi, double* rhsi, double* rhslowi, 
			    double rowlbi, double rowubi) ;
//@}

/*! \name Copy helpers

  Copy function templates for simple vectors and fixed-size objects, and
  specializations for various complex structures.
*/
//@{
  template<class T> static void copy(const T* src, T* dst, int n) ;
  template<class T> static T* copy(const T* src, int n) ;
  template<class T> static T* copy(const T* src) ;
/*
  Specializations for more complicated structures.
*/
  static basis_struct* copy_basis(const basis_struct* src) ;
  static void copy_basis(const basis_struct* src, basis_struct* dst) ;
  static consys_struct* copy_consys(const consys_struct* src) ;
  static lpprob_struct* copy_lpprob(const lpprob_struct* src) ;
  static void copy_consysmtx(const consys_struct* src, consys_struct* dst) ;
//@}

#ifndef _MSC_VER
/*! \name Copy verification functions

  Copy verification functions, to check that two structures are identical.
*/
//@{
  static void assert_same(const basis_struct& b1, 
  			  const basis_struct& b2, bool exact) ;
  static void assert_same(const consys_struct& c1, const 
    consys_struct& c2, bool exact) ;
  static void assert_same(const conbnd_struct& c1, const 
    conbnd_struct& c2, bool exact) ;
  static void assert_same(double d1, double d2, bool exact) ;
  static void assert_same(int i1, int i2, bool exact) ;
  static void assert_same(const lpprob_struct& l1, 
    const lpprob_struct& l2, bool exact) ;
  static void assert_same(const OsiDylpSolverInterface& o1, 
    const OsiDylpSolverInterface& o2, bool exact) ;
  template<class T> static void assert_same(const T& t1, 
    const T& t2, bool exact) ;
  template<class T> static void assert_same(const T* t1, const T* t2,
    int n, bool exact) ;
//@}
#endif

/*! \name Vector helper functions */
//@{
  static int idx(int i) ;
  template<class T> static T* inv_vec(T* data) ;
  static int inv(int i) ;

  static pkvec_struct* packed_vector(
    const OsiShallowPackedVector vector, int dimension) ;
//@}

/*! \name File i/o helper routines */
//@{
  static std::string make_filename(const std::string filename, 
			    const char *ext1, const char *ext2) ;
//@}

} ;



/*! \class OsiDylpWarmStartBasis
    \brief The dylp warm start object

  This derived class is necessary because dylp does not always work with the
  full constraint system. The warm start object needs to contain a list of
  the active constraints in addition to the status information included in
  OsiWarmStartBasis.

  Constraint status is coded using the OsiWarmStartBasis::Status codes. Active
  constraints are coded as atUpperBound or atLowerBound, inactive as isFree.
*/

#include "OsiWarmStartBasis.hpp"

class OsiDylpWarmStartBasis : public OsiWarmStartBasis

{ public:

/*! \name Functions to access constraint status */
//@{

  /*! \brief Return the number of active constraints */

  int getNumConstraint() const { return (numConstraints_) ; }

  /*! \brief Return the constraint status vector */

  const char *getConstraintStatus () const { return (constraintStatus_) ; }

  char *getConstraintStatus_nc () { return (constraintStatus_) ; }

  /*! \brief Return the status of a single constraint. */

  Status getConStatus (int i) const

  { const int st = (constraintStatus_[i>>2] >> ((i&3)<<1)) & 3 ;
    return (static_cast<OsiWarmStartBasis::Status>(st)) ; }

  /*! \brief Set the status of a single constraint */

  void setConStatus (int i, Status st)

  { char &st_byte = constraintStatus_[i>>2] ;
    st_byte &= ~(3 << ((i&3)<<1)) ;
    st_byte |= (st << ((i&3)<<1)) ; }

//@}

/*! \name Constructors, Destructors, and related functions */
//@{

  /*! \brief Default constructor (empty object) */

  OsiDylpWarmStartBasis ()

    : OsiWarmStartBasis(),
      numConstraints_(0),
      constraintStatus_(0)

  { /* intentionally left blank */ }


  /*! \brief Copy constructor */

  OsiDylpWarmStartBasis (const OsiDylpWarmStartBasis &ws)

    : OsiWarmStartBasis(ws),
      numConstraints_(ws.numConstraints_),
      constraintStatus_(0)

  { constraintStatus_ = new char[(numConstraints_+3)/4] ;
    CoinDisjointCopyN(ws.constraintStatus_,
		      (numConstraints_+3)/4,constraintStatus_) ; }


  /*! \brief Construct by copying status arrays */

  OsiDylpWarmStartBasis
  (int ns, int na, const char *sStat, const char *aStat, const char *cStat)

    : OsiWarmStartBasis(ns,na,sStat,aStat),
      numConstraints_(na),
      constraintStatus_(0)

  { constraintStatus_ = new char[(na+3)/4] ;
    CoinDisjointCopyN(cStat,(na+3)/4,constraintStatus_) ; }


  /*! \brief Destructor */

  ~OsiDylpWarmStartBasis ()

  { delete[] constraintStatus_ ; }


  /*! \brief Resize and initialize the warm start object.  */

  void setSize (int ns, int na)

  { OsiWarmStartBasis::setSize(ns,na) ;
    delete[] constraintStatus_ ;
    constraintStatus_ = new char[(na+3)/4] ;
    CoinFillN(constraintStatus_,(na+3)/4,(char)0) ;
    numConstraints_ = na ; }


  /*! \brief Assignment of status arrays (parameters destroyed) */

  void assignBasisStatus
  (int ns, int na, char *&sStat, char *&aStat, char *&cStat)

  { OsiWarmStartBasis::assignBasisStatus(ns,na,sStat,sStat) ;
    delete[] constraintStatus_ ;
    numConstraints_ = na ;
    constraintStatus_ = cStat ;
    cStat = 0 ; }


  /*! \brief Assignment of structure (parameter preserved) */

  OsiDylpWarmStartBasis& operator= (const OsiDylpWarmStartBasis &rhs)

  { if (this != &rhs)
    { OsiWarmStartBasis::operator= (rhs) ;
      delete[] constraintStatus_ ;
      constraintStatus_ = new char[(numConstraints_+3)/4] ;
      numConstraints_ = rhs.numConstraints_ ;
      CoinDisjointCopyN(rhs.constraintStatus_,(numConstraints_+3)/4,
			constraintStatus_) ; }
    
    return *this ; }

//@}

  private:

/*! \name Constraint status private data members */
//@{

  int numConstraints_ ;
  char *constraintStatus_ ;

//@}

} ;

/*
  OsiDylpSolverInterfaceTest.cpp
*/

void OsiDylpSolverInterfaceUnitTest(const std::string & mpsDir) ;

#endif // OsiDylpSolverInterface_H
#endif // COIN_USE_DYLP
