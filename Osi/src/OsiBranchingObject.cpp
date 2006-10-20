// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiBranchingObject.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

// Default Constructor
OsiObject::OsiObject() 
  :infeasibility_(0.0),
   whichWay_(0),
   priority_(1000),
   numberWays_(2)
{
}


// Destructor 
OsiObject::~OsiObject ()
{
}

// Copy constructor 
OsiObject::OsiObject ( const OsiObject & rhs)
{
  infeasibility_ = rhs.infeasibility_;
  whichWay_ = rhs.whichWay_;
  priority_ = rhs.priority_;
  numberWays_ = rhs.numberWays_;
}

// Assignment operator 
OsiObject & 
OsiObject::operator=( const OsiObject& rhs)
{
  if (this!=&rhs) {
    infeasibility_ = rhs.infeasibility_;
    whichWay_ = rhs.whichWay_;
    priority_ = rhs.priority_;
    numberWays_ = rhs.numberWays_;
  }
  return *this;
}
// Return "up" estimate (default 1.0e-5)
double 
OsiObject::upEstimate() const
{
  return 1.0e-5;
}
// Return "down" estimate (default 1.0e-5)
double 
OsiObject::downEstimate() const
{
  return 1.0e-5;
}
// Column number if single column object -1 otherwise
int 
OsiObject::columnNumber() const
{
  return -1;
}
// Infeasibility - large is 0.5
double 
OsiObject::infeasibility(const OsiSolverInterface * solver, int & preferredWay) const
{
  OsiBranchingInformation info(solver);
  return infeasibility(&info,preferredWay);
}
// This does NOT set mutable stuff
double 
OsiObject::checkInfeasibility(const OsiBranchingInformation * info) const
{
  int way;
  double saveInfeasibility = infeasibility_;
  int saveWhichWay = whichWay_;
  double value = infeasibility(info,way);
  infeasibility_ = saveInfeasibility;
  whichWay_ = saveWhichWay;
  return value;
}
/* For the variable(s) referenced by the object,
   look at the current solution and set bounds to match the solution.
   Returns measure of how much it had to move solution to make feasible
*/
double 
OsiObject::feasibleRegion(OsiSolverInterface * solver) const 
{
  OsiBranchingInformation info(solver);
  return feasibleRegion(solver,&info);
}
// Default Constructor 
OsiBranchingObject::OsiBranchingObject()
{
  originalObject_=NULL;
  branchIndex_=0;
  value_=0.0;
  numberBranches_=2;
}

// Useful constructor
OsiBranchingObject::OsiBranchingObject (OsiSolverInterface * solver,
					 double value)
{
  originalObject_=NULL;
  branchIndex_=0;
  value_=value;
  numberBranches_=2;
}

// Copy constructor 
OsiBranchingObject::OsiBranchingObject ( const OsiBranchingObject & rhs)
{
  originalObject_=rhs.originalObject_;
  branchIndex_=rhs.branchIndex_;
  value_=rhs.value_;
  numberBranches_=rhs.numberBranches_;
}

// Assignment operator 
OsiBranchingObject & 
OsiBranchingObject::operator=( const OsiBranchingObject& rhs)
{
  if (this != &rhs) {
    originalObject_=rhs.originalObject_;
    branchIndex_=rhs.branchIndex_;
    value_=rhs.value_;
    numberBranches_=rhs.numberBranches_;
  }
  return *this;
}

// Destructor 
OsiBranchingObject::~OsiBranchingObject ()
{
}
// For debug
int 
OsiBranchingObject::columnNumber() const
{
  if (originalObject_)
    return originalObject_->columnNumber();
  else
    return -1;
}
/** Default Constructor

*/
OsiBranchingInformation::OsiBranchingInformation ()
  : objectiveValue_(COIN_DBL_MAX),
    cutoff_(COIN_DBL_MAX),
    direction_(COIN_DBL_MAX),
    integerTolerance_(1.0e-7),
    primalTolerance_(1.0e-7),
    timeRemaining_(COIN_DBL_MAX),
    solver_(NULL),
    lower_(NULL),
    solution_(NULL),
    upper_(NULL),
    hotstartSolution_(NULL),
    numberSolutions_(0),
    numberBranchingSolutions_(0),
    depth_(0)
{
}

/** Useful constructor
*/
OsiBranchingInformation::OsiBranchingInformation (const OsiSolverInterface * solver)
  : timeRemaining_(COIN_DBL_MAX),
    solver_(solver),
    hotstartSolution_(NULL),
    numberSolutions_(0),
    numberBranchingSolutions_(0),
    depth_(0)
{
  direction_ = solver_->getObjSense();
  objectiveValue_ = solver_->getObjValue();
  objectiveValue_ *= direction_;
  solver_->getDblParam(OsiDualObjectiveLimit,cutoff_) ;
  cutoff_ *= direction_;
  integerTolerance_ = solver_->getIntegerTolerance();
  solver_->getDblParam(OsiPrimalTolerance,primalTolerance_) ;
  lower_ = solver_->getColLower();
  solution_ = solver_->getColSolution();
  upper_ = solver_->getColUpper();
}
// Copy constructor 
OsiBranchingInformation::OsiBranchingInformation ( const OsiBranchingInformation & rhs)
{
  objectiveValue_ = rhs.objectiveValue_;
  cutoff_ = rhs.cutoff_;
  direction_ = rhs.direction_;
  integerTolerance_ = rhs.integerTolerance_;
  primalTolerance_ = rhs.primalTolerance_;
  timeRemaining_ = rhs.timeRemaining_;
  solver_ = rhs.solver_;
  lower_ = rhs.lower_;
  solution_ = rhs.solution_;
  upper_ = rhs.upper_;
  hotstartSolution_ = rhs.hotstartSolution_;
  numberSolutions_ = rhs.numberSolutions_;
  numberBranchingSolutions_ = rhs.numberBranchingSolutions_;
  depth_ = rhs.depth_;
}

// Clone
OsiBranchingInformation *
OsiBranchingInformation::clone() const
{
  return new OsiBranchingInformation(*this);
}

// Assignment operator 
OsiBranchingInformation & 
OsiBranchingInformation::operator=( const OsiBranchingInformation& rhs)
{
  if (this!=&rhs) {
    objectiveValue_ = rhs.objectiveValue_;
    cutoff_ = rhs.cutoff_;
    direction_ = rhs.direction_;
    integerTolerance_ = rhs.integerTolerance_;
    primalTolerance_ = rhs.primalTolerance_;
    timeRemaining_ = rhs.timeRemaining_;
    lower_ = rhs.lower_;
    solution_ = rhs.solution_;
    upper_ = rhs.upper_;
    hotstartSolution_ = rhs.hotstartSolution_;
    numberSolutions_ = rhs.numberSolutions_;
    numberBranchingSolutions_ = rhs.numberBranchingSolutions_;
    depth_ = rhs.depth_;
  }
  return *this;
}

// Destructor 
OsiBranchingInformation::~OsiBranchingInformation ()
{
}
// Default Constructor 
OsiTwoWayBranchingObject::OsiTwoWayBranchingObject()
  :OsiBranchingObject()
{
  firstBranch_=0;
}

// Useful constructor
OsiTwoWayBranchingObject::OsiTwoWayBranchingObject (OsiSolverInterface * solver, 
						      const OsiObject * object,
						      int way , double value)
  :OsiBranchingObject(solver,value)
{
  originalObject_ = object;
  firstBranch_=way;
}
  

// Copy constructor 
OsiTwoWayBranchingObject::OsiTwoWayBranchingObject ( const OsiTwoWayBranchingObject & rhs) :OsiBranchingObject(rhs)
{
  firstBranch_=rhs.firstBranch_;
}

// Assignment operator 
OsiTwoWayBranchingObject & 
OsiTwoWayBranchingObject::operator=( const OsiTwoWayBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiBranchingObject::operator=(rhs);
    firstBranch_=rhs.firstBranch_;
  }
  return *this;
}

// Destructor 
OsiTwoWayBranchingObject::~OsiTwoWayBranchingObject ()
{
}

/********* Simple Integers *******************************/
/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
OsiSimpleInteger::OsiSimpleInteger ()
  : OsiObject(),
    columnNumber_(-1),
    originalLower_(0.0),
    originalUpper_(1.0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
OsiSimpleInteger::OsiSimpleInteger (const OsiSolverInterface * solver, int iColumn)
  : OsiObject()
{
  columnNumber_ = iColumn ;
  originalLower_ = solver->getColLower()[columnNumber_] ;
  originalUpper_ = solver->getColUpper()[columnNumber_] ;
}

  
// Useful constructor - passed solver index and original bounds
OsiSimpleInteger::OsiSimpleInteger ( int iColumn, double lower, double upper)
  : OsiObject()
{
  columnNumber_ = iColumn ;
  originalLower_ = lower;
  originalUpper_ = upper;
}

// Copy constructor 
OsiSimpleInteger::OsiSimpleInteger ( const OsiSimpleInteger & rhs)
  :OsiObject(rhs)

{
  columnNumber_ = rhs.columnNumber_;
  originalLower_ = rhs.originalLower_;
  originalUpper_ = rhs.originalUpper_;
}

// Clone
OsiObject *
OsiSimpleInteger::clone() const
{
  return new OsiSimpleInteger(*this);
}

// Assignment operator 
OsiSimpleInteger & 
OsiSimpleInteger::operator=( const OsiSimpleInteger& rhs)
{
  if (this!=&rhs) {
    OsiObject::operator=(rhs);
    columnNumber_ = rhs.columnNumber_;
    originalLower_ = rhs.originalLower_;
    originalUpper_ = rhs.originalUpper_;
  }
  return *this;
}

// Destructor 
OsiSimpleInteger::~OsiSimpleInteger ()
{
}
/* Reset variable bounds to their original values.
   
Bounds may be tightened, so it may be good to be able to reset them to
their original values.
*/
void 
OsiSimpleInteger::resetBounds(const OsiSolverInterface * solver) 
{
  originalLower_ = solver->getColLower()[columnNumber_] ;
  originalUpper_ = solver->getColUpper()[columnNumber_] ;
}
// Redoes data when sequence numbers change
void 
OsiSimpleInteger::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
  int i;
  for (i=0;i<numberColumns;i++) {
    if (originalColumns[i]==columnNumber_)
      break;
  }
  if (i<numberColumns)
    columnNumber_=i;
  else
    abort(); // should never happen
}

// Infeasibility - large is 0.5
double 
OsiSimpleInteger::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  double nearest = floor(value+(1.0-0.5));
  if (nearest>value) { 
    whichWay=1;
  } else {
    whichWay=0;
  }
  infeasibility_ = fabs(value-nearest);
  whichWay_=whichWay;
  if (infeasibility_<=info->integerTolerance_) 
    return 0.0;
  else
    return infeasibility_;
}

// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
double
OsiSimpleInteger::feasibleRegion(OsiSolverInterface * solver,
				 const OsiBranchingInformation * info) const
{
  double value = info->solution_[columnNumber_];
  double newValue = CoinMax(value, info->lower_[columnNumber_]);
  newValue = CoinMin(newValue, info->upper_[columnNumber_]);
  newValue = floor(newValue+0.5);
  solver->setColLower(columnNumber_,newValue);
  solver->setColUpper(columnNumber_,newValue);
  return fabs(value-newValue);
}
/* Column number if single column object -1 otherwise,
   so returns >= 0
   Used by heuristics
*/
int 
OsiSimpleInteger::columnNumber() const
{
  return columnNumber_;
}
// Creates a branching object
OsiBranchingObject * 
OsiSimpleInteger::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const 
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  assert (fabs(value-nearest)>info->integerTolerance_);
#endif
  OsiBranchingObject * branch = new OsiIntegerBranchingObject(solver,this,way,
					     value);
  return branch;
}
// Return "down" estimate
double 
OsiSimpleInteger::downEstimate() const
{
  if (whichWay_)
    return 1.0-infeasibility_;
  else
    return infeasibility_;
}
// Return "up" estimate
double 
OsiSimpleInteger::upEstimate() const
{
  if (!whichWay_)
    return 1.0-infeasibility_;
  else
    return infeasibility_;
}

// Default Constructor 
OsiIntegerBranchingObject::OsiIntegerBranchingObject()
  :OsiTwoWayBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
}

// Useful constructor
OsiIntegerBranchingObject::OsiIntegerBranchingObject (OsiSolverInterface * solver, 
						      const OsiSimpleInteger * object,
						      int way , double value)
  :OsiTwoWayBranchingObject(solver,object, way, value)
{
  int iColumn = object->columnNumber();
  down_[0] = solver->getColLower()[iColumn];
  down_[1] = floor(value_);
  up_[0] = ceil(value_);
  up_[1] = solver->getColUpper()[iColumn];
}
  

// Copy constructor 
OsiIntegerBranchingObject::OsiIntegerBranchingObject ( const OsiIntegerBranchingObject & rhs) :OsiTwoWayBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
}

// Assignment operator 
OsiIntegerBranchingObject & 
OsiIntegerBranchingObject::operator=( const OsiIntegerBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiTwoWayBranchingObject::operator=(rhs);
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
  }
  return *this;
}
OsiBranchingObject * 
OsiIntegerBranchingObject::clone() const
{ 
  return (new OsiIntegerBranchingObject(*this));
}


// Destructor 
OsiIntegerBranchingObject::~OsiIntegerBranchingObject ()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of branchIndex_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
  Returns change in guessed objective on next branch
*/
double
OsiIntegerBranchingObject::branch(OsiSolverInterface * solver)
{
  const OsiSimpleInteger * obj =
    dynamic_cast <const OsiSimpleInteger *>(originalObject_) ;
  assert (obj);
  int iColumn = obj->columnNumber();
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  if (way<0) {
#ifdef OSI_DEBUG
  { double olb,oub ;
    olb = solver->getColLower()[iColumn] ;
    oub = solver->getColUpper()[iColumn] ;
    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
#endif
    solver->setColLower(iColumn,down_[0]);
    solver->setColUpper(iColumn,down_[1]);
  } else {
#ifdef OSI_DEBUG
  { double olb,oub ;
    olb = solver->getColLower()[iColumn] ;
    oub = solver->getColUpper()[iColumn] ;
    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
#endif
    solver->setColLower(iColumn,up_[0]);
    solver->setColUpper(iColumn,up_[1]);
  }
  branchIndex_++;
  return 0.0;
}
// Print what would happen  
void
OsiIntegerBranchingObject::print(const OsiSolverInterface * solver)
{
  const OsiSimpleInteger * obj =
    dynamic_cast <const OsiSimpleInteger *>(originalObject_) ;
  assert (obj);
  int iColumn = obj->columnNumber();
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  if (way<0) {
  { double olb,oub ;
    olb = solver->getColLower()[iColumn] ;
    oub = solver->getColUpper()[iColumn] ;
    printf("OsiInteger would branch down on var %d : [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
  } else {
  { double olb,oub ;
    olb = solver->getColLower()[iColumn] ;
    oub = solver->getColUpper()[iColumn] ;
    printf("OsiInteger would branch up on var %d : [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
  }
}
// Default Constructor 
OsiSOS::OsiSOS ()
  : OsiObject(),
    members_(NULL),
    weights_(NULL),
    numberMembers_(0),
    sosType_(-1),
    integerValued_(false)
{
}

// Useful constructor (which are indices)
OsiSOS::OsiSOS (const OsiSolverInterface * solver,  int numberMembers,
	   const int * which, const double * weights, int type)
  : numberMembers_(numberMembers),
    sosType_(type)
{
  integerValued_ = type==1; // not strictly true - should check problem
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    weights_ = new double[numberMembers_];
    memcpy(members_,which,numberMembers_*sizeof(int));
    if (weights) {
      memcpy(weights_,weights,numberMembers_*sizeof(double));
    } else {
      for (int i=0;i<numberMembers_;i++)
        weights_[i]=i;
    }
    // sort so weights increasing
    CoinSort_2(weights_,weights_+numberMembers_,members_);
    double last = -COIN_DBL_MAX;
    int i;
    for (i=0;i<numberMembers_;i++) {
      double possible = CoinMax(last+1.0e-10,weights_[i]);
      weights_[i] = possible;
      last=possible;
    }
  } else {
    members_ = NULL;
    weights_ = NULL;
  }
  assert (sosType_>0&&sosType_<3);
}

// Copy constructor 
OsiSOS::OsiSOS ( const OsiSOS & rhs)
  :OsiObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  sosType_ = rhs.sosType_;
  integerValued_ = rhs.integerValued_;
  if (numberMembers_) {
    members_ = new int[numberMembers_];
    weights_ = new double[numberMembers_];
    memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
    memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
  } else {
    members_ = NULL;
    weights_ = NULL;
  }
}

// Clone
OsiObject *
OsiSOS::clone() const
{
  return new OsiSOS(*this);
}

// Assignment operator 
OsiSOS & 
OsiSOS::operator=( const OsiSOS& rhs)
{
  if (this!=&rhs) {
    OsiObject::operator=(rhs);
    delete [] members_;
    delete [] weights_;
    numberMembers_ = rhs.numberMembers_;
    sosType_ = rhs.sosType_;
    integerValued_ = rhs.integerValued_;
    if (numberMembers_) {
      members_ = new int[numberMembers_];
      weights_ = new double[numberMembers_];
      memcpy(members_,rhs.members_,numberMembers_*sizeof(int));
      memcpy(weights_,rhs.weights_,numberMembers_*sizeof(double));
    } else {
      members_ = NULL;
      weights_ = NULL;
    }
  }
  return *this;
}

// Destructor 
OsiSOS::~OsiSOS ()
{
  delete [] members_;
  delete [] weights_;
}

// Infeasibility - large is 0.5
double 
OsiSOS::infeasibility(const OsiBranchingInformation * info,int & whichWay) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  const OsiSolverInterface * solver = info->solver_;
  const double * solution = solver->getColSolution();
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  //double largestValue=0.0;
  double integerTolerance = info->integerTolerance_;
  double weight = 0.0;
  double sum =0.0;

  // check bounds etc
  double lastWeight=-1.0e100;
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    if (lastWeight>=weights_[j]-1.0e-7)
      throw CoinError("Weights too close together in SOS","infeasibility","OsiSOS");
    double value = CoinMax(0.0,solution[iColumn]);
    sum += value;
    if (value>integerTolerance&&upper[iColumn]) {
      // Possibly due to scaling a fixed variable might slip through
      if (value>upper[iColumn]) {
        value=upper[iColumn];
#ifdef COIN_DEVELOP
	printf("** Variable %d (%d) has value %g and upper bound of %g\n",
	       iColumn,j,value,upper[iColumn]);
#endif
      } 
      weight += weights_[j]*value;
      if (firstNonZero<0)
        firstNonZero=j;
      lastNonZero=j;
    }
  }
  whichWay=1;
  if (lastNonZero-firstNonZero>=sosType_) {
    // find where to branch
    assert (sum>0.0);
    weight /= sum;
    //int iWhere;
    //for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    //if (weight<weights_[iWhere+1])
    //break;
    // probably best to use pseudo duals
    double value = lastNonZero-firstNonZero+1;
    value *= 0.5/((double) numberMembers_);
    return value;
  } else {
    return 0.0; // satisfied
  }
}

// This looks at solution and sets bounds to contain solution
double
OsiSOS::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  const double * solution = info->solution_;
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double sum =0.0;
  // Find largest one or pair
  double movement=0.0;
  if (sosType_==1) {
    for (j=0;j<numberMembers_;j++) {
      int iColumn = members_[j];
      double value = CoinMax(0.0,solution[iColumn]);
      if (value>sum&&upper[iColumn]) {
	firstNonZero=j;
	sum=value;
      }
    }
    lastNonZero=firstNonZero;
  } else {
    // type 2
    for (j=1;j<numberMembers_;j++) {
      int iColumn = members_[j];
      int jColumn = members_[j-1];
      double value1 = CoinMax(0.0,solution[iColumn]);
      double value0 = CoinMax(0.0,solution[jColumn]);
      double value = value0+value1;
      if (value>sum) {
	if (upper[iColumn]||upper[jColumn]) {
	  firstNonZero=upper[jColumn] ? j-1 : j;
	  lastNonZero=upper[iColumn] ? j : j-1;
	  sum=value;
	}
      }
    }
  }
  for (j=0;j<numberMembers_;j++) {
    if (j<firstNonZero||j>lastNonZero) {
      int iColumn = members_[j];
      double value = CoinMax(0.0,solution[iColumn]);
      movement += value;
      solver->setColUpper(iColumn,0.0);
    }
  }
  return movement;
}
// Redoes data when sequence numbers change
void 
OsiSOS::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
  int n2=0;
  for (int j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    int i;
    for (i=0;i<numberColumns;i++) {
      if (originalColumns[i]==iColumn)
        break;
    }
    if (i<numberColumns) {
      members_[n2]=i;
      weights_[n2++]=weights_[j];
    }
  }
  if (n2<numberMembers_) {
    printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2);
    numberMembers_=n2;
  }
}
// Return "up" estimate (default 1.0e-5)
double 
OsiSOS::upEstimate() const
{
  return 1.0e-5;
}
// Return "down" estimate (default 1.0e-5)
double 
OsiSOS::downEstimate() const
{
  return 1.0e-5;
}


// Creates a branching object
OsiBranchingObject * 
OsiSOS::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
  int j;
  const double * solution = info->solution_;
  double tolerance = info->primalTolerance_;
  const double * upper = info->upper_;
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  for (j=0;j<numberMembers_;j++) {
    int iColumn = members_[j];
    if (upper[iColumn]) {
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (firstNonFixed<0)
	firstNonFixed=j;
      lastNonFixed=j;
      if (value>tolerance) {
	weight += weights_[j]*value;
	if (firstNonZero<0)
	  firstNonZero=j;
	lastNonZero=j;
      }
    }
  }
  assert (lastNonZero-firstNonZero>=sosType_) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  double separator=0.0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  if (sosType_==1) {
    // SOS 1
    separator = 0.5 *(weights_[iWhere]+weights_[iWhere+1]);
  } else {
    // SOS 2
    if (iWhere==firstNonFixed)
      iWhere++;;
    if (iWhere==lastNonFixed-1)
      iWhere = lastNonFixed-2;
    separator = weights_[iWhere+1];
  }
  // create object
  OsiBranchingObject * branch;
  branch = new OsiSOSBranchingObject(solver,this,way,separator);
  return branch;
}
// Default Constructor 
OsiSOSBranchingObject::OsiSOSBranchingObject()
  :OsiTwoWayBranchingObject()
{
}

// Useful constructor
OsiSOSBranchingObject::OsiSOSBranchingObject (OsiSolverInterface * solver,
					      const OsiSOS * set,
					      int way ,
					      double separator)
  :OsiTwoWayBranchingObject(solver, set,way,separator)
{
}

// Copy constructor 
OsiSOSBranchingObject::OsiSOSBranchingObject ( const OsiSOSBranchingObject & rhs) :OsiTwoWayBranchingObject(rhs)
{
}

// Assignment operator 
OsiSOSBranchingObject & 
OsiSOSBranchingObject::operator=( const OsiSOSBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiTwoWayBranchingObject::operator=(rhs);
  }
  return *this;
}
OsiBranchingObject * 
OsiSOSBranchingObject::clone() const
{ 
  return (new OsiSOSBranchingObject(*this));
}


// Destructor 
OsiSOSBranchingObject::~OsiSOSBranchingObject ()
{
}
double
OsiSOSBranchingObject::branch(OsiSolverInterface * solver)
{
  const OsiSOS * set =
    dynamic_cast <const OsiSOS *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  branchIndex_++;
  int numberMembers = set->numberMembers();
  const int * which = set->members();
  const double * weights = set->weights();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  // *** for way - up means fix all those in down section
  if (way<0) {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > value_)
	break;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) 
      solver->setColUpper(which[i],0.0);
  } else {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] >= value_)
	break;
      else
	solver->setColUpper(which[i],0.0);
    }
    assert (i<numberMembers);
  }
  return 0.0;
}
// Print what would happen  
void
OsiSOSBranchingObject::print(const OsiSolverInterface * solver)
{
  const OsiSOS * set =
    dynamic_cast <const OsiSOS *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  int numberMembers = set->numberMembers();
  const int * which = set->members();
  const double * weights = set->weights();
  //const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int first=numberMembers;
  int last=-1;
  int numberFixed=0;
  int numberOther=0;
  int i;
  for ( i=0;i<numberMembers;i++) {
    double bound = upper[which[i]];
    if (bound) {
      first = CoinMin(first,i);
      last = CoinMax(last,i);
    }
  }
  // *** for way - up means fix all those in down section
  if (way<0) {
    printf("SOS Down");
    for ( i=0;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (weights[i] > value_)
	break;
      else if (bound)
	numberOther++;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (bound)
	numberFixed++;
    }
  } else {
    printf("SOS Up");
    for ( i=0;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (weights[i] >= value_)
	break;
      else if (bound)
	numberFixed++;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      double bound = upper[which[i]];
      if (bound)
	numberOther++;
    }
  }
  printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
	 value_,which[first],weights[first],which[last],weights[last],numberFixed,numberOther);
}
  
