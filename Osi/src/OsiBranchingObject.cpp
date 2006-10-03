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
  return infeasibility(solver,&info,preferredWay);
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
  solver_=NULL;
  originalObject_=NULL;
  branchIndex_=0;
  value_=0.0;
  numberBranches_=2;
}

// Useful constructor
OsiBranchingObject::OsiBranchingObject (OsiSolverInterface * solver,
					 double value)
{
  solver_= solver;
  originalObject_=NULL;
  branchIndex_=0;
  value_=value;
  numberBranches_=2;
}

// Copy constructor 
OsiBranchingObject::OsiBranchingObject ( const OsiBranchingObject & rhs)
{
  solver_=rhs.solver_;
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
    solver_=rhs.solver_;
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
/** Default Constructor

*/
OsiBranchingInformation::OsiBranchingInformation ()
  : objectiveValue_(COIN_DBL_MAX),
    cutoff_(COIN_DBL_MAX),
    direction_(COIN_DBL_MAX),
    integerTolerance_(1.0e-7),
    primalTolerance_(1.0e-7),
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
  :  hotstartSolution_(NULL),
     numberSolutions_(0),
     numberBranchingSolutions_(0),
     depth_(0)
{
  direction_ = solver->getObjSense();
  objectiveValue_ = solver->getObjValue();
  objectiveValue_ *= direction_;
  solver->getDblParam(OsiDualObjectiveLimit,cutoff_) ;
  cutoff_ *= direction_;
  integerTolerance_ = solver->getIntegerTolerance();
  solver->getDblParam(OsiPrimalTolerance,primalTolerance_) ;
  lower_ = solver->getColLower();
  solution_ = solver->getColSolution();
  upper_ = solver->getColUpper();
}
// Copy constructor 
OsiBranchingInformation::OsiBranchingInformation ( const OsiBranchingInformation & rhs)
{
  objectiveValue_ = rhs.objectiveValue_;
  cutoff_ = rhs.cutoff_;
  direction_ = rhs.direction_;
  integerTolerance_ = rhs.integerTolerance_;
  primalTolerance_ = rhs.primalTolerance_;
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
OsiSimpleInteger::infeasibility(const OsiSolverInterface * solver, 
				const OsiBranchingInformation * info, int & preferredWay) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  double nearest = floor(value+(1.0-0.5));
  if (nearest>value) 
    preferredWay=1;
  else
    preferredWay=-1;
  double weight = fabs(value-nearest);
  if (fabs(value-nearest)<=info->integerTolerance_) 
    return 0.0;
  else
    return weight;
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
OsiSimpleInteger::createBranch(OsiSolverInterface * solver, int way) const 
{
  const double * solution = solver->getColSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  assert (upper[columnNumber_]>lower[columnNumber_]);
#ifndef NDEBUG
  double nearest = floor(value+0.5);
  double integerTolerance = solver->getIntegerTolerance();
  assert (fabs(value-nearest)>integerTolerance);
#endif
  OsiBranchingObject * branch = new OsiIntegerBranchingObject(solver,this,way,
					     value);
  return branch;
}

// Default Constructor 
OsiIntegerBranchingObject::OsiIntegerBranchingObject()
  :OsiBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
  firstBranch_=0;
}

// Useful constructor
OsiIntegerBranchingObject::OsiIntegerBranchingObject (OsiSolverInterface * solver, 
						      const OsiSimpleInteger * object,
						      int way , double value)
  :OsiBranchingObject(solver,value)
{
  originalObject_ = object;
  int iColumn = object->columnNumber();
  down_[0] = solver_->getColLower()[iColumn];
  down_[1] = floor(value_);
  up_[0] = ceil(value_);
  up_[1] = solver_->getColUpper()[iColumn];
  firstBranch_=way;
}
  

// Copy constructor 
OsiIntegerBranchingObject::OsiIntegerBranchingObject ( const OsiIntegerBranchingObject & rhs) :OsiBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
  firstBranch_=rhs.firstBranch_;
}

// Assignment operator 
OsiIntegerBranchingObject & 
OsiIntegerBranchingObject::operator=( const OsiIntegerBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiBranchingObject::operator=(rhs);
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
    firstBranch_=rhs.firstBranch_;
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
OsiIntegerBranchingObject::branch(bool normalBranch)
{
  const OsiSimpleInteger * obj =
    dynamic_cast <const OsiSimpleInteger *>(originalObject_) ;
  assert (obj);
  int iColumn = obj->columnNumber();
  int way = (!branchIndex_) ? firstBranch_ : -firstBranch_;
  if (way<0) {
#ifdef OSI_DEBUG
  { double olb,oub ;
    olb = solver_->getColLower()[iColumn] ;
    oub = solver_->getColUpper()[iColumn] ;
    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
#endif
    solver_->setColLower(iColumn,down_[0]);
    solver_->setColUpper(iColumn,down_[1]);
  } else {
#ifdef OSI_DEBUG
  { double olb,oub ;
    olb = solver_->getColLower()[iColumn] ;
    oub = solver_->getColUpper()[iColumn] ;
    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
#endif
    solver_->setColLower(iColumn,up_[0]);
    solver_->setColUpper(iColumn,up_[1]);
  }
  branchIndex_++;
  return 0.0;
}
  
