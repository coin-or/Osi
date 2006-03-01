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
#include "OsiAuxInfo.hpp"
#include "CoinHelperFunctions.hpp"

// Default Constructor
OsiAuxInfo::OsiAuxInfo(void * appData) 
  : appData_(appData)
{
}

// Destructor 
OsiAuxInfo::~OsiAuxInfo ()
{
}

// Clone
OsiAuxInfo *
OsiAuxInfo::clone() const
{
  return new OsiAuxInfo(*this);
}

// Copy constructor 
OsiAuxInfo::OsiAuxInfo(const OsiAuxInfo & rhs)
:
  appData_(rhs.appData_)
{
}
OsiAuxInfo &
OsiAuxInfo::operator=(const OsiAuxInfo &rhs)
{
  if (this != &rhs) {
    appData_ = rhs.appData_;
  }
  return *this;
}
// Default Constructor
OsiBabSolver::OsiBabSolver(int solverType) 
  :OsiAuxInfo(),
   solver_(NULL),
   solverType_(solverType),
   bestObjectiveValue_(1.0e100),
   bestSolution_(NULL),
   mipBound_(-1.0e100)
{
}

// Destructor 
OsiBabSolver::~OsiBabSolver ()
{
  delete [] bestSolution_;
}

// Clone
OsiAuxInfo *
OsiBabSolver::clone() const
{
  return new OsiBabSolver(*this);
}

// Copy constructor 
OsiBabSolver::OsiBabSolver(const OsiBabSolver & rhs)
:
  OsiAuxInfo(rhs),
  solver_(rhs.solver_),
  solverType_(rhs.solverType_),
  bestObjectiveValue_(rhs.bestObjectiveValue_),
  bestSolution_(NULL),
  mipBound_(rhs.mipBound_)
{
  if (rhs.bestSolution_) {
    assert (solver_);
    bestSolution_ = CoinCopyOfArray(rhs.bestSolution_,solver_->getNumCols());
  }
}
OsiBabSolver &
OsiBabSolver::operator=(const OsiBabSolver &rhs)
{
  if (this != &rhs) {
    OsiAuxInfo::operator=(rhs);
    delete [] bestSolution_;
    solver_ = rhs.solver_;
    solverType_ = rhs.solverType_;
    bestObjectiveValue_ = rhs.bestObjectiveValue_;
    bestSolution_ = NULL;
    mipBound_ = rhs.mipBound_;
    if (rhs.bestSolution_) {
      assert (solver_);
      bestSolution_ = CoinCopyOfArray(rhs.bestSolution_,solver_->getNumCols());
    }
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int
OsiBabSolver::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!solver_)
    return 0;
  if (bestObjectiveValue_<solutionValue&&bestSolution_) {
    int numberColumns = solver_->getNumCols();
    // new solution
    memcpy(betterSolution,bestSolution_,numberColumns*sizeof(double));
    solutionValue = bestObjectiveValue_;
    // free up
    delete [] bestSolution_;
    bestSolution_=NULL;
    bestObjectiveValue_=1.0e100;
    return 1;
  } else {
    return 0;
  }
}
// set solution
void
OsiBabSolver::setSolution(const OsiSolverInterface * solver)
{
  if (solver)
    solver_ = solver;
  assert (solver_);
  // just in case size has changed
  delete [] bestSolution_;
  bestSolution_ = CoinCopyOfArray(solver_->getColSolution(),solver_->getNumCols());
  bestObjectiveValue_ = solver_->getObjValue()*solver_->getObjSense();
}
// Get objective  (well mip bound)
double 
OsiBabSolver::mipBound() const
{
  assert (solver_);
  if (solverType_!=3)
    return solver_->getObjSense()*solver_->getObjValue();
  else
    return mipBound_;
}
// Returns true if node feasible
bool 
OsiBabSolver::mipFeasible() const
{
  assert (solver_);
  if (solverType_!=3)
    return solver_->isProvenOptimal();
  else
    return mipBound_<1.0e50;
}
