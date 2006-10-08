// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <string>
#include <cassert>
#include <cfloat>
#include "OsiSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CoinSort.hpp"
#include "OsiChooseVariable.hpp"
using namespace std;

OsiChooseVariable::OsiChooseVariable() :
  goodObjectiveValue_(COIN_DBL_MAX),
  upChange_(0.0),
  downChange_(0.0),
  goodSolution_(NULL),
  list_(NULL),
  useful_(NULL),
  solver_(NULL),
  status_(-1),
  bestObjectIndex_(-1),
  numberUnsatisfied_(0),
  numberStrong_(0),
  numberOnList_(0),
  numberStrongDone_(0),
  trustStrongForBound_(true),
  trustStrongForSolution_(true)
{
}

OsiChooseVariable::OsiChooseVariable(const OsiSolverInterface * solver) :
  goodObjectiveValue_(COIN_DBL_MAX),
  upChange_(0.0),
  downChange_(0.0),
  goodSolution_(NULL),
  solver_(solver),
  status_(-1),
  bestObjectIndex_(-1),
  numberUnsatisfied_(0),
  numberStrong_(0),
  numberOnList_(0),
  numberStrongDone_(0),
  trustStrongForBound_(true),
  trustStrongForSolution_(true)
{
  // create useful arrays
  int numberObjects = solver_->numberObjects();
  list_ = new int [numberObjects];
  useful_ = new double [numberObjects];
}

OsiChooseVariable::OsiChooseVariable(const OsiChooseVariable & rhs) 
{  
  goodObjectiveValue_ = rhs.goodObjectiveValue_;
  upChange_ = rhs.upChange_;
  downChange_ = rhs.downChange_;
  status_ = rhs.status_;
  bestObjectIndex_ = rhs.bestObjectIndex_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
  numberStrong_ = rhs.numberStrong_;
  numberOnList_ = rhs.numberOnList_;
  numberStrongDone_ = rhs.numberStrongDone_;
  trustStrongForBound_ = rhs.trustStrongForBound_;
  trustStrongForSolution_ = rhs.trustStrongForSolution_;
  solver_ = rhs.solver_;
  if (solver_) {
    int numberObjects = solver_->numberObjects();
    int numberColumns = solver_->getNumCols();
    if (rhs.goodSolution_) {
      goodSolution_ = CoinCopyOfArray(rhs.goodSolution_,numberColumns);
    } else {
      goodSolution_ = NULL;
    }
    list_ = CoinCopyOfArray(rhs.list_,numberObjects);
    useful_ = CoinCopyOfArray(rhs.useful_,numberObjects);
  } else {
    goodSolution_ = NULL;
    list_ = NULL;
    useful_ = NULL;
  }
}

OsiChooseVariable &
OsiChooseVariable::operator=(const OsiChooseVariable & rhs)
{
  if (this != &rhs) {
    delete [] goodSolution_;
    delete [] list_;
    delete [] useful_;
    goodObjectiveValue_ = rhs.goodObjectiveValue_;
    upChange_ = rhs.upChange_;
    downChange_ = rhs.downChange_;
    status_ = rhs.status_;
    bestObjectIndex_ = rhs.bestObjectIndex_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
    numberStrong_ = rhs.numberStrong_;
    numberOnList_ = rhs.numberOnList_;
    numberStrongDone_ = rhs.numberStrongDone_;
    trustStrongForBound_ = rhs.trustStrongForBound_;
    trustStrongForSolution_ = rhs.trustStrongForSolution_;
    solver_ = rhs.solver_;
    if (solver_) {
      int numberObjects = solver_->numberObjects();
      int numberColumns = solver_->getNumCols();
      if (rhs.goodSolution_) {
	goodSolution_ = CoinCopyOfArray(rhs.goodSolution_,numberColumns);
      } else {
	goodSolution_ = NULL;
      }
      list_ = CoinCopyOfArray(rhs.list_,numberObjects);
      useful_ = CoinCopyOfArray(rhs.useful_,numberObjects);
    } else {
      goodSolution_ = NULL;
      list_ = NULL;
      useful_ = NULL;
    }
  }
  return *this;
}


OsiChooseVariable::~OsiChooseVariable ()
{
  delete [] goodSolution_;
  delete [] list_;
  delete [] useful_;
}

// Clone
OsiChooseVariable *
OsiChooseVariable::clone() const
{
  return new OsiChooseVariable(*this);
}
// Set solver and redo arrays
void 
OsiChooseVariable::setSolver (const OsiSolverInterface * solver) 
{
  solver_ = solver;
  delete [] list_;
  delete [] useful_;
  // create useful arrays
  int numberObjects = solver_->numberObjects();
  list_ = new int [numberObjects];
  useful_ = new double [numberObjects];
}


// Initialize
int 
OsiChooseVariable::setupList ( OsiBranchingInformation *info, bool initialize)
{
  if (initialize) {
    status_=-2;
    delete [] goodSolution_;
    bestObjectIndex_=-1;
    numberStrongDone_=0;
    goodSolution_ = NULL;
    goodObjectiveValue_ = COIN_DBL_MAX;
  }
  numberOnList_=0;
  numberUnsatisfied_=0;
  int numberObjects = solver_->numberObjects();
  assert (numberObjects);
  double check = 0.0;
  int checkIndex=0;
  int bestPriority=INT_MAX;
  // pretend one strong even if none
  int maximumStrong= numberStrong_ ? numberStrong_ : 1;
  int putOther = numberObjects;
  int i;
  for (i=0;i<maximumStrong;i++) {
    list_[i]=-1;
    useful_[i]=0.0;
  }
  OsiObject ** object = info->solver_->objects();
  for (int i=0;i<numberObjects;i++) {
    int way;
    double value = object[i]->infeasibility(info,way);
    if (value>0.0) {
      numberUnsatisfied_++;
      int priorityLevel = object[i]->priority();
      // Better priority? Flush choices.
      if (priorityLevel<bestPriority) {
	for (int j=0;j<maximumStrong;j++) {
	  if (list_[j]>=0) {
	    int iObject = list_[j];
	    list_[j]=-1;
	    useful_[j]=0.0;
	    list_[--putOther]=iObject;
	  }
	}
	bestPriority = priorityLevel;
	check=0.0;
      } 
      if (priorityLevel==bestPriority) {
	if (value>check) {
	  //add to list
	  int iObject = list_[checkIndex];
	  if (iObject>=0)
	    list_[--putOther]=iObject;  // to end
	  list_[checkIndex]=i;
	  useful_[checkIndex]=value;
	  // find worst
	  check=COIN_DBL_MAX;
	  for (int j=0;j<maximumStrong;j++) {
	    if (list_[j]>=0) {
	      if (useful_[j]<check) {
		check=useful_[j];
		checkIndex=j;
	      }
	    } else {
	      check=0.0;
	      checkIndex = j;
	      break;
	    }
	  }
	} else {
	  // to end
	  list_[--putOther]=i;
	}
      }
    }
  }
  // Get list
  numberOnList_=0;
  for (i=0;i<maximumStrong;i++) {
    if (list_[i]>=0) {
      list_[numberOnList_]=list_[i];
      useful_[numberOnList_++]=-useful_[i];
    }
  }
  if (numberOnList_) {
    // Sort 
    CoinSort_2(useful_,useful_+numberOnList_,list_);
    // move others
    i = numberOnList_;
    for (;putOther<numberObjects;putOther++) 
      list_[i++]=list_[putOther];
    assert (i==numberUnsatisfied_);
    if (!numberStrong_)
      numberOnList_=0;
  } 
  return numberUnsatisfied_;
}
/* Choose a variable
   Returns - 
   -1 Node is infeasible
   0  Normal termination - we have a candidate
   1  All looks satisfied - no candidate
   2  We can change the bound on a variable - but we also have a strong branching candidate
   3  We can change the bound on a variable - but we have a non-strong branching candidate
   4  We can change the bound on a variable - no other candidates
   We can pick up branch from whichObject() and whichWay()
   We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
   If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
*/
int 
OsiChooseVariable::chooseVariable( OsiBranchingInformation *info)
{
  if (numberUnsatisfied_) {
    bestObjectIndex_=list_[0];
    return 0;
  } else {
    return 1;
  }
}
// Given a candidate (bestObjectIndex_) fill in useful information e.g. estimates
void 
OsiChooseVariable::updateInformation( OsiBranchingInformation *info)
{
  if (bestObjectIndex_>=0) {
    assert (bestObjectIndex_<solver_->numberObjects());
    OsiObject ** object = info->solver_->objects();
    upChange_ = object[bestObjectIndex_]->upEstimate();
    downChange_ = object[bestObjectIndex_]->downEstimate();
  }  
}
		
