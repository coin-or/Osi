// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "ClpPresolve.hpp"
static double totalTime=0.0;

//#############################################################################
// Solve methods
//#############################################################################
void OsiClpSolverInterface::initialSolve()
{
  ClpSimplex solver;
  double time1 = CoinCpuTime();
  solver.borrowModel(*modelPtr_);
  // Treat as if user simplex not enabled
  int saveSolveType=solver.solveType();
  bool doingPrimal = solver.algorithm()>0;
  if (saveSolveType==2) {
    disableSimplexInterface();
    solver.setSolveType(1);
  }
  // Set message handler to have same levels etc
  solver.passInMessageHandler(handler_);
  // set reasonable defaults
  bool takeHint;
  OsiHintStrength strength;
  // Switch off printing if asked to
  bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  int saveMessageLevel=messageHandler()->logLevel();
  if (strength!=OsiHintIgnore&&takeHint) {
    if (saveMessageLevel)
      solver.messageHandler()->setLogLevel(saveMessageLevel-1);
  }
  // scaling
  // save initial state
  const double * rowScale1 = solver.rowScale();
  if (modelPtr_->solveType()==1) {
    gotHint = (getHintParam(OsiDoScale,takeHint,strength));
    assert (gotHint);
    if (strength==OsiHintIgnore||takeHint)
      solver.scaling(1);
    else
      solver.scaling(0);
  } else {
    solver.scaling(0);
  }
  //solver.setDualBound(1.0e6);
  //solver.setDualTolerance(1.0e-7);

  //ClpDualRowSteepest steep;
  //solver.setDualRowPivotAlgorithm(steep);
  //solver.setPrimalTolerance(1.0e-8);
  //ClpPrimalColumnSteepest steepP;
  //solver.setPrimalColumnPivotAlgorithm(steepP);
  /*
    If basis then do primal (as user could do dual with resolve)
    If not then see if dual feasible (and allow for gubs etc?)
   */
  bool doPrimal = (basis_.numberBasicStructurals()>0);
  setBasis(basis_,&solver);

  // sort out hints;
  // algorithm 0 whatever, -1 force dual, +1 force primal
  int algorithm = 0;
  gotHint = (getHintParam(OsiDoDualInInitial,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
    algorithm = takeHint ? -1 : 1;
  // crash 0 do lightweight if all slack, 1 do, -1 don't
  int doCrash=0;
  gotHint = (getHintParam(OsiDoCrash,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
    doCrash = takeHint ? 1 : -1;
	 
  // presolve
  gotHint = (getHintParam(OsiDoPresolveInInitial,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    ClpPresolve pinfo;
    ClpSimplex * model2 = pinfo.presolvedModel(solver,1.0e-8);
    if (!model2) {
      // problem found to be infeasible - whats best?
      model2 = &solver;
    }
      
    // change from 200 (unless changed)
    if (model2->factorization()->maximumPivots()==200)
      model2->factorization()->maximumPivots(100+model2->numberRows()/50);
    int savePerturbation = model2->perturbation();
    if (savePerturbation==100)
      model2->setPerturbation(50);
    if (!doPrimal) {
      // faster if bounds tightened
      //int numberInfeasibilities = model2->tightenPrimalBounds();
      model2->tightenPrimalBounds();
      // look further
      bool crashResult=false;
      if (doCrash>0)
	crashResult =  (solver.crash(1000.0,1)>0);
      else if (doCrash==0&&algorithm>0)
	crashResult =  (solver.crash(1000.0,1)>0);
      doPrimal=crashResult;
    }
    if (algorithm<0)
      doPrimal=false;
    else if (algorithm>0)
      doPrimal=true;
    if (!doPrimal) {
      //if (numberInfeasibilities)
      //std::cout<<"** Analysis indicates model infeasible"
      //       <<std::endl;
      // up dual bound for safety
      //model2->setDualBound(1.0e11);
      model2->dual();
      // check if clp thought it was in a loop
      if (model2->status()==3&&
	  model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->primal();
      }
    } else {
      // up infeasibility cost for safety
      //model2->setInfeasibilityCost(1.0e10);
      model2->primal();
      // check if clp thought it was in a loop
      if (model2->status()==3
	  &&model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->dual();
      }
    }
    model2->setPerturbation(savePerturbation);
    if (model2!=&solver) {
      pinfo.postsolve(true);
    
      delete model2;
      //printf("Resolving from postsolved model\n");
      // later try without (1) and check duals before solve
      solver.primal(1);
    }
    lastAlgorithm_=1; // primal
    //if (solver.numberIterations())
    //printf("****** iterated %d\n",solver.numberIterations());
  } else {
    if (doPrimal) {
      if (doCrash>0)
	solver.crash(1000.0,2);
      else if (doCrash==0)
	solver.crash(1000.0,0);
    }
    if (algorithm<0)
      doPrimal=false;
    else if (algorithm>0)
      doPrimal=true;
    if (!doPrimal) {
      //printf("doing dual\n");
      solver.dual();
      lastAlgorithm_=2; // dual
      // check if clp thought it was in a loop
      if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
	// switch algorithm
	solver.primal();
	lastAlgorithm_=1; // primal
      }
    } else {
      //printf("doing primal\n");
      solver.primal();
      lastAlgorithm_=1; // primal
      // check if clp thought it was in a loop
      if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
	// switch algorithm
	solver.dual();
	lastAlgorithm_=2; // dual
      }
    }
  }
  basis_ = getBasis(&solver);
  //basis_.print();
  solver.messageHandler()->setLogLevel(saveMessageLevel);
  const double * rowScale2 = solver.rowScale();
  if (!rowScale1&&rowScale2) {
    // need to release memory
    solver.setRowScale(NULL);
    solver.setColumnScale(NULL);
  }
  solver.returnModel(*modelPtr_);
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  time1 = CoinCpuTime()-time1;
  totalTime += time1;
  //std::cout<<time1<<" seconds - total "<<totalTime<<std::endl;
}

//-----------------------------------------------------------------------------
void OsiClpSolverInterface::resolve()
{
  //#define BORROW_MODEL
#ifdef BORROW_MODEL
  ClpSimplex solver;
  solver.borrowModel(*modelPtr_);
  int saveSolveType=solver.solveType();
  bool doingPrimal = solver.algorithm()>0;
  if (saveSolveType==2) {
    disableSimplexInterface();
    solver.setSolveType(1);
  }
  // Set message handler to have same levels etc
  solver.passInMessageHandler(handler_);
  //basis_.print();
  setBasis(basis_,&solver);
  // set reasonable defaults
  bool takeHint;
  OsiHintStrength strength;
  // Switch off printing if asked to
  bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  int saveMessageLevel=messageHandler()->logLevel();
  if (strength!=OsiHintIgnore&&takeHint) {
    if (saveMessageLevel)
      solver.messageHandler()->setLogLevel(saveMessageLevel-1);
  }
  // scaling
  if (modelPtr_->solveType()==1) {
    gotHint = (getHintParam(OsiDoScale,takeHint,strength));
    assert (gotHint);
    if (strength==OsiHintIgnore||takeHint)
      solver.scaling(1);
    else
      solver.scaling(0);
  } else {
    solver.scaling(0);
  }
  ClpDualRowSteepest steep;
  solver.setDualRowPivotAlgorithm(steep);
  // sort out hints;
  // algorithm -1 force dual, +1 force primal
  int algorithm = -1;
  gotHint = (getHintParam(OsiDoDualInResolve,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
    algorithm = takeHint ? -1 : 1;
  //solver.saveModel("save.bad");
  // presolve
  gotHint = (getHintParam(OsiDoPresolveInResolve,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    ClpPresolve pinfo;
    ClpSimplex * model2 = pinfo.presolvedModel(solver,1.0e-8);
    if (!model2) {
      // problem found to be infeasible - whats best?
      model2 = &solver;
    }
    // change from 200
    model2->factorization()->maximumPivots(100+model2->numberRows()/50);
    if (algorithm<0) {
      // up dual bound for safety
      //model2->setDualBound(1.0e10);
      model2->dual();
      // check if clp thought it was in a loop
      if (model2->status()==3&&
	  model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->primal();
      }
    } else {
      // up infeasibility cost for safety
      //model2->setInfeasibilityCost(1.0e10);
      model2->primal();
      // check if clp thought it was in a loop
      if (model2->status()==3
	  &&model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->dual();
      }
    }
    if (model2!=&solver) {
      pinfo.postsolve(true);
    
      delete model2;
      // later try without (1) and check duals before solve
      solver.primal(1);
      lastAlgorithm_=1; // primal
    }
    //if (solver.numberIterations())
    //printf("****** iterated %d\n",solver.numberIterations());
  } else {
    if (algorithm<0) {
      //printf("doing dual\n");
      solver.dual();
      lastAlgorithm_=2; // dual
      // check if clp thought it was in a loop
      if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
	// switch algorithm
	solver.primal();
	lastAlgorithm_=1; // primal
	if (solver.status()==3&&
	    solver.numberIterations()<solver.maximumIterations()) {
	  printf("in trouble - try all slack\n");
	  CoinWarmStartBasis allSlack;
	  setBasis(allSlack,&solver);
	  solver.primal();
	  if (solver.status()==3&&
	      solver.numberIterations()<solver.maximumIterations()) {
	    printf("Real real trouble - treat as infeasible\n");
	    solver.setProblemStatus(1);
	  }
	}
      }
    } else {
      //printf("doing primal\n");
      solver.primal();
      lastAlgorithm_=1; // primal
      // check if clp thought it was in a loop
      if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
	// switch algorithm
	solver.dual();
	lastAlgorithm_=2; // dual
      }
    }
  }
  basis_ = getBasis(&solver);
  //basis_.print();
  solver.messageHandler()->setLogLevel(saveMessageLevel);
  solver.returnModel(*modelPtr_);
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  //solver.setSolveType(saveSolveType);
#else
  int saveSolveType=modelPtr_->solveType();
  bool doingPrimal = modelPtr_->algorithm()>0;
  if (saveSolveType==2) {
    disableSimplexInterface();
  }
  //modelPtr_->setSolveType(1);
  // Set message handler to have same levels etc
  bool oldDefault;
  CoinMessageHandler * saveHandler = modelPtr_->pushMessageHandler(handler_,oldDefault);
  //basis_.print();
  setBasis(basis_,modelPtr_);
  // set reasonable defaults
  bool takeHint;
  OsiHintStrength strength;
  // Switch off printing if asked to
  bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  int saveMessageLevel=messageHandler()->logLevel();
  if (strength!=OsiHintIgnore&&takeHint) {
    if (saveMessageLevel)
      modelPtr_->messageHandler()->setLogLevel(saveMessageLevel-1);
  }
  // scaling
  if (modelPtr_->solveType()==1) {
    gotHint = (getHintParam(OsiDoScale,takeHint,strength));
    assert (gotHint);
    if (strength==OsiHintIgnore||takeHint)
      modelPtr_->scaling(1);
    else
      modelPtr_->scaling(0);
  } else {
    modelPtr_->scaling(0);
  }
  ClpDualRowSteepest steep;
  modelPtr_->setDualRowPivotAlgorithm(steep);
  // sort out hints;
  // algorithm -1 force dual, +1 force primal
  int algorithm = -1;
  gotHint = (getHintParam(OsiDoDualInResolve,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
    algorithm = takeHint ? -1 : 1;
  //modelPtr_->saveModel("save.bad");
  // presolve
  gotHint = (getHintParam(OsiDoPresolveInResolve,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    ClpPresolve pinfo;
    ClpSimplex * model2 = pinfo.presolvedModel(*modelPtr_,1.0e-8);
    if (!model2) {
      // problem found to be infeasible - whats best?
      model2 = modelPtr_;
    }
    // change from 200
    model2->factorization()->maximumPivots(100+model2->numberRows()/50);
    if (algorithm<0) {
      // up dual bound for safety
      //model2->setDualBound(1.0e10);
      model2->dual();
      // check if clp thought it was in a loop
      if (model2->status()==3&&
	  model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->primal();
      }
    } else {
      // up infeasibility cost for safety
      //model2->setInfeasibilityCost(1.0e10);
      model2->primal();
      // check if clp thought it was in a loop
      if (model2->status()==3
	  &&model2->numberIterations()<model2->maximumIterations()) {
	// switch algorithm
	model2->dual();
      }
    }
    if (model2!=modelPtr_) {
      pinfo.postsolve(true);
    
      delete model2;
      // later try without (1) and check duals before solve
      modelPtr_->primal(1);
      lastAlgorithm_=1; // primal
    }
    //if (modelPtr_->numberIterations())
    //printf("****** iterated %d\n",modelPtr_->numberIterations());
  } else {
    if (algorithm<0) {
      //printf("doing dual\n");
      modelPtr_->dual();
      lastAlgorithm_=2; // dual
      // check if clp thought it was in a loop
      if (modelPtr_->status()==3&&modelPtr_->numberIterations()<modelPtr_->maximumIterations()) {
	// switch algorithm
	modelPtr_->primal();
	lastAlgorithm_=1; // primal
	if (modelPtr_->status()==3&&
	    modelPtr_->numberIterations()<modelPtr_->maximumIterations()) {
	  printf("in trouble - try all slack\n");
	  CoinWarmStartBasis allSlack;
	  setBasis(allSlack,modelPtr_);
	  modelPtr_->primal();
	  if (modelPtr_->status()==3&&
	      modelPtr_->numberIterations()<modelPtr_->maximumIterations()) {
	    printf("Real real trouble - treat as infeasible\n");
	    modelPtr_->setProblemStatus(1);
	  }
	}
      }
    } else {
      //printf("doing primal\n");
      modelPtr_->primal();
      lastAlgorithm_=1; // primal
      // check if clp thought it was in a loop
      if (modelPtr_->status()==3&&modelPtr_->numberIterations()<modelPtr_->maximumIterations()) {
	// switch algorithm
	modelPtr_->dual();
	lastAlgorithm_=2; // dual
      }
    }
  }
  basis_ = getBasis(modelPtr_);
  //basis_.print();
  //modelPtr_->messageHandler()->setLogLevel(saveMessageLevel);
  modelPtr_->popMessageHandler(saveHandler,oldDefault);
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  //modelPtr_->setSolveType(saveSolveType);
#endif
}
#ifndef NDEBUG
// For errors to make sure print to screen
// only called in debug mode
static void indexError(int index,
			std::string methodName)
{
  std::cerr<<"Illegal index "<<index<<" in OsiClpSolverInterface::"<<methodName<<std::endl;
  throw CoinError("Illegal index",methodName,"OsiClpSolverInterface");
}
#endif
//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiClpSolverInterface::setIntParam(OsiIntParam key, int value)
{
   std::map<OsiIntParam, ClpIntParam>::const_iterator clpkey =
      intParamMap_.find(key);
   if (clpkey != intParamMap_.end() ) {
      return modelPtr_->setIntParam(clpkey->second, value);
   }
   return false;
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setDblParam(OsiDblParam key, double value)
{
   std::map<OsiDblParam, ClpDblParam>::const_iterator clpkey =
      dblParamMap_.find(key);
   if (clpkey != dblParamMap_.end() ) {
     if (key==OsiDualObjectiveLimit||key==OsiPrimalObjectiveLimit)
       value = modelPtr_->optimizationDirection()*value;
      return modelPtr_->setDblParam(clpkey->second, value);
   }
   return false;
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
   std::map<OsiStrParam, ClpStrParam>::const_iterator clpkey =
      strParamMap_.find(key);
   if (clpkey != strParamMap_.end() ) {
      return modelPtr_->setStrParam(clpkey->second, value);
   }
   return false;
}


//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getIntParam(OsiIntParam key, int& value) const 
{
   std::map<OsiIntParam, ClpIntParam>::const_iterator clpkey =
      intParamMap_.find(key);
   if (clpkey != intParamMap_.end() ) {
      return modelPtr_->getIntParam(clpkey->second, value);
   }
   return false;
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
   std::map<OsiDblParam, ClpDblParam>::const_iterator clpkey =
      dblParamMap_.find(key);
   if (clpkey != dblParamMap_.end() ) {
      bool condition =  modelPtr_->getDblParam(clpkey->second, value);
      if (key==OsiDualObjectiveLimit||key==OsiPrimalObjectiveLimit)
	value = modelPtr_->optimizationDirection()*value;
      return condition;
   }
   return false;
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  if ( key==OsiSolverName ) {
    value = "clp";
    return true;
  }
   std::map<OsiStrParam, ClpStrParam>::const_iterator clpkey =
      strParamMap_.find(key);
   if (clpkey != strParamMap_.end() ) {
      return modelPtr_->getStrParam(clpkey->second, value);
   }
   return false;
}


//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiClpSolverInterface::isAbandoned() const
{
  // not sure about -1 (should not happen)
  return (modelPtr_->status()==4||modelPtr_->status()==-1);
}

bool OsiClpSolverInterface::isProvenOptimal() const
{

  const int stat = modelPtr_->status();
  return (stat == 0);
}

bool OsiClpSolverInterface::isProvenPrimalInfeasible() const
{

  const int stat = modelPtr_->status();
  if (stat != 1)
     return false;
  return true;
}

bool OsiClpSolverInterface::isProvenDualInfeasible() const
{
  const int stat = modelPtr_->status();
  return stat == 2;
}
/* 
   NOTE - Coding if limit > 1.0e30 says that 1.0e29 is loose bound
   so all maximization tests are changed 
*/
bool OsiClpSolverInterface::isPrimalObjectiveLimitReached() const
{
  double limit = 0.0;
  modelPtr_->getDblParam(ClpPrimalObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = modelPtr_->objectiveValue();
  int maxmin = (int) modelPtr_->optimizationDirection();

  switch (lastAlgorithm_) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj < limit) /*minim*/ : (-obj < limit) /*maxim*/;
   case 2: // dual simplex
     if (modelPtr_->status() == 0) // optimal
	return maxmin > 0 ? (obj < limit) /*minim*/ : (-obj < limit) /*maxim*/;
     return false;
   case 1: // primal simplex
     return maxmin > 0 ? (obj < limit) /*minim*/ : (-obj < limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiClpSolverInterface::isDualObjectiveLimitReached() const
{

  double limit = 0.0;
  modelPtr_->getDblParam(ClpDualObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = modelPtr_->objectiveValue();
  int maxmin = (int) modelPtr_->optimizationDirection();

  switch (lastAlgorithm_) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj > limit) /*minim*/ : (-obj > limit) /*maxim*/;
   case 1: // primal simplex
     if (modelPtr_->status() == 0) // optimal
	return maxmin > 0 ? (obj > limit) /*minim*/ : (-obj > limit) /*maxim*/;
     return false;
   case 2: // dual simplex
     if (modelPtr_->status() != 0 && modelPtr_->status() != 3)
	// over dual limit
	return true;
     return maxmin > 0 ? (obj > limit) /*minim*/ : (-obj > limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiClpSolverInterface::isIterationLimitReached() const
{
  const int stat = modelPtr_->status();
  return (stat == 3);
}

//#############################################################################
// WarmStart related methods
//#############################################################################
CoinWarmStart *OsiClpSolverInterface::getEmptyWarmStart () const
  { return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartBasis())) ; }

CoinWarmStart* OsiClpSolverInterface::getWarmStart() const
{

  return new CoinWarmStartBasis(basis_);
}

//-----------------------------------------------------------------------------

bool OsiClpSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{

  const CoinWarmStartBasis* ws =
    dynamic_cast<const CoinWarmStartBasis*>(warmstart);

  if (ws) {
    basis_ = CoinWarmStartBasis(*ws);
    return true;
  } else if (!warmstart) {
    // create from current basis
    basis_ = getBasis(modelPtr_);
    return true;
  } else {
    return false;
  }
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiClpSolverInterface::markHotStart()
{
  delete ws_;
  ws_ = dynamic_cast<CoinWarmStartBasis*>(getWarmStart());
  int numberRows = modelPtr_->numberRows();
  rowActivity_= new double[numberRows];
  memcpy(rowActivity_,modelPtr_->primalRowSolution(),
	 numberRows*sizeof(double));
  int numberColumns = modelPtr_->numberColumns();
  columnActivity_= new double[numberColumns];
  memcpy(columnActivity_,modelPtr_->primalColumnSolution(),
	 numberColumns*sizeof(double));

}

void OsiClpSolverInterface::solveFromHotStart()
{
  setWarmStart(ws_);
  int numberRows = modelPtr_->numberRows();
  memcpy(modelPtr_->primalRowSolution(),
	 rowActivity_,numberRows*sizeof(double));
  int numberColumns = modelPtr_->numberColumns();
  memcpy(modelPtr_->primalColumnSolution(),columnActivity_,
	 numberColumns*sizeof(double));
  bool takeHint;
  OsiHintStrength strength;
  // Switch off printing if asked to
  bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  int saveMessageLevel=messageHandler()->logLevel();
  if (strength!=OsiHintIgnore&&takeHint) {
    if (saveMessageLevel)
      messageHandler()->setLogLevel(saveMessageLevel-1);
  }
  messageHandler()->setLogLevel(saveMessageLevel);

  modelPtr_->getIntParam(ClpMaxNumIteration,itlimOrig_);
  int itlim;
  modelPtr_->getIntParam(ClpMaxNumIterationHotStart, itlim);
  modelPtr_->setIntParam(ClpMaxNumIteration,itlim);
  resolve();
  modelPtr_->setIntParam(ClpMaxNumIteration,itlimOrig_);
}

void OsiClpSolverInterface::unmarkHotStart()
{
  delete ws_;
  ws_ = NULL;
  delete [] rowActivity_;
  delete [] columnActivity_;
  rowActivity_=NULL;
  columnActivity_=NULL;
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
const char * OsiClpSolverInterface::getRowSense() const
{
  extractSenseRhsRange();
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiClpSolverInterface::getRightHandSide() const
{
  extractSenseRhsRange();
  return rhs_;
}
//------------------------------------------------------------------
const double * OsiClpSolverInterface::getRowRange() const
{
  extractSenseRhsRange();
  return rowrange_;
}
//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------
bool OsiClpSolverInterface::isContinuous(int colNumber) const
{
  if ( integerInformation_==NULL ) return true;
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (colNumber<0||colNumber>=n) {
    indexError(colNumber,"isContinuous");
  }
#endif
  if ( integerInformation_[colNumber]==0 ) return true;
  return false;
}
//------------------------------------------------------------------

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------
const CoinPackedMatrix * OsiClpSolverInterface::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL ) {
    matrixByRow_ = new CoinPackedMatrix(); 
    matrixByRow_->reverseOrderedCopyOf(*modelPtr_->matrix());
    matrixByRow_->removeGaps();
#if 0
    CoinPackedMatrix back;
    std::cout<<"start check"<<std::endl;
    back.reverseOrderedCopyOf(*matrixByRow_);
    modelPtr_->matrix()->isEquivalent2(back);
    std::cout<<"stop check"<<std::endl;
#endif
  }
  return matrixByRow_;
}

const CoinPackedMatrix * OsiClpSolverInterface::getMatrixByCol() const
{
  return modelPtr_->matrix();
}

//------------------------------------------------------------------
std::vector<double*> OsiClpSolverInterface::getDualRays(int maxNumRays) const
{
  return std::vector<double*>(1, modelPtr_->infeasibilityRay());
}
//------------------------------------------------------------------
std::vector<double*> OsiClpSolverInterface::getPrimalRays(int maxNumRays) const
{
  return std::vector<double*>(1, modelPtr_->unboundedRay());
}
//------------------------------------------------------------------

//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  double * lower = modelPtr_->columnLower();
  double * upper = modelPtr_->columnUpper();
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
#endif
  while (indexFirst != indexLast) {
    const int iCol=*indexFirst++;
#ifndef NDEBUG
    if (iCol<0||iCol>=n) {
      indexError(iCol,"setColSetBounds");
    }
#endif
    lower[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
    upper[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
  }
  if (modelPtr_->solveType()==2) {
    // directly into code as well
    double * lower = modelPtr_->lowerRegion(1);
    double * upper = modelPtr_->upperRegion(1);
    while (indexFirst != indexLast) {
      const int iCol=*indexFirst++;
      lower[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
      upper[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
    }
    
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  // *TEST*
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (i<0||i>=n) {
    indexError(i,"setRowType");
  }
#endif
  double lower, upper;
  convertSenseToBound(sense, rightHandSide, range, lower, upper);
  setRowBounds(i, lower, upper);
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
#endif
  double * lower = modelPtr_->rowLower();
  double * upper = modelPtr_->rowUpper();
  const int len = indexLast - indexFirst;
  while (indexFirst != indexLast) {
    const int iRow=*indexFirst++;
#ifndef NDEBUG
    if (iRow<0||iRow>=n) {
      indexError(iRow,"setRowSetBounds");
    }
#endif
    lower[iRow]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
  }
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    indexFirst -= len;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      convertBoundToSense(lower[iRow], upper[iRow],
			  rowsense_[iRow], rhs_[iRow], rowrange_[iRow]);
    }
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
#endif
  double * lower = modelPtr_->rowLower();
  double * upper = modelPtr_->rowUpper();
  const int len = indexLast - indexFirst;
  while (indexFirst != indexLast) {
    const int iRow= *indexFirst++;
#ifndef NDEBUG
    if (iRow<0||iRow>=n) {
      indexError(iRow,"isContinuous");
    }
#endif
      if (rangeList){
	convertSenseToBound(*senseList++, *rhsList++, *rangeList++,
			    lower[iRow], upper[iRow]);
      } else {
	convertSenseToBound(*senseList++, *rhsList++, 0,
			    lower[iRow], upper[iRow]);
      }
  }
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    indexFirst -= len;
    senseList -= len;
    rhsList -= len;
    if (rangeList)
       rangeList -= len;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      rowsense_[iRow] = *senseList++;
      rhs_[iRow] = *rhsList++;
      if (rangeList)
	 rowrange_[iRow] = *rangeList++;
    }
  }
}
//#############################################################################
void
OsiClpSolverInterface::setContinuous(int index)
{

  if (integerInformation_) {
#ifndef NDEBUG
    int n = modelPtr_->numberColumns();
    if (index<0||index>=n) {
      indexError(index,"setContinuous");
    }
#endif
    integerInformation_[index]=0;
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setInteger(int index)
{
  if (!integerInformation_) {
    integerInformation_ = new char[modelPtr_->numberColumns()];
    CoinFillN ( integerInformation_, modelPtr_->numberColumns(),(char) 0);
  }
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (index<0||index>=n) {
    indexError(index,"setInteger");
  }
#endif
  integerInformation_[index]=1;
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setContinuous(const int* indices, int len)
{
  if (integerInformation_) {
#ifndef NDEBUG
    int n = modelPtr_->numberColumns();
#endif
    int i;
    for (i=0; i<len;i++) {
      int colNumber = indices[i];
#ifndef NDEBUG
      if (colNumber<0||colNumber>=n) {
	indexError(colNumber,"setContinuous");
      }
#endif
      integerInformation_[colNumber]=0;
    }
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setInteger(const int* indices, int len)
{
  if (!integerInformation_) {
    integerInformation_ = new char[modelPtr_->numberColumns()];
    CoinFillN ( integerInformation_, modelPtr_->numberColumns(),(char) 0);
  }
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
#endif
  int i;
  for (i=0; i<len;i++) {
    int colNumber = indices[i];
#ifndef NDEBUG
    if (colNumber<0||colNumber>=n) {
      indexError(colNumber,"setInteger");
    }
#endif
    integerInformation_[colNumber]=1;
  }
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setColSolution(const double * cs) 
{
  CoinDisjointCopyN(cs,modelPtr_->numberColumns(),
		    modelPtr_->primalColumnSolution());
  if (modelPtr_->solveType()==2) {
    // directly into code as well
    CoinDisjointCopyN(cs,modelPtr_->numberColumns(),
		      modelPtr_->solutionRegion(1));
  }
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowPrice(const double * rs) 
{
  CoinDisjointCopyN(rs,modelPtr_->numberRows(),
		    modelPtr_->dualRowSolution());
  if (modelPtr_->solveType()==2) {
    // directly into code as well (? sign )
    CoinDisjointCopyN(rs,modelPtr_->numberRows(),
		      modelPtr_->djRegion(0));
  }
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiClpSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  int numberColumns = modelPtr_->numberColumns();
  modelPtr_->resize(modelPtr_->numberRows(),numberColumns+1);
  linearObjective_ = modelPtr_->objective();
  basis_.resize(modelPtr_->numberRows(),numberColumns+1);
  setColBounds(numberColumns,collb,colub);
  setObjCoeff(numberColumns,obj);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendCol(vec);
  if (integerInformation_) {
    char * temp = new char[numberColumns+1];
    memcpy(temp,integerInformation_,numberColumns*sizeof(char));
    delete [] integerInformation_;
    integerInformation_ = temp;
    integerInformation_[numberColumns]=0;
  }
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  int numberColumns = modelPtr_->numberColumns();
  modelPtr_->resize(modelPtr_->numberRows(),numberColumns+numcols);
  linearObjective_ = modelPtr_->objective();
  basis_.resize(modelPtr_->numberRows(),numberColumns+numcols);
  double * lower = modelPtr_->columnLower()+numberColumns;
  double * upper = modelPtr_->columnUpper()+numberColumns;
  double * objective = modelPtr_->objective()+numberColumns;
  int iCol;
  for (iCol = 0; iCol < numcols; iCol++) {
    lower[iCol]= forceIntoRange(collb[iCol], -OsiClpInfinity, OsiClpInfinity);
    upper[iCol]= forceIntoRange(colub[iCol], -OsiClpInfinity, OsiClpInfinity);
    objective[iCol] = obj[iCol];
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendCols(numcols,cols);
  if (integerInformation_) {
    char * temp = new char[numberColumns+numcols];
    memcpy(temp,integerInformation_,numberColumns*sizeof(char));
    delete [] integerInformation_;
    integerInformation_ = temp;
    for (iCol = 0; iCol < numcols; iCol++) 
      integerInformation_[numberColumns+iCol]=0;
  }
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  modelPtr_->deleteColumns(num,columnIndices);
  basis_.deleteColumns(num,columnIndices);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  setRowBounds(numberRows,rowlb,rowub);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  double rowlb, rowub;
  convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
  setRowBounds(numberRows,rowlb,rowub);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    lower[iRow]= forceIntoRange(rowlb[iRow], -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(rowub[iRow], -OsiClpInfinity, OsiClpInfinity);
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRows(numrows,rows);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    double rowlb, rowub;
    convertSenseToBound(rowsen[iRow], rowrhs[iRow], rowrng[iRow], 
			rowlb, rowub);
    lower[iRow]= forceIntoRange(rowlb, -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(rowub, -OsiClpInfinity, OsiClpInfinity);
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRows(numrows,rows);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  modelPtr_->deleteRows(num,rowIndices);
  basis_.deleteRows(num,rowIndices);
  freeCachedResults();
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiClpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  modelPtr_->loadProblem(matrix, collb, colub, obj, rowlb, rowub);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();

}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     double*& rowlb, double*& rowub)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  modelPtr_->loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
  linearObjective_ = modelPtr_->objective();
  
  freeCachedResults();
  delete matrix;   matrix = NULL;
  delete[] collb;  collb = NULL;
  delete[] colub;  colub = NULL;
  delete[] obj;    obj = NULL;
  delete[] rowlb;  rowlb = NULL;
  delete[] rowub;  rowub = NULL;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  assert( rowsen != NULL );
  assert( rowrhs != NULL );
  int numrows = matrix.getNumRows();
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
  }
  modelPtr_->loadProblem(matrix, collb, colub, obj, rowlb, rowub);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();
  delete [] rowlb;
  delete [] rowub;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     char*& rowsen, double*& rowrhs,
				     double*& rowrng)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
  linearObjective_ = modelPtr_->objective();
  delete matrix;   matrix = NULL;
  delete[] collb;  collb = NULL;
  delete[] colub;  colub = NULL;
  delete[] obj;    obj = NULL;
  delete[] rowsen; rowsen = NULL;
  delete[] rowrhs; rowrhs = NULL;
  delete[] rowrng; rowrng = NULL;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const int numcols, const int numrows,
				   const CoinBigIndex * start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  modelPtr_->loadProblem(numcols, numrows, start,  index,
	    value, collb, colub, obj,
	    rowlb,  rowub);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();
}
//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const int numcols, const int numrows,
				   const CoinBigIndex * start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  delete [] integerInformation_;
  integerInformation_=NULL;
  assert( rowsen != NULL );
  assert( rowrhs != NULL );
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
  }
  modelPtr_->loadProblem(numcols, numrows, start,  index,
			 value, collb, colub, obj,
			 rowlb,  rowub);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();
  delete[] rowlb;
  delete[] rowub;
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------

void OsiClpSolverInterface::writeMps(const char * filename,
				     const char * extension,
				     double objSense) const
{
  std::string f(filename);
  std::string e(extension);
  std::string fullname;
  if (e!="") {
    fullname = f + "." + e;
  } else {
    // no extension so no trailing period
    fullname = f;
  }
  // Fall back on Osi version - without names
  OsiSolverInterface::writeMpsNative(fullname.c_str(), 
				     NULL, NULL,0,2,objSense);
}

int 
OsiClpSolverInterface::writeMpsNative(const char *filename, 
		  const char ** rowNames, const char ** columnNames,
		  int formatType,int numberAcross,double objSense) const 
{
  return OsiSolverInterface::writeMpsNative(filename, rowNames, columnNames,
			       formatType, numberAcross,objSense);
}

//#############################################################################
// CLP specific public interfaces
//#############################################################################

ClpSimplex * OsiClpSolverInterface::getModelPtr() const
{
  freeCachedResults();
  return modelPtr_;
}

//------------------------------------------------------------------- 

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::OsiClpSolverInterface ()
:
OsiSolverInterface(),
linearObjective_(NULL),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
rowActivity_(NULL),
columnActivity_(NULL),
matrixByRow_(NULL),
integerInformation_(NULL)
{
  modelPtr_=NULL;
  notOwned_=false;
  reset();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * OsiClpSolverInterface::clone(bool CopyData) const
{
   if (CopyData) {
      return new OsiClpSolverInterface(*this);
   } else {
      return new OsiClpSolverInterface();
   }
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::OsiClpSolverInterface (
                  const OsiClpSolverInterface & rhs)
:
OsiSolverInterface(rhs),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
rowActivity_(NULL),
columnActivity_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL)
{
  if ( rhs.modelPtr_  ) 
    modelPtr_ = new ClpSimplex(*rhs.modelPtr_);
  else
    modelPtr_ = new ClpSimplex();
  linearObjective_ = modelPtr_->objective();
  if ( rhs.ws_ ) 
    ws_ = new CoinWarmStartBasis(*rhs.ws_);
  basis_ = rhs.basis_;
  if (rhs.integerInformation_) {
    int numberColumns = modelPtr_->numberColumns();
    integerInformation_ = new char[numberColumns];
    memcpy(integerInformation_,rhs.integerInformation_,
	   numberColumns*sizeof(char));
  }
  saveData_ = rhs.saveData_;
  fillParamMaps();
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
}

// Borrow constructor - only delete one copy
OsiClpSolverInterface::OsiClpSolverInterface (ClpSimplex * rhs)
:
OsiSolverInterface(),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
rowActivity_(NULL),
columnActivity_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL)
{
  modelPtr_ = rhs;
  linearObjective_ = modelPtr_->objective();
  if (rhs) {
    notOwned_=true;

    if (rhs->integerInformation()) {
      int numberColumns = modelPtr_->numberColumns();
      integerInformation_ = new char[numberColumns];
      memcpy(integerInformation_,rhs->integerInformation(),
	     numberColumns*sizeof(char));
    }
  }
  fillParamMaps();
}
    
// Releases so won't error
void 
OsiClpSolverInterface::releaseClp()
{
  modelPtr_=NULL;
  notOwned_=false;
}
    

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::~OsiClpSolverInterface ()
{
  freeCachedResults();
  if (!notOwned_)
    delete modelPtr_;
  delete ws_;
  delete [] rowActivity_;
  delete [] columnActivity_;
  delete [] integerInformation_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiClpSolverInterface &
OsiClpSolverInterface::operator=(const OsiClpSolverInterface& rhs)
{
  if (this != &rhs) {    
    OsiSolverInterface::operator=(rhs);
    freeCachedResults();
    if (!notOwned_)
      delete modelPtr_;
    delete ws_;
    if ( rhs.modelPtr_  ) 
      modelPtr_ = new ClpSimplex(*rhs.modelPtr_);
    notOwned_=false;
    linearObjective_ = modelPtr_->objective();
    
    if ( rhs.ws_ ) 
      ws_ = new CoinWarmStartBasis(*rhs.ws_);
    delete [] rowActivity_;
    delete [] columnActivity_;
    rowActivity_=NULL;
    columnActivity_=NULL;
    basis_ = rhs.basis_;
    intParamMap_ = rhs.intParamMap_;
    dblParamMap_ = rhs.dblParamMap_;
    strParamMap_ = rhs.strParamMap_;
    messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiClpSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  const CoinPackedVector & row=rowCut.row();
  addRow(row ,  rowCut.lb(),rowCut.ub());
}
/* Apply a collection of row cuts which are all effective.
   applyCuts seems to do one at a time which seems inefficient.
*/
void 
OsiClpSolverInterface::applyRowCuts(int numberCuts, const OsiRowCut * cuts)
{
  int i;
  if (!numberCuts)
    return;

  const CoinPackedVectorBase * * rows
    =     new const CoinPackedVectorBase * [numberCuts];
  double * rowlb = new double [numberCuts];
  double * rowub = new double [numberCuts];
  for (i=0;i<numberCuts;i++) {
    rowlb[i] = cuts[i].lb();
    rowub[i] = cuts[i].ub();
    rows[i] = &cuts[i].row();
#ifdef TAKEOUT
    if (rows[i]->getNumElements()==10||rows[i]->getNumElements()==15)
      printf("ApplyCuts %d size %d\n",getNumRows()+i,rows[i]->getNumElements());
#endif
  }
  addRows(numberCuts,rows,rowlb,rowub);
  delete [] rows;
  delete [] rowlb;
  delete [] rowub;

}

//-----------------------------------------------------------------------------

void OsiClpSolverInterface::applyColCut( const OsiColCut & cc )
{
  double * lower = modelPtr_->columnLower();
  double * upper = modelPtr_->columnUpper();
  const CoinPackedVector & lbs = cc.lbs();
  const CoinPackedVector & ubs = cc.ubs();
  int i;

  for ( i=0; i<lbs.getNumElements(); i++ ) {
    int iCol = lbs.getIndices()[i];
    double value = lbs.getElements()[i];
    if ( value > lower[iCol] )
      lower[iCol]= value;
  }
  for ( i=0; i<ubs.getNumElements(); i++ ) {
    int iCol = ubs.getIndices()[i];
    double value = ubs.getElements()[i];
    if ( value < upper[iCol] )
      upper[iCol]= value;
  }
}
//#############################################################################
// Private methods
//#############################################################################


//------------------------------------------------------------------- 

void OsiClpSolverInterface::freeCachedResults() const
{  
  delete [] rowsense_;
  delete [] rhs_;
  delete [] rowrange_;
  delete matrixByRow_;
  delete ws_;
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
  matrixByRow_=NULL;
  ws_ = NULL;
  if (modelPtr_&&modelPtr_->clpMatrix())
    modelPtr_->clpMatrix()->refresh(modelPtr_); // make sure all clean
}

//------------------------------------------------------------------
void OsiClpSolverInterface::extractSenseRhsRange() const
{
  if (rowsense_ == NULL) {
    // all three must be NULL
    assert ((rhs_ == NULL) && (rowrange_ == NULL));
    
    int nr=modelPtr_->numberRows();
    if ( nr!=0 ) {
      rowsense_ = new char[nr];
      rhs_ = new double[nr];
      rowrange_ = new double[nr];
      std::fill(rowrange_,rowrange_+nr,0.0);
      
      const double * lb = modelPtr_->rowLower();
      const double * ub = modelPtr_->rowUpper();
      
      int i;
      for ( i=0; i<nr; i++ ) {
        convertBoundToSense(lb[i], ub[i], rowsense_[i], rhs_[i], rowrange_[i]);
      }
    }
  }
}
// Set language
void 
OsiClpSolverInterface::newLanguage(CoinMessages::Language language)
{
  modelPtr_->newLanguage(language);
  OsiSolverInterface::newLanguage(language);
}
//#############################################################################

void
OsiClpSolverInterface::fillParamMaps()
{
   intParamMap_[OsiMaxNumIteration]         = ClpMaxNumIteration;
   intParamMap_[OsiMaxNumIterationHotStart] = ClpMaxNumIterationHotStart;
   intParamMap_[OsiLastIntParam]            = ClpLastIntParam;

   dblParamMap_[OsiDualObjectiveLimit]   = ClpDualObjectiveLimit;
   dblParamMap_[OsiPrimalObjectiveLimit] = ClpPrimalObjectiveLimit;
   dblParamMap_[OsiDualTolerance]        = ClpDualTolerance;
   dblParamMap_[OsiPrimalTolerance]      = ClpPrimalTolerance;
   dblParamMap_[OsiObjOffset]            = ClpObjOffset;
   dblParamMap_[OsiLastDblParam]         = ClpLastDblParam;

   strParamMap_[OsiProbName]     = ClpProbName;
   strParamMap_[OsiLastStrParam] = ClpLastStrParam;
}
// Warm start
CoinWarmStartBasis
OsiClpSolverInterface::getBasis(ClpSimplex * model) const
{
  int iRow,iColumn;
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  CoinWarmStartBasis basis;
  basis.setSize(numberColumns,numberRows);

  if (model->statusExists()) {
    int lookup[]={0,1,2,3,0,3};
    for (iRow=0;iRow<numberRows;iRow++) {
      int iStatus = model->getRowStatus(iRow);
      iStatus = lookup[iStatus];
      basis.setArtifStatus(iRow,(CoinWarmStartBasis::Status) iStatus);
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int iStatus = model->getColumnStatus(iColumn);
      iStatus = lookup[iStatus];
      basis.setStructStatus(iColumn,(CoinWarmStartBasis::Status) iStatus);
    }
  }
  //basis.print();
  return basis;
}
// Sets up basis
void 
OsiClpSolverInterface::setBasis ( const CoinWarmStartBasis & basis,
				  ClpSimplex * model)
{
  // transform basis to status arrays
  int iRow,iColumn;
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  if (!model->statusExists()) {
    /*
      get status arrays
      ClpBasis would seem to have overheads and we will need
      extra bits anyway.
    */
    model->createStatus();
  }
  CoinWarmStartBasis basis2 = basis;
  // resize if necessary
  basis2.resize(numberRows,numberColumns);
  // move status
  model->createStatus();
  for (iRow=0;iRow<numberRows;iRow++) {
    model->setRowStatus(iRow,
		 (ClpSimplex::Status) basis2.getArtifStatus(iRow));
  }
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    model->setColumnStatus(iColumn,
		    (ClpSimplex::Status) basis2.getStructStatus(iColumn));
  }
}
/* Read an mps file from the given filename (defaults to Osi reader) - returns
   number of errors (see OsiMpsReader class) */
int 
OsiClpSolverInterface::readMps(const char *filename,
			       const char *extension ) 
{
  // Get rid of integer stuff
  delete [] integerInformation_;
  integerInformation_=NULL;
  
  int numberErrors = OsiSolverInterface::readMps(filename,extension);
  // move across integer information
  int numberColumns = modelPtr_->numberColumns();
  int i;
  char * info = new char [numberColumns];
  int numberIntegers=0;
  for (i=0;i<numberColumns;i++) {
    if (isInteger(i)) {
      info[i]=1;
      numberIntegers++;
    } else {
      info[i]=0;
    }
  }
  if (numberIntegers)
    modelPtr_->copyInIntegerInformation(info);
  delete [] info;
  return numberErrors;
}
// Get pointer to array[getNumCols()] of primal solution vector
const double * 
OsiClpSolverInterface::getColSolution() const 
{ 
  if (modelPtr_->solveType()!=2) {
    return modelPtr_->primalColumnSolution();
  } else {
    // simplex interface
    return modelPtr_->solutionRegion(1);
  }
}
  
// Get pointer to array[getNumRows()] of dual prices
const double * 
OsiClpSolverInterface::getRowPrice() const
{ 
  if (modelPtr_->solveType()!=2) {
    return modelPtr_->dualRowSolution();
  } else {
    // simplex interface
    //return modelPtr_->djRegion(0);
    return modelPtr_->dualRowSolution();
  }
}
  
// Get a pointer to array[getNumCols()] of reduced costs
const double * 
OsiClpSolverInterface::getReducedCost() const 
{ 
  if (modelPtr_->solveType()!=2) {
    return modelPtr_->dualColumnSolution();
  } else {
    // simplex interface
    return modelPtr_->djRegion(1);
  }
}

/* Get pointer to array[getNumRows()] of row activity levels (constraint
   matrix times the solution vector */
const double * 
OsiClpSolverInterface::getRowActivity() const 
{ 
  if (modelPtr_->solveType()!=2) {
    return modelPtr_->primalRowSolution();
  } else {
    // simplex interface
    return modelPtr_->solutionRegion(0);
  }
}
double 
OsiClpSolverInterface::getObjValue() const 
{
#if 0
   // This does not pass unitTest when getObjValue is called before solve.
   return modelPtr_->objectiveValue();
#else
   return OsiSolverInterface::getObjValue();
#endif
}

/* Set an objective function coefficient */
void 
OsiClpSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setObjCoeff");
  }
#endif
  linearObjective_[elementIndex] = elementValue;
  if (modelPtr_->solveType()==2) {
    // simplex interface
    modelPtr_->costRegion(1)[elementIndex] = elementValue;
  }
}

/* Set a single column lower bound<br>
   Use -DBL_MAX for -infinity. */
void 
OsiClpSolverInterface::setColLower( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColLower");
  }
#endif
  modelPtr_->columnLower()[elementIndex] = elementValue;
  if (modelPtr_->solveType()==2) {
    // simplex interface
    modelPtr_->lowerRegion(1)[elementIndex] = elementValue;
  }
}
      
/* Set a single column upper bound<br>
   Use DBL_MAX for infinity. */
void 
OsiClpSolverInterface::setColUpper( int elementIndex, double elementValue )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColUpper");
  }
#endif
  modelPtr_->columnUpper()[elementIndex] = elementValue;
  if (modelPtr_->solveType()==2) {
    // simplex interface
    modelPtr_->upperRegion(1)[elementIndex] = elementValue;
  }
}

/* Set a single column lower and upper bound */
void 
OsiClpSolverInterface::setColBounds( int elementIndex,
				     double lower, double upper )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColBounds");
  }
#endif
  modelPtr_->columnLower()[elementIndex] = lower;
  modelPtr_->columnUpper()[elementIndex] = upper;
  if (modelPtr_->solveType()==2) {
    // simplex interface
    modelPtr_->lowerRegion(1)[elementIndex] = lower;
    modelPtr_->upperRegion(1)[elementIndex] = upper;
  }
}
/*Enables normal operation of subsequent functions.
  This method is supposed to ensure that all typical things (like
  reduced costs, etc.) are updated when individual pivots are executed
  and can be queried by other methods 
*/
void 
OsiClpSolverInterface::enableSimplexInterface(bool doingPrimal)
{
  assert (modelPtr_->solveType()==1);
  modelPtr_->setSolveType(2);
  if (doingPrimal)
    modelPtr_->setAlgorithm(1);
  else
    modelPtr_->setAlgorithm(-1);
  modelPtr_->scaling(0);
  // Do initialization
  saveData_ = modelPtr_->saveData();
  // set infeasibility cost up
  modelPtr_->setInfeasibilityCost(1.0e12);
  // probably should save and restore?
  ClpDualRowDantzig dantzig;
  modelPtr_->setDualRowPivotAlgorithm(dantzig);
  ClpPrimalColumnDantzig dantzigP;
  modelPtr_->setPrimalColumnPivotAlgorithm(dantzigP);
  int returnCode=modelPtr_->startup(0);
  assert (!returnCode);
}

//Undo whatever setting changes the above method had to make
void 
OsiClpSolverInterface::disableSimplexInterface()
{
  assert (modelPtr_->solveType()==2);
  // declare optimality anyway  (for message handler)
  modelPtr_->setProblemStatus(0);
  modelPtr_->setSolveType(1);
  modelPtr_->finish();
  modelPtr_->restoreData(saveData_);
  basis_ = getBasis(modelPtr_);
  modelPtr_->setSolveType(1);
}
/* The following two methods may be replaced by the
   methods of OsiSolverInterface using OsiWarmStartBasis if:
   1. OsiWarmStartBasis resize operation is implemented
   more efficiently and
   2. It is ensured that effects on the solver are the same
   
   Returns a basis status of the structural/artificial variables 
*/
void 
OsiClpSolverInterface::getBasisStatus(int* cstat, int* rstat)
{
  int i, n;
  n=modelPtr_->numberRows();
  int lookup[]={0,1,2,3,0,3};
  for (i=0;i<n;i++) 
    rstat[i] = lookup[modelPtr_->getRowStatus(i)];
  n=modelPtr_->numberColumns();
  for (i=0;i<n;i++)
    cstat[i] = lookup[modelPtr_->getColumnStatus(i)];
}

//Set the status of structural/artificial variables 
int 
OsiClpSolverInterface::setBasisStatus(const int* cstat, const int* rstat)
{
  modelPtr_->createStatus();
  int i, n;
  double * lower, * upper, * solution;
  n=modelPtr_->numberRows();
  lower = modelPtr_->rowLower();
  upper = modelPtr_->rowUpper();
  solution = modelPtr_->primalRowSolution();
  for (i=0;i<n;i++) {
    int status = rstat[i];
    if (status<0||status>3)
      status = 3;
    if (lower[i]<-1.0e50&&upper[i]>1.0e50&&status!=1)
      status = 0; // set free if should be
    else if (lower[i]<-1.0e50&&status==3)
      status = 2; // can't be at lower bound
    else if (upper[i]>1.0e50&&status==2)
      status = 3; // can't be at upper bound
    switch (status) {
      // free or superbasic
    case 0:
      if (lower[i]<-1.0e50&&upper[i]>1.0e50) {
	modelPtr_->setRowStatus(i,ClpSimplex::isFree);
	if (fabs(solution[i])>1.0e20)
	  solution[i]=0.0;
      } else {
	modelPtr_->setRowStatus(i,ClpSimplex::superBasic);
	if (fabs(solution[i])>1.0e20)
	  solution[i]=0.0;
      }
      break;
    case 1:
      // basic
      modelPtr_->setRowStatus(i,ClpSimplex::basic);
      break;
    case 2:
      // at upper bound
      solution[i]=upper[i];
      if (upper[i]>lower[i])
	modelPtr_->setRowStatus(i,ClpSimplex::atUpperBound);
      else
	modelPtr_->setRowStatus(i,ClpSimplex::isFixed);
      break;
    case 3:
      // at lower bound
      solution[i]=lower[i];
      if (upper[i]>lower[i])
	modelPtr_->setRowStatus(i,ClpSimplex::atLowerBound);
      else
	modelPtr_->setRowStatus(i,ClpSimplex::isFixed);
      break;
    }
  }
  n=modelPtr_->numberColumns();
  lower = modelPtr_->columnLower();
  upper = modelPtr_->columnUpper();
  solution = modelPtr_->primalColumnSolution();
  for (i=0;i<n;i++) {
    int status = cstat[i];
    if (status<0||status>3)
      status = 3;
    if (lower[i]<-1.0e50&&upper[i]>1.0e50&&status!=1)
      status = 0; // set free if should be
    else if (lower[i]<-1.0e50&&status==3)
      status = 2; // can't be at lower bound
    else if (upper[i]>1.0e50&&status==2)
      status = 3; // can't be at upper bound
    switch (status) {
      // free or superbasic
    case 0:
      if (lower[i]<-1.0e50&&upper[i]>1.0e50) {
	modelPtr_->setColumnStatus(i,ClpSimplex::isFree);
	if (fabs(solution[i])>1.0e20)
	  solution[i]=0.0;
      } else {
	modelPtr_->setColumnStatus(i,ClpSimplex::superBasic);
	if (fabs(solution[i])>1.0e20)
	  solution[i]=0.0;
      }
      break;
    case 1:
      // basic
      modelPtr_->setColumnStatus(i,ClpSimplex::basic);
      break;
    case 2:
      // at upper bound
      solution[i]=upper[i];
      if (upper[i]>lower[i])
	modelPtr_->setColumnStatus(i,ClpSimplex::atUpperBound);
      else
	modelPtr_->setColumnStatus(i,ClpSimplex::isFixed);
      break;
    case 3:
      // at lower bound
      solution[i]=lower[i];
      if (upper[i]>lower[i])
	modelPtr_->setColumnStatus(i,ClpSimplex::atLowerBound);
      else
	modelPtr_->setColumnStatus(i,ClpSimplex::isFixed);
      break;
    }
  }
  // say first time
  modelPtr_->statusOfProblem(true);
  return 0;
}

/* Perform a pivot by substituting a colIn for colOut in the basis. 
   The status of the leaving variable is given in statOut. Where
   1 is to upper bound, -1 to lower bound
*/
int 
OsiClpSolverInterface::pivot(int colIn, int colOut, int outStatus)
{
  assert (modelPtr_->solveType()==2);
  // convert to Clp style (what about flips?)
  if (colIn<0) 
    colIn = modelPtr_->numberColumns()+(-1-colIn);
  if (colOut<0) 
    colOut = modelPtr_->numberColumns()+(-1-colOut);
  // in clp direction of out is reversed
  outStatus = - outStatus;
  // set in clp
  modelPtr_->setDirectionOut(outStatus);
  modelPtr_->setSequenceIn(colIn);
  modelPtr_->setSequenceOut(colOut);
  // do pivot
  modelPtr_->pivot();
  return 0;
}

/* Obtain a result of the primal pivot 
   Outputs: colOut -- leaving column, outStatus -- its status,
   t -- step size, and, if dx!=NULL, *dx -- primal ray direction.
   Inputs: colIn -- entering column, sign -- direction of its change (+/-1).
   Both for colIn and colOut, artificial variables are index by
   the negative of the row index minus 1.
   Return code (for now): 0 -- leaving variable found, 
   -1 -- everything else?
   Clearly, more informative set of return values is required 
*/
int 
OsiClpSolverInterface::primalPivotResult(int colIn, int sign, 
					 int& colOut, int& outStatus, 
					 double& t, CoinPackedVector* dx)
{
  assert (modelPtr_->solveType()==2);
  // convert to Clp style
  if (colIn<0) 
    colIn = modelPtr_->numberColumns()+(-1-colIn);
  // set in clp
  modelPtr_->setDirectionIn(sign);
  modelPtr_->setSequenceIn(colIn);
  modelPtr_->setSequenceOut(-1);
  int returnCode = modelPtr_->primalPivotResult();
  t = modelPtr_->theta();
  int numberColumns = modelPtr_->numberColumns();
  if (dx) {
    double * ray = modelPtr_->unboundedRay();
    if  (ray)
      dx->setFullNonZero(numberColumns,ray);
    else
      printf("No ray?\n");
    delete [] ray;
  }
  outStatus = - modelPtr_->directionOut();
  colOut = modelPtr_->sequenceOut();
  if (colOut>= numberColumns) 
    colOut = -1-(colOut - numberColumns);
  return returnCode;
}

/* Obtain a result of the dual pivot (similar to the previous method)
   Differences: entering variable and a sign of its change are now
   the outputs, the leaving variable and its statuts -- the inputs
   If dx!=NULL, then *dx contains dual ray
   Return code: same
*/
int 
OsiClpSolverInterface::dualPivotResult(int& colIn, int& sign, 
			      int colOut, int outStatus, 
			      double& t, CoinPackedVector* dx)
{
  assert (modelPtr_->solveType()==2);
  abort();
  return 0;
}

//Get the reduced gradient for the cost vector c 
void 
OsiClpSolverInterface::getReducedGradient(
					  double* columnReducedCosts, 
					  double * duals,
					  const double * c)
{
  assert (modelPtr_->solveType()==2);
  // could do this faster with coding inside Clp
  // save current costs
  int numberColumns = modelPtr_->numberColumns();
  double * save = new double [numberColumns];
  memcpy(save,modelPtr_->costRegion(),numberColumns*sizeof(double));
  memcpy(modelPtr_->costRegion(),c,numberColumns*sizeof(double));
  modelPtr_->computeDuals(NULL);
  memcpy(modelPtr_->costRegion(),save,numberColumns*sizeof(double));
  delete [] save;
  int numberRows = modelPtr_->numberRows();
  memcpy(duals,modelPtr_->dualRowSolution(),numberRows*sizeof(double));
  memcpy(columnReducedCosts,modelPtr_->djRegion(1),
	 numberColumns*sizeof(double));
}

/* Set a new objective and apply the old basis so that the
   reduced costs are properly updated  */
void OsiClpSolverInterface::setObjectiveAndRefresh(double* c)
{
  assert (modelPtr_->solveType()==2);
  int numberColumns = modelPtr_->numberColumns();
  memcpy(modelPtr_->objective(),c,numberColumns*sizeof(double));
  memcpy(modelPtr_->costRegion(),c,numberColumns*sizeof(double));
  modelPtr_->computeDuals(NULL);
}

//Get a row of the tableau (slack part in slack if not NULL)
void 
OsiClpSolverInterface::getBInvARow(int row, double* z, double * slack)
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvARow");
  }
#endif
  assert (modelPtr_->solveType()==2);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  CoinIndexedVector * columnArray0 = modelPtr_->columnArray(0);
  CoinIndexedVector * columnArray1 = modelPtr_->columnArray(1);
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
  // put +1 in row
  rowArray1->insert(row,1.0);
  factorization->updateColumnTranspose(rowArray0,rowArray1);
  // put row of tableau in rowArray1 and columnArray0
  modelPtr_->clpMatrix()->transposeTimes(modelPtr_,1.0,
			    rowArray1,columnArray1,columnArray0);
  memcpy(z,columnArray0->denseVector(),
	 modelPtr_->numberColumns()*sizeof(double));
  if (slack) {
    int n = modelPtr_->numberRows();
    double * array = rowArray1->denseVector();
    for (int i=0;i<n;i++) {
      // clp stores slacks as -1.0  (Does not seem to matter - basics should be 1.0)
      slack[i] =  array[i];
    }
    //memcpy(slack,rowArray1->denseVector(),
    //   modelPtr_->numberRows()*sizeof(double));
  }
  // don't need to clear everything always, but doesn't cost
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
}

//Get a row of the basis inverse
void 
OsiClpSolverInterface::getBInvRow(int row, double* z)

{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvRow");
  }
#endif
  assert (modelPtr_->solveType()==2);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // put +1 in row
  rowArray1->insert(row,1.0);
  factorization->updateColumnTranspose(rowArray0,rowArray1);
  memcpy(z,rowArray1->denseVector(),modelPtr_->numberRows()*sizeof(double));
  rowArray1->clear();
}

//Get a column of the tableau
void 
OsiClpSolverInterface::getBInvACol(int col, double* vec)
{
  assert (modelPtr_->solveType()==2);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // get column of matrix
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (col<0||col>=n) {
    indexError(col,"getBInvACol");
  }
#endif
  modelPtr_->unpack(rowArray1,col);
  factorization->updateColumn(rowArray0,rowArray1,false);
  memcpy(vec,rowArray1->denseVector(),modelPtr_->numberRows()*sizeof(double));
  rowArray1->clear();
}

//Get a column of the basis inverse
void 
OsiClpSolverInterface::getBInvCol(int col, double* vec)
{
  assert (modelPtr_->solveType()==2);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (col<0||col>=n) {
    indexError(col,"getBInvCol");
  }
#endif
  // put +1 in row
  rowArray1->insert(col,1.0);
  factorization->updateColumn(rowArray0,rowArray1,false);
  memcpy(vec,rowArray1->denseVector(),modelPtr_->numberRows()*sizeof(double));
  rowArray1->clear();
}

/* Get basic indices (order of indices corresponds to the
   order of elements in a vector retured by getBInvACol() and
   getBInvCol()).
*/
void 
OsiClpSolverInterface::getBasics(int* index)
{
  assert (modelPtr_->solveType()==2);
  assert (index);
  memcpy(index,modelPtr_->pivotVariable(),
	 modelPtr_->numberRows()*sizeof(int));
}
// Resets as if default constructor
void 
OsiClpSolverInterface::reset()
{
  setInitialData(); // clear base class
  freeCachedResults();
  if (!notOwned_)
    delete modelPtr_;
  delete ws_;
  delete [] rowActivity_;
  delete [] columnActivity_;
  delete [] integerInformation_;
  rowActivity_ = NULL;
  columnActivity_ = NULL;
  integerInformation_ = NULL;
  basis_ = CoinWarmStartBasis();
  itlimOrig_=9999999;
  lastAlgorithm_=0;
  notOwned_=false;
  modelPtr_ = new ClpSimplex();
  // This is also deleted by Clp --tkr 7/31/03
  // delete linearObjective_;
  linearObjective_ = NULL;
  fillParamMaps();
}


// Invoke solver's built-in enumeration algorithm
void 
OsiClpSolverInterface::branchAndBound() {
  throw CoinError("Sorry, Clp doesn't implement B&B.  This is because the original Simple Branch and Bound\n\
code has been much developed and moved to Sbb.  If you want a native Branch and Bound code use Sbb",
		  "branchAndBound","OsiClpSolverInterface");
};
