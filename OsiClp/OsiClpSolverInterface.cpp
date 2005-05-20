// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpNonLinearCost.hpp"
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
  int saveOptions = solver.specialOptions();
  solver.setSpecialOptions(saveOptions|64|32768); // go as far as possible
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
    if (strength==OsiHintIgnore||takeHint) {
      if (!solver.scalingFlag())
	solver.scaling(3);
    } else {
      solver.scaling(0);
    }
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
  // Allow for specialOptions_==1+8 forcing saving factorization
  int startFinishOptions=0;
  if (specialOptions_==(1+8)) {
    startFinishOptions =1+2+4; // allow re-use of factorization
  }

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
  // doPrimal set true if any structurals in basis so switch off crash
  if (doPrimal)
    doCrash = -1;
	 
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
      model2->dual(0);
      // check if clp thought it was in a loop
      if (model2->status()==3&&!model2->hitMaximumIterations()) {
	// switch algorithm
	model2->primal();
      }
    } else {
      // up infeasibility cost for safety
      //model2->setInfeasibilityCost(1.0e10);
      model2->primal(1);
      // check if clp thought it was in a loop
      if (model2->status()==3&&!model2->hitMaximumIterations()) {
	// switch algorithm
	model2->dual();
      }
    }
    model2->setPerturbation(savePerturbation);
    if (model2!=&solver) {
      int numberIterations = model2->numberIterations();
      pinfo.postsolve(true);
    
      delete model2;
      //printf("Resolving from postsolved model\n");
      // later try without (1) and check duals before solve
      solver.primal(1);
      solver.setNumberIterations(solver.numberIterations()+numberIterations);
    }
    lastAlgorithm_=1; // primal
    //if (solver.numberIterations())
    //printf("****** iterated %d\n",solver.numberIterations());
  } else {
    // do we want crash
    if (doCrash>0)
      solver.crash(1000.0,2);
    else if (doCrash==0)
      solver.crash(1000.0,0);
    if (algorithm<0)
      doPrimal=false;
    else if (algorithm>0)
      doPrimal=true;
    if (!doPrimal) {
      //printf("doing dual\n");
      solver.dual(0);
      lastAlgorithm_=2; // dual
      // check if clp thought it was in a loop
      if (solver.status()==3&&!solver.hitMaximumIterations()) {
	// switch algorithm
	solver.primal(0);
	lastAlgorithm_=1; // primal
      }
    } else {
      //printf("doing primal\n");
      solver.primal(0);
      lastAlgorithm_=1; // primal
      // check if clp thought it was in a loop
      if (solver.status()==3&&!solver.hitMaximumIterations()) {
	// switch algorithm
	solver.dual(0);
	lastAlgorithm_=2; // dual
      }
    }
  }
  // If scaled feasible but unscaled infeasible take action
  if (!solver.status()&&cleanupScaling_) {
    solver.cleanup(cleanupScaling_);
  }
  basis_ = getBasis(&solver);
  //basis_.print();
  solver.messageHandler()->setLogLevel(saveMessageLevel);
  const double * rowScale2 = solver.rowScale();
  solver.setSpecialOptions(saveOptions);
  if (!rowScale1&&rowScale2) {
    // need to release memory
    solver.setRowScale(NULL);
    solver.setColumnScale(NULL);
  }
  solver.returnModel(*modelPtr_);
  if (startFinishOptions) {
    int save = modelPtr_->logLevel();
    modelPtr_->setLogLevel(0);
    modelPtr_->dual(0,startFinishOptions);
    modelPtr_->setLogLevel(save);
  }
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  // mark so we can pick up objective value quickly
  modelPtr_->upperIn_=0.0;
  time1 = CoinCpuTime()-time1;
  totalTime += time1;
  //std::cout<<time1<<" seconds - total "<<totalTime<<std::endl;
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::resolve()
{
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = (getHintParam(OsiDoInBranchAndCut,takeHint,strength));
  assert (gotHint);
  // mark so we can pick up objective value quickly
  modelPtr_->upperIn_=0.0;
  if ((specialOptions_>>10)==2) {
    // Quick check to see if optimal
    modelPtr_->checkSolutionInternal();
    if (modelPtr_->problemStatus()==0) {
      return;
    }
  }
  int saveSolveType=modelPtr_->solveType();
  bool doingPrimal = modelPtr_->algorithm()>0;
  if (saveSolveType==2) {
    disableSimplexInterface();
  }
  int saveOptions = modelPtr_->specialOptions();
  int startFinishOptions=0;
  if (specialOptions_>=0) {
    if((specialOptions_&1)==0) {
      startFinishOptions=0;
      modelPtr_->setSpecialOptions(saveOptions|(64|1024|32768));
    } else {
      startFinishOptions=1+4;
      if ((specialOptions_&8)!=0)
        startFinishOptions +=2; // allow re-use of factorization
      if((specialOptions_&4)==0||!takeHint) 
        modelPtr_->setSpecialOptions(saveOptions|(64|128|512|1024|4096|32768));
      else
        modelPtr_->setSpecialOptions(saveOptions|(64|128|512|1024|2048|4096|32768));
    }
  } else {
    modelPtr_->setSpecialOptions(saveOptions|64|32768);
  }
  //printf("options %d size %d\n",modelPtr_->specialOptions(),modelPtr_->numberColumns());
  //modelPtr_->setSolveType(1);
  // Set message handler to have same levels etc
  bool oldDefault;
  CoinMessageHandler * saveHandler = modelPtr_->pushMessageHandler(handler_,oldDefault);
  //printf("basis before dual\n");
  //basis_.print();
  setBasis(basis_,modelPtr_);
  // set reasonable defaults
  // Switch off printing if asked to
  gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
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
    if (strength==OsiHintIgnore||takeHint) {
      if (!modelPtr_->scalingFlag())
	modelPtr_->scaling(3);
    } else {
      modelPtr_->scaling(0);
    }
  } else {
    modelPtr_->scaling(0);
  }
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
    if ((modelPtr_->specialOptions()&1024)!=0) {
      pinfo.setDoDual(false);
      pinfo.setDoTripleton(false);
      pinfo.setDoDupcol(false);
      pinfo.setDoDuprow(false);
      pinfo.setDoSingletonColumn(false);
    }
    ClpSimplex * model2 = pinfo.presolvedModel(*modelPtr_,1.0e-8);
    if (!model2) {
      // problem found to be infeasible - whats best?
      model2 = modelPtr_;
    }
    // return number of rows
    int * stats = (int *) getApplicationData();
    if (stats) {
      stats[0]=model2->numberRows();
      stats[1]=model2->numberColumns();
    }
    //printf("rows %d -> %d, columns %d -> %d\n",
    //     modelPtr_->numberRows(),model2->numberRows(),
    //     modelPtr_->numberColumns(),model2->numberColumns());
    // change from 200
    model2->factorization()->maximumPivots(100+model2->numberRows()/50);
    if (algorithm<0) {
      model2->dual();
      // check if clp thought it was in a loop
      if (model2->status()==3&&!model2->hitMaximumIterations()) {
	// switch algorithm
	model2->primal();
      }
    } else {
      model2->primal(1);
      // check if clp thought it was in a loop
      if (model2->status()==3&&!model2->hitMaximumIterations()) {
	// switch algorithm
	model2->dual();
      }
    }
    if (model2!=modelPtr_) {
      int numberIterations = model2->numberIterations();
      int finalStatus=model2->status();
      pinfo.postsolve(true);
    
      delete model2;
      // later try without (1) and check duals before solve
      if (finalStatus!=3&&(finalStatus||modelPtr_->status()==-1)) {
        modelPtr_->primal(1);
        modelPtr_->setNumberIterations(modelPtr_->numberIterations()+numberIterations);
        lastAlgorithm_=1; // primal
        //if (modelPtr_->numberIterations())
        //printf("****** iterated %d\n",modelPtr_->numberIterations());
      }
    }
  } else {
    //modelPtr_->setLogLevel(63);
    //modelPtr_->setDualTolerance(1.0e-7);
    if (algorithm<0) {
      //printf("doing dual\n");
      int savePerturbation = modelPtr_->perturbation();
      if (specialOptions_>=0&&(specialOptions_&2)!=0)
	modelPtr_->setPerturbation(100);
      //modelPtr_->messageHandler()->setLogLevel(1);
      
      //writeMpsNative("bad",NULL,NULL,2,1,1.0);
      if ((modelPtr_->specialOptions()&1024)==0)
        modelPtr_->dual(0,startFinishOptions);
      else
        crunch();
      if (!modelPtr_->problemStatus()&&0) {
        int numberColumns = modelPtr_->numberColumns();
        const double * columnLower = modelPtr_->columnLower();
        const double * columnUpper = modelPtr_->columnUpper();
        int nBad=0;
        for (int i=0;i<numberColumns;i++) {
          if (columnLower[i]==columnUpper[i]&&modelPtr_->getColumnStatus(i)==ClpSimplex::basic) {
            nBad++;
            modelPtr_->setColumnStatus(i,ClpSimplex::isFixed);
          }
        }
        if (nBad) {
          modelPtr_->primal(1);
          printf("%d fixed basic - %d iterations\n",nBad,modelPtr_->numberIterations());
        }
      }
      assert (modelPtr_->objectiveValue()<1.0e100);
      modelPtr_->setPerturbation(savePerturbation);
      lastAlgorithm_=2; // dual
      // check if clp thought it was in a loop
      if (modelPtr_->status()==3&&!modelPtr_->hitMaximumIterations()) {
	modelPtr_->setSpecialOptions(saveOptions);
	// switch algorithm
	//modelPtr_->messageHandler()->setLogLevel(63);
	// Allow for catastrophe
	int saveMax = modelPtr_->maximumIterations();
	int numberIterations = modelPtr_->numberIterations();
	int numberRows = modelPtr_->numberRows();
	int numberColumns = modelPtr_->numberColumns();
	if (modelPtr_->maximumIterations()>100000+numberIterations)
	  modelPtr_->setMaximumIterations(numberIterations + 1000 + 2*numberRows+numberColumns);
	modelPtr_->primal(0,startFinishOptions);
	modelPtr_->setMaximumIterations(saveMax);
	lastAlgorithm_=1; // primal
        if (modelPtr_->status()==3&&!modelPtr_->hitMaximumIterations()) {
	  printf("in trouble - try all slack\n");
	  CoinWarmStartBasis allSlack;
	  setBasis(allSlack,modelPtr_);
	  modelPtr_->dual();
          if (modelPtr_->status()==3&&!modelPtr_->hitMaximumIterations()) {
	    if (modelPtr_->numberPrimalInfeasibilities()) {
	      printf("Real real trouble - treat as infeasible\n");
	      modelPtr_->setProblemStatus(1);
	    } else {
	      printf("Real real trouble - treat as optimal\n");
	      modelPtr_->setProblemStatus(0);
	    }
	  }
	}
      }
      assert (modelPtr_->objectiveValue()<1.0e100);
    } else {
      //printf("doing primal\n");
      modelPtr_->primal(0,startFinishOptions);
      lastAlgorithm_=1; // primal
      // check if clp thought it was in a loop
      if (modelPtr_->status()==3&&!modelPtr_->hitMaximumIterations()) {
	// switch algorithm
	modelPtr_->dual();
	lastAlgorithm_=2; // dual
      }
    }
  }
  // If scaled feasible but unscaled infeasible take action
  //if (!modelPtr_->status()&&cleanupScaling_) {
  if (cleanupScaling_) {
    modelPtr_->cleanup(cleanupScaling_);
  }
  basis_ = getBasis(modelPtr_);
  //printf("basis after dual\n");
  //basis_.print();
  //modelPtr_->messageHandler()->setLogLevel(saveMessageLevel);
  modelPtr_->popMessageHandler(saveHandler,oldDefault);
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  //modelPtr_->setSolveType(saveSolveType);
  modelPtr_->setSpecialOptions(saveOptions); // restore
}
/* Sets up solver for repeated use by Osi interface.
   The normal usage does things like keeping factorization around so can be used.
   Will also do things like keep scaling and row copy of matrix if
   matrix does not change.
   adventure:
   0 - safe stuff as above
   1 - will take more risks - if it does not work then bug which will be fixed
   2 - don't bother doing most extreme termination checks e.g. don't bother
       re-factorizing if less than 20 iterations.
   3 - Actually safer than 1 (mainly just keeps factorization)

   printOut - -1 always skip round common messages instead of doing some work
               0 skip if normal defaults
               1 leaves
  */
void 
OsiClpSolverInterface::setupForRepeatedUse(int senseOfAdventure, int printOut)
{
  // First try
  switch (senseOfAdventure) {
  case 0:
    specialOptions_=8;
    break;
  case 1:
    specialOptions_=1+2+8;
    break;
  case 2:
    specialOptions_=1+2+4+8;
    break;
  case 3:
    specialOptions_=1+8;
    break;
  }
  bool stopPrinting=false;
  if (printOut<0) {
    stopPrinting=true;
  } else if (!printOut) {
    bool takeHint;
    OsiHintStrength strength;
    bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
    assert (gotHint);
    int messageLevel=messageHandler()->logLevel();
    if (strength!=OsiHintIgnore&&takeHint) 
      messageLevel--;
    stopPrinting = (messageLevel<=0);
  }
  if (stopPrinting) {
    CoinMessages * messagesPointer = modelPtr_->messagesPointer();
    // won't even build messages 
    messagesPointer->setDetailMessages(100,10000,NULL);
  }
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

  const int stat = modelPtr_->status();
  if (stat == 1)
  return true;
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
     if (stat == 0) // optimal
	return maxmin > 0 ? (obj > limit) /*minim*/ : (-obj > limit) /*maxim*/;
     return false;
   case 2: // dual simplex
     if (stat != 0 && stat != 3)
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
  if ((specialOptions_&1024)==0) {
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
  } else {
    // called from CbcNode
    int numberColumns = modelPtr_->numberColumns();
    int numberRows = modelPtr_->numberRows();
    // Get space for crunch and strong branching (too much)
    int size = (1+4*(numberRows+numberColumns))*sizeof(double);
    // and for save of original column bounds
    size += 2*numberColumns*sizeof(double);
    size += (1+4*numberRows+2*numberColumns)*sizeof(int);
    size += numberRows+numberColumns;
    assert (spareArrays_==NULL);
    spareArrays_ = new char[size];
    double * arrayD = (double *) spareArrays_;
    arrayD[0]=modelPtr_->objectiveValue()* modelPtr_->optimizationDirection();
    double * saveSolution = arrayD+1;
    double * saveLower = saveSolution + (numberRows+numberColumns);
    double * saveUpper = saveLower + (numberRows+numberColumns);
    double * saveObjective = saveUpper + (numberRows+numberColumns);
    double * saveLowerOriginal = saveObjective + (numberRows+numberColumns);
    double * saveUpperOriginal = saveLowerOriginal + numberColumns;
    arrayD = saveUpperOriginal + numberColumns;
    int * savePivot = (int *) arrayD;
    int * whichRow = savePivot+numberRows;
    int * whichColumn = whichRow + 3*numberRows;
    int * arrayI = whichColumn + 2*numberColumns;
    //unsigned char * saveStatus = (unsigned char *) (arrayI+1);
    // Use dual region
    double * rhs = modelPtr_->dualRowSolution();
    int nBound;
    ClpSimplex * small = ((ClpSimplexOther *) modelPtr_)->crunch(rhs,whichRow,whichColumn,nBound,true);
    int clpOptions = modelPtr_->specialOptions();
    if((specialOptions_&1)==0) {
      small->setSpecialOptions(clpOptions|(64|1024));
    } else {
      if((specialOptions_&4)==0) 
        small->setSpecialOptions(clpOptions|(64|128|512|1024|4096));
      else
        small->setSpecialOptions(clpOptions|(64|128|512|1024|2048|4096));
    }
    arrayI[0]=nBound;
    assert (small);
    assert (smallModel_==NULL);
    smallModel_=small;
    smallModel_->setLogLevel(0);
    // Setup for strong branching
    assert (factorization_==NULL);
    factorization_ = ((ClpSimplexDual *)smallModel_)->setupForStrongBranching(spareArrays_,numberRows,
                                                          numberColumns);
    int numberColumns2 = smallModel_->numberColumns();
    memcpy( saveLowerOriginal, smallModel_->columnLower(),numberColumns2*sizeof(double));
    memcpy( saveUpperOriginal, smallModel_->columnUpper(),numberColumns2*sizeof(double));
    int itlim;
    modelPtr_->getIntParam(ClpMaxNumIterationHotStart, itlim);
    smallModel_->setIntParam(ClpMaxNumIteration,itlim);
    if (whichRange_&&whichRange_[0]) {
      // get ranging information
      int numberToDo = whichRange_[0];
      int * which = new int [numberToDo];
      // Convert column numbers
      int * backColumn = whichColumn+numberColumns;
      for (int i=0;i<numberToDo;i++) {
        int iColumn = whichRange_[i+1];
        which[i]=backColumn[iColumn];
      }
      double * downRange=new double [numberToDo];
      double * upRange=new double [numberToDo];
      int * whichDown = new int [numberToDo];
      int * whichUp = new int [numberToDo];
      small->gutsOfSolution(NULL,NULL,false);
      ((ClpSimplexOther *) small)->dualRanging(numberToDo,which,
                         upRange, whichUp, downRange, whichDown);
      delete [] whichDown;
      delete [] whichUp;
      delete [] which;
      rowActivity_=upRange;
      columnActivity_=downRange;
    }
  }
}

void OsiClpSolverInterface::solveFromHotStart()
{
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  if (smallModel_==NULL) {
    setWarmStart(ws_);
    memcpy(modelPtr_->primalRowSolution(),
           rowActivity_,numberRows*sizeof(double));
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
  } else {
    double * arrayD = (double *) spareArrays_;
    double saveObjectiveValue = arrayD[0];
    double * saveSolution = arrayD+1;
    int numberRows2 = smallModel_->numberRows();
    int numberColumns2 = smallModel_->numberColumns();
    int number = numberRows2+numberColumns2;
    memcpy(smallModel_->solutionRegion(),saveSolution,number*sizeof(double));
    double * saveLower = saveSolution + (numberRows+numberColumns);
    memcpy(smallModel_->lowerRegion(),saveLower,number*sizeof(double));
    double * saveUpper = saveLower + (numberRows+numberColumns);
    memcpy(smallModel_->upperRegion(),saveUpper,number*sizeof(double));
    double * saveObjective = saveUpper + (numberRows+numberColumns);
    memcpy(smallModel_->costRegion(),saveObjective,number*sizeof(double));
    double * saveLowerOriginal = saveObjective + (numberRows+numberColumns);
    double * saveUpperOriginal = saveLowerOriginal + numberColumns;
    arrayD = saveUpperOriginal + numberColumns;
    int * savePivot = (int *) arrayD;
    memcpy(smallModel_->pivotVariable(),savePivot,numberRows2*sizeof(int));
    int * whichRow = savePivot+numberRows;
    int * whichColumn = whichRow + 3*numberRows;
    int * arrayI = whichColumn + 2*numberColumns;
    unsigned char * saveStatus = (unsigned char *) (arrayI+1);
    memcpy(smallModel_->statusArray(),saveStatus,number);
    smallModel_->setFactorization(*factorization_);
    //int * backColumn = whichColumn+numberColumns;
    const double * lowerBig = modelPtr_->columnLower();
    const double * upperBig = modelPtr_->columnUpper();
    // make sure whatsChanged_ has 1 set
    //smallModel_->setWhatsChanged(1);
    double * lowerSmall = smallModel_->lowerRegion();
    double * upperSmall = smallModel_->upperRegion();
    double * lowerSmallReal = smallModel_->columnLower();
    double * upperSmallReal = smallModel_->columnUpper();
    int i;
    double rhsScale = smallModel_->rhsScale();
    const double * columnScale = (smallModel_->scalingFlag()>0) ? smallModel_->columnScale() : NULL;
    // and do bounds in case dual needs them
    for (i=0;i<numberColumns2;i++) {
      int iColumn = whichColumn[i];
      if (lowerBig[iColumn]>saveLowerOriginal[i]) {
        double value = lowerBig[iColumn];
        lowerSmallReal[i]=value;
        value *= rhsScale;
        if (columnScale)
          value /= columnScale[i];
        lowerSmall[i]=value;
      }
      if (upperBig[iColumn]<saveUpperOriginal[i]) {
        double value = upperBig[iColumn];
        upperSmallReal[i]=value;
        value *= rhsScale;
        if (columnScale)
          value /= columnScale[i];
        upperSmall[i]=value;
      }
    }
    // Start of fast iterations
    bool alwaysFinish= ((specialOptions_&32)==0&&true) ? true : false;
    //smallModel_->setLogLevel(1);
    int status = ((ClpSimplexDual *)smallModel_)->fastDual(alwaysFinish);
    
    int problemStatus = smallModel_->problemStatus();
    double objectiveValue =smallModel_->objectiveValue() * modelPtr_->optimizationDirection();
    CoinAssert (smallModel_->problemStatus()||smallModel_->objectiveValue()<1.0e50);
    // make sure plausible
    double obj = CoinMax(objectiveValue,saveObjectiveValue);
    if (problemStatus==10||problemStatus<0) {
      // was trying to clean up or something odd
      if (problemStatus==10)
        lastAlgorithm_=1; // so won't fail on cutoff (in CbcNode)
      status=1;
    }
    if (status) {
      // not finished - might be optimal
      smallModel_->checkPrimalSolution(smallModel_->solutionRegion(0),
                                     smallModel_->solutionRegion(1));
      //smallModel_->gutsOfSolution(NULL,NULL,0);
      objectiveValue =smallModel_->objectiveValue() *
        modelPtr_->optimizationDirection();
      obj = CoinMax(objectiveValue,saveObjectiveValue);
      if (!smallModel_->numberDualInfeasibilities()) { 
        double limit = 0.0;
        modelPtr_->getDblParam(ClpDualObjectiveLimit, limit);
        if (smallModel_->secondaryStatus()==1&&!problemStatus&&obj<limit) {
#if 0
          // switch off
          ClpSimplex temp = *smallModel_;
          temp.dual();
          if (temp.problemStatus()==0&&temp.objectiveValue()<limit) {
            printf("inf obj %g, true %g - offsets %g %g\n",smallModel_->objectiveValue(),
                   temp.objectiveValue(),
                   smallModel_->objectiveOffset(),temp.objectiveOffset());
          }
          lastAlgorithm_=1;
          obj=limit;
          problemStatus=10;
#else
          obj=limit;
          problemStatus=3;
#endif
        }
        if (!smallModel_->numberPrimalInfeasibilities()&&obj<limit) { 
          problemStatus=0;
#if 0
          ClpSimplex temp = *smallModel_;
          temp.dual();
          if (temp.numberIterations())
            printf("temp iterated\n");
          assert (temp.problemStatus()==0&&temp.objectiveValue()<limit);
#endif
        } else if (problemStatus==10) {
          problemStatus=3;
        } else if (!smallModel_->numberPrimalInfeasibilities()) {
          problemStatus=1; // infeasible
        } 
      } else {
        // can't say much
        lastAlgorithm_=1; // so won't fail on cutoff (in CbcNode)
        problemStatus=3;
      }
    } else if (!problemStatus) {
      if (smallModel_->isDualObjectiveLimitReached()) 
        problemStatus=1; // infeasible
    }
    if (status&&!problemStatus) {
      problemStatus=3; // can't be sure
      lastAlgorithm_=1;
    }
    if (problemStatus<0)
      problemStatus=3;
    modelPtr_->setProblemStatus(problemStatus);
    modelPtr_->setObjectiveValue(obj*modelPtr_->optimizationDirection());
    modelPtr_->setSumDualInfeasibilities(smallModel_->sumDualInfeasibilities());
    modelPtr_->setNumberDualInfeasibilities(smallModel_->numberDualInfeasibilities());
    modelPtr_->setSumPrimalInfeasibilities(smallModel_->sumPrimalInfeasibilities());
    modelPtr_->setNumberPrimalInfeasibilities(smallModel_->numberPrimalInfeasibilities());
    double * solution = modelPtr_->primalColumnSolution();
    const double * solution2 = smallModel_->solutionRegion();
    const double * scale = smallModel_->columnScale();
    if (modelPtr_->scalingFlag()<=0||!scale) {
      for (i=0;i<numberColumns2;i++) {
        int iColumn = whichColumn[i];
        solution[iColumn]= solution2[i];
        lowerSmallReal[i]=saveLowerOriginal[i];
        upperSmallReal[i]=saveUpperOriginal[i];
      }
    } else {
      for (i=0;i<numberColumns2;i++) {
        int iColumn = whichColumn[i];
        solution[iColumn]= solution2[i]*scale[i];
        lowerSmallReal[i]=saveLowerOriginal[i];
        upperSmallReal[i]=saveUpperOriginal[i];
      }
    }
#if 1
    if (status&&!problemStatus) {
      memset(modelPtr_->primalRowSolution(),0,numberRows*sizeof(double));
      modelPtr_->clpMatrix()->times(1.0,solution,modelPtr_->primalRowSolution());
      modelPtr_->checkSolutionInternal();
      //modelPtr_->setLogLevel(1);
      //modelPtr_->allSlackBasis();
      //modelPtr_->primal(1);
      //memset(modelPtr_->primalRowSolution(),0,numberRows*sizeof(double));
      //modelPtr_->clpMatrix()->times(1.0,solution,modelPtr_->primalRowSolution());
      //modelPtr_->checkSolutionInternal();
      assert (!modelPtr_->problemStatus());
    }
#endif
    modelPtr_->setNumberIterations(smallModel_->numberIterations());
    // and back bounds
    memcpy(smallModel_->lowerRegion(),saveLower,number*sizeof(double));
    memcpy(smallModel_->upperRegion(),saveUpper,number*sizeof(double));
  }
}

void OsiClpSolverInterface::unmarkHotStart()
{
  delete [] rowActivity_;
  delete [] columnActivity_;
  rowActivity_=NULL;
  columnActivity_=NULL;
  if (smallModel_==NULL) {
    delete ws_;
    ws_ = NULL;
  } else {
    delete smallModel_;
    delete factorization_;
    delete [] spareArrays_;
    smallModel_=NULL;
    spareArrays_=NULL;
    factorization_=NULL;
  }
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
    matrixByRow_->setExtraGap(0.0);
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
  modelPtr_->setContinuous(index);
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
  modelPtr_->setInteger(index);
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
      modelPtr_->setContinuous(colNumber);
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
    modelPtr_->setInteger(colNumber);
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
/* Add a column (primal variable) to the problem. */
void 
OsiClpSolverInterface::addCol(int numberElements, const int * rows, const double * elements,
			   const double collb, const double colub,   
			   const double obj) 
{
  CoinPackedVector column(numberElements, rows, elements);
  addCol(column,collb,colub,obj);
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
    if (lower[iCol]<-1.0e27)
      lower[iCol]=-COIN_DBL_MAX;
    if (upper[iCol]>1.0e27)
      upper[iCol]=COIN_DBL_MAX;
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
    if (lower[iRow]<-1.0e27)
      lower[iRow]=-COIN_DBL_MAX;
    if (upper[iRow]>1.0e27)
      upper[iRow]=COIN_DBL_MAX;
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
    if (lower[iRow]<-1.0e27)
      lower[iRow]=-COIN_DBL_MAX;
    if (upper[iRow]>1.0e27)
      upper[iRow]=COIN_DBL_MAX;
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
  basis_=CoinWarmStartBasis();
  if (ws_) {
     delete ws_;
     ws_ = 0;
  }
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     double*& rowlb, double*& rowub)
{
  // Get rid of integer information (modelPtr will get rid of its copy)
  loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
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
  assert( rowsen != NULL );
  assert( rowrhs != NULL );
  int numrows = matrix.getNumRows();
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
  }
  loadProblem(matrix, collb, colub, obj, rowlb, rowub);
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
  loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
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
  basis_=CoinWarmStartBasis();
  if (ws_) {
     delete ws_;
     ws_ = 0;
  }
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
  assert( rowsen != NULL );
  assert( rowrhs != NULL );
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
  }
  loadProblem(numcols, numrows, start,  index, value, collb, colub, obj,
	      rowlb,  rowub);
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
smallModel_(NULL),
factorization_(NULL),
spareArrays_(NULL),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL),
cleanupScaling_(0),
specialOptions_(-1)
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
smallModel_(NULL),
factorization_(NULL),
spareArrays_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL)
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
  cleanupScaling_ = rhs.cleanupScaling_;
  specialOptions_ = rhs.specialOptions_;
  fillParamMaps();
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
}

// Borrow constructor - only delete one copy
OsiClpSolverInterface::OsiClpSolverInterface (ClpSimplex * rhs,
					      bool reallyOwn)
:
OsiSolverInterface(),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
rowActivity_(NULL),
columnActivity_(NULL),
smallModel_(NULL),
factorization_(NULL),
spareArrays_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL),
cleanupScaling_(0),
specialOptions_(-1)
{
  modelPtr_ = rhs;
  linearObjective_ = modelPtr_->objective();
  if (rhs) {
    notOwned_=!reallyOwn;

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
  assert(smallModel_==NULL);
  assert(factorization_==NULL);
  assert(spareArrays_==NULL);
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
    saveData_ = rhs.saveData_;
    cleanupScaling_ = rhs.cleanupScaling_;
    specialOptions_ = rhs.specialOptions_;
    basis_ = rhs.basis_;
    if (rhs.integerInformation_) {
      int numberColumns = modelPtr_->numberColumns();
      integerInformation_ = new char[numberColumns];
      memcpy(integerInformation_,rhs.integerInformation_,
	     numberColumns*sizeof(char));
    }
    if ( rhs.ws_ ) 
      ws_ = new CoinWarmStartBasis(*rhs.ws_);
    delete [] rowActivity_;
    delete [] columnActivity_;
    rowActivity_=NULL;
    columnActivity_=NULL;
    assert(smallModel_==NULL);
    assert(factorization_==NULL);
    assert(spareArrays_==NULL);
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
/* Apply a collection of row cuts which are all effective.
   applyCuts seems to do one at a time which seems inefficient.
*/
void 
OsiClpSolverInterface::applyRowCuts(int numberCuts, const OsiRowCut ** cuts)
{
  int i;
  if (!numberCuts)
    return;

  const CoinPackedVectorBase * * rows
    =     new const CoinPackedVectorBase * [numberCuts];
  double * rowlb = new double [numberCuts];
  double * rowub = new double [numberCuts];
  for (i=0;i<numberCuts;i++) {
    rowlb[i] = cuts[i]->lb();
    rowub[i] = cuts[i]->ub();
    rows[i] = &cuts[i]->row();
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
    // Flip slacks
    int lookupA[]={0,1,3,2,0,2};
    for (iRow=0;iRow<numberRows;iRow++) {
      int iStatus = model->getRowStatus(iRow);
      iStatus = lookupA[iStatus];
      basis.setArtifStatus(iRow,(CoinWarmStartBasis::Status) iStatus);
    }
    int lookupS[]={0,1,2,3,0,3};
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int iStatus = model->getColumnStatus(iColumn);
      iStatus = lookupS[iStatus];
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
  // For rows lower and upper are flipped
  for (iRow=0;iRow<numberRows;iRow++) {
    int stat = basis2.getArtifStatus(iRow);
    if (stat>1)
      stat = 5 - stat; // so 2->3 and 3->2
    model->setRowStatus(iRow, (ClpSimplex::Status) stat);
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
  if (numberIntegers) {
    modelPtr_->copyInIntegerInformation(info);
  }
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
  if (modelPtr_->numberIterations()||modelPtr_->upperIn_!=-COIN_DBL_MAX) {
    // This does not pass unitTest when getObjValue is called before solve.
    return modelPtr_->objectiveValue();
  } else {
    return OsiSolverInterface::getObjValue();
  }
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
  modelPtr_->setObjectiveCoefficient(elementIndex,elementValue);
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
  modelPtr_->setColumnLower(elementIndex,elementValue);
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
  modelPtr_->setColumnUpper(elementIndex,elementValue);
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
  modelPtr_->setColumnBounds(elementIndex,lower,upper);
}
void OsiClpSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  const int * indexFirst2=indexFirst;
  while (indexFirst2 != indexLast) {
    const int iColumn=*indexFirst2++;
    if (iColumn<0||iColumn>=n) {
      indexError(iColumn,"setColSetBounds");
    }
  }
#endif
  modelPtr_->setColSetBounds(indexFirst,indexLast,boundList);
}
//------------------------------------------------------------------
/* Set a single row lower bound<br>
   Use -DBL_MAX for -infinity. */
void 
OsiClpSolverInterface::setRowLower( int elementIndex, double elementValue ) {
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowLower");
  }
#endif
  modelPtr_->setRowLower(elementIndex , elementValue);
}
      
/* Set a single row upper bound<br>
   Use DBL_MAX for infinity. */
void 
OsiClpSolverInterface::setRowUpper( int elementIndex, double elementValue ) {
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowUpper");
  }
#endif
  modelPtr_->setRowUpper(elementIndex , elementValue);
}
    
/* Set a single row lower and upper bound */
void 
OsiClpSolverInterface::setRowBounds( int elementIndex,
	      double lower, double upper ) {
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowBounds");
  }
#endif
  modelPtr_->setRowBounds(elementIndex,lower,upper);
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (i<0||i>=n) {
    indexError(i,"setRowType");
  }
#endif
  double lower, upper;
  convertSenseToBound(sense, rightHandSide, range, lower, upper);
  setRowBounds(i, lower, upper);
  // If user is using sense then set
  if (rowsense_) {
    rowsense_[i] = sense;
    rhs_[i] = rightHandSide;
    rowrange_[i] = range;
  }
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  const int * indexFirst2=indexFirst;
  while (indexFirst2 != indexLast) {
    const int iColumn=*indexFirst2++;
    if (iColumn<0||iColumn>=n) {
      indexError(iColumn,"setColumnSetBounds");
    }
  }
#endif
  modelPtr_->setRowSetBounds(indexFirst,indexLast,boundList);
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    double * lower = modelPtr_->rowLower();
    double * upper = modelPtr_->rowUpper();
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
  const int len = indexLast - indexFirst;
  while (indexFirst != indexLast) {
    const int iRow= *indexFirst++;
#ifndef NDEBUG
    if (iRow<0||iRow>=n) {
      indexError(iRow,"isContinuous");
    }
#endif
    double lowerValue;
    double upperValue;
    if (rangeList){
      convertSenseToBound(*senseList++, *rhsList++, *rangeList++,
			  lowerValue, upperValue);
    } else {
      convertSenseToBound(*senseList++, *rhsList++, 0,
			  lowerValue, upperValue);
    }
    modelPtr_->setRowBounds(iRow,lowerValue,upperValue);
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
  // Do initialization
  saveData_ = modelPtr_->saveData();
  saveData_.scalingFlag_=modelPtr_->scalingFlag();
  modelPtr_->scaling(0);
  specialOptions_ = -1;
  // set infeasibility cost up
  modelPtr_->setInfeasibilityCost(1.0e12);
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
  // message will not appear anyway
  int saveMessageLevel=modelPtr_->messageHandler()->logLevel();
  modelPtr_->messageHandler()->setLogLevel(0);
  modelPtr_->finish();
  modelPtr_->messageHandler()->setLogLevel(saveMessageLevel);
  modelPtr_->restoreData(saveData_);
  modelPtr_->scaling(saveData_.scalingFlag_);
  ClpDualRowSteepest steepest;
  modelPtr_->setDualRowPivotAlgorithm(steepest);
  ClpPrimalColumnSteepest steepestP;
  modelPtr_->setPrimalColumnPivotAlgorithm(steepestP);
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
  // Flip slacks
  int lookupA[]={0,1,3,2,0,2};
  for (i=0;i<n;i++) 
    rstat[i] = lookupA[modelPtr_->getRowStatus(i)];
  n=modelPtr_->numberColumns();
  int lookupS[]={0,1,2,3,0,3};
  for (i=0;i<n;i++)
    cstat[i] = lookupS[modelPtr_->getColumnStatus(i)];
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
  // For rows lower and upper are flipped
  int lookupA[]={0,1,3,2};
  for (i=0;i<n;i++) {
    int status = lookupA[rstat[i]];
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
  // Save 
  basis_ = getBasis(modelPtr_);
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
  if (modelPtr_->nonLinearCost()) {
    modelPtr_->nonLinearCost()->refreshCosts(c);
  }
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
  assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  CoinIndexedVector * columnArray0 = modelPtr_->columnArray(0);
  CoinIndexedVector * columnArray1 = modelPtr_->columnArray(1);
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  // put +1 in row 
  // But swap if pivot variable was slack as clp stores slack as -1.0
  const int * pivotVariable = modelPtr_->pivotVariable();
  const double * rowScale = modelPtr_->rowScale();
  const double * columnScale = modelPtr_->columnScale();
  int pivot = pivotVariable[row];
  double value;
  // And if scaled then adjust
  if (!rowScale) {
    if (pivot<numberColumns)
      value = 1.0;
    else
      value = -1.0;
  } else {
    if (pivot<numberColumns)
      value = columnScale[pivot];
    else
      value = -1.0/rowScale[pivot-numberColumns];
  }
  rowArray1->insert(row,value);
  modelPtr_->factorization()->updateColumnTranspose(rowArray0,rowArray1);
  // put row of tableau in rowArray1 and columnArray0
  modelPtr_->clpMatrix()->transposeTimes(modelPtr_,1.0,
                                         rowArray1,columnArray1,columnArray0);
  if (!rowScale) {
    memcpy(z,columnArray0->denseVector(),
           numberColumns*sizeof(double));
  } else {
    double * array = columnArray0->denseVector();
    for (int i=0;i<numberColumns;i++)
      z[i] = array[i]/columnScale[i];
  }
  if (slack) {
    if (!rowScale) {
      memcpy(slack,rowArray1->denseVector(),
             numberRows*sizeof(double));
    } else {
      double * array = rowArray1->denseVector();
      for (int i=0;i<numberRows;i++)
        slack[i] = array[i]*rowScale[i];
    }
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
  assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // put +1 in row
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double value = (modelPtr_->pivotVariable()[row]<modelPtr_->numberColumns()) ? 1.0 : -1.0;
  // What about scaling ?
  rowArray1->insert(row,value);
  factorization->updateColumnTranspose(rowArray0,rowArray1);
  memcpy(z,rowArray1->denseVector(),modelPtr_->numberRows()*sizeof(double));
  rowArray1->clear();
}

//Get a column of the tableau
void 
OsiClpSolverInterface::getBInvACol(int col, double* vec)
{
  assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // get column of matrix
#ifndef NDEBUG
  int n = modelPtr_->numberColumns()+modelPtr_->numberRows();
  if (col<0||col>=n) {
    indexError(col,"getBInvACol");
  }
#endif
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const int * pivotVariable = modelPtr_->pivotVariable();
  const double * rowScale = modelPtr_->rowScale();
  const double * columnScale = modelPtr_->columnScale();
  if (!rowScale) {
    if (col<numberColumns) {
      modelPtr_->unpack(rowArray1,col);
    } else {
      rowArray1->insert(col-numberColumns,1.0);
    }
  } else {
    if (col<numberColumns) {
      modelPtr_->unpack(rowArray1,col);
      double multiplier = 1.0/columnScale[col];
      int number = rowArray1->getNumElements();
      int * index = rowArray1->getIndices();
      double * array = rowArray1->denseVector();
      for (int i=0;i<number;i++) {
	int iRow = index[i];
	// make sure not packed
	assert (array[iRow]);
	array[iRow] *= multiplier;
      }
    } else {
      rowArray1->insert(col-numberColumns,rowScale[col-numberColumns]);
    }
  }
  modelPtr_->factorization()->updateColumn(rowArray0,rowArray1,false);
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double * array = rowArray1->denseVector();
  if (!rowScale) {
    for (int i=0;i<numberRows;i++) {
      double multiplier = (pivotVariable[i]<numberColumns) ? 1.0 : -1.0;
      vec[i] = multiplier * array[i];
    }
  } else {
    for (int i=0;i<numberRows;i++) {
      int pivot = pivotVariable[i];
      if (pivot<numberColumns)
	vec[i] = array[i] * columnScale[pivot];
      else
	vec[i] = - array[i] / rowScale[pivot-numberColumns];
    }
  }
  rowArray1->clear();
}

//Get a column of the basis inverse
void 
OsiClpSolverInterface::getBInvCol(int col, double* vec)
{
  assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
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
  assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  assert (index);
  assert (modelPtr_->pivotVariable());
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
  assert(smallModel_==NULL);
  assert(factorization_==NULL);
  assert(spareArrays_==NULL);
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
// Set a hint parameter
bool 
OsiClpSolverInterface::setHintParam(OsiHintParam key, bool yesNo,
                                    OsiHintStrength strength,
                                    void * otherInformation)
{
  if ( OsiSolverInterface::setHintParam(key,yesNo,strength,otherInformation)) {
    // special coding for branch and cut
    if (yesNo&&strength == OsiHintDo&&key==OsiDoInBranchAndCut) {
      if ( specialOptions_==-1) 
        setupForRepeatedUse(0,0);
      // set normal
      specialOptions_ &= 1023;
      if (otherInformation!=NULL) {
        int * array = (int *) otherInformation;
        if (array[0]>=0||array[0]<=2)
          specialOptions_ |= array[0]<<10;
      }
    }
    return true;
  } else {
    return false;
  }
}


// Invoke solver's built-in enumeration algorithm
void 
OsiClpSolverInterface::branchAndBound() {
  throw CoinError("Sorry, Clp doesn't implement B&B.  This is because the original Simple Branch and Bound\n\
code has been much developed and moved to Sbb which is now Cbc\n\
.  If you want a native Branch and Bound code use Cbc",
		  "branchAndBound","OsiClpSolverInterface");
}
// Crunch down model
void 
OsiClpSolverInterface::crunch()
{
  int numberColumns = modelPtr_->numberColumns();
  int numberRows = modelPtr_->numberRows();
  // Use dual region
  double * rhs = modelPtr_->dualRowSolution();
  int * whichRow = new int[3*numberRows];
  int * whichColumn = new int[2*numberColumns];
  int nBound;
  ClpSimplex * small = ((ClpSimplexOther *) modelPtr_)->crunch(rhs,whichRow,whichColumn,nBound);
  if (small) {
    small->dual();
    if (small->problemStatus()==0) {
      modelPtr_->setProblemStatus(0);
      ((ClpSimplexOther *) modelPtr_)->afterCrunch(*small,whichRow,whichColumn,nBound);
    } else {
      modelPtr_->setProblemStatus(1);
    }
    delete small;
  } else {
    modelPtr_->setProblemStatus(1);
  }
  delete [] whichRow;
  delete [] whichColumn;
}
