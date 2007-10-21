// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>
#ifndef CLP_FAST_CODE
// If NDEBUG defined then Osi unit test fails but code will be slightly faster
#ifdef NDEBUG
#undef NDEBUG
#endif
#endif
 
#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinModel.hpp"
#include "CoinMpsIO.hpp"
#include "CoinSort.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSimplexPrimal.hpp"
#include "ClpSimplexDual.hpp"
#include "ClpNonLinearCost.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiBranchingObject.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "ClpPresolve.hpp"
#include "CoinLpIO.hpp"
static double totalTime=0.0;
//#define SAVE_MODEL 1
#ifdef SAVE_MODEL
static int resolveTry=0;
static int loResolveTry=0;
static int hiResolveTry=9999999;
#endif
//#############################################################################
// Solve methods
//#############################################################################
void OsiClpSolverInterface::initialSolve()
{
  ClpSimplex solver(true);
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
  // get original log levels
  int saveMessageLevel=modelPtr_->logLevel();
  int messageLevel=messageHandler()->logLevel();
  // Set message handler
  solver.passInMessageHandler(handler_);
  // But keep log level
  solver.messageHandler()->setLogLevel(saveMessageLevel);
  // See if user set factorization frequency
  int userFactorizationFrequency = modelPtr_->factorization()->maximumPivots();
  // borrowModel does not move
  solver.factorization()->maximumPivots(userFactorizationFrequency);
  // set reasonable defaults
  bool takeHint;
  OsiHintStrength strength;
  // Switch off printing if asked to
  bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    if (messageLevel>0)
      messageLevel--;
  }
  if (messageLevel<saveMessageLevel)
    solver.messageHandler()->setLogLevel(messageLevel);
  // Allow for specialOptions_==1+8 forcing saving factorization
  int startFinishOptions=0;
  if ((specialOptions_&9)==(1+8)) {
    startFinishOptions =1+2+4; // allow re-use of factorization
  }
  bool defaultHints=true;
  {
    int hint;
    for (hint=OsiDoPresolveInInitial;hint<OsiLastHintParam;hint++) {
      if (hint!=OsiDoReducePrint&&
          hint!=OsiDoInBranchAndCut) {
        bool yesNo;
        OsiHintStrength strength;
        getHintParam((OsiHintParam) hint,yesNo,strength);
        if (yesNo) {
          defaultHints=false;
          break;
        }
        if (strength != OsiHintIgnore) {
          defaultHints=false;
          break;
        }
      }
    }
  }
  /*
    If basis then do primal (as user could do dual with resolve)
    If not then see if dual feasible (and allow for gubs etc?)
  */
  bool doPrimal = (basis_.numberBasicStructurals()>0);
  setBasis(basis_,&solver);
  if ((!defaultHints||doPrimal)&&!solveOptions_.getSpecialOption(6)) {
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
      } else {
	model2->setSpecialOptions(solver.specialOptions());
      }
      
      // change from 200 (unless changed)
      if (modelPtr_->factorization()->maximumPivots()==200)
        model2->factorization()->maximumPivots(100+model2->numberRows()/50);
      else
        model2->factorization()->maximumPivots(userFactorizationFrequency);
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
        bool stopped = model2->status()==3;
        pinfo.postsolve(true);
        
        delete model2;
        //printf("Resolving from postsolved model\n");
        // later try without (1) and check duals before solve
	if (!stopped)
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
      OsiClpDisasterHandler handler(this);
      handler.setSimplex(&solver); // as "borrowed"
      bool inCbcOrOther = (modelPtr_->specialOptions()&0x03000000)!=0;
      if (!doPrimal) 
	handler.setWhereFrom(4);
      else
	handler.setWhereFrom(6);
      if (inCbcOrOther)
	solver.setDisasterHandler(&handler);
      if (!doPrimal) {
        //printf("doing dual\n");
        solver.dual(0);
	if (inCbcOrOther) {
	  if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	    printf("dual trouble a\n");
#endif
	    // try just going back in
	    handler.setPhase(1);
	    solver.dual();
	    if (handler.inTrouble()) {
#ifdef COIN_DEVELOP
	      printf("dual trouble b\n");
#endif
	      // try primal with original basis
	      handler.setPhase(2);
	      setBasis(basis_,&solver);
	      solver.primal();
	    }
	    if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	      printf("disaster - treat as infeasible\n");
#endif
	      solver.setProblemStatus(1);
	    }
	  }
	  // reset
	  solver.setDisasterHandler(NULL);
	}
        lastAlgorithm_=2; // dual
        // check if clp thought it was in a loop
        if (solver.status()==3&&!solver.hitMaximumIterations()) {
          // switch algorithm
          solver.primal(0);
          lastAlgorithm_=1; // primal
        }
      } else {
        //printf("doing primal\n");
        solver.primal(1);
	if (inCbcOrOther) {
	  if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	    printf("primal trouble a\n");
#endif
	    // try just going back in (but with dual)
	    handler.setPhase(1);
	    solver.dual();
	    if (handler.inTrouble()) {
#ifdef COIN_DEVELOP
	      printf("primal trouble b\n");
#endif
	      // try primal with original basis
	      handler.setPhase(2);
	      setBasis(basis_,&solver);
	      solver.dual();
	    }
	    if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	      printf("disaster - treat as infeasible\n");
#endif
	      solver.setProblemStatus(1);
	    }
	  }
	  // reset
	  solver.setDisasterHandler(NULL);
	}
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
    const double * rowScale2 = solver.rowScale();
    solver.setSpecialOptions(saveOptions);
    if (!rowScale1&&rowScale2) {
      // need to release memory
      solver.setRowScale(NULL);
      solver.setColumnScale(NULL);
    }
  } else {
    // User doing nothing and all slack basis
    ClpSolve options=solveOptions_;
    bool yesNo;
    OsiHintStrength strength;
    getHintParam(OsiDoInBranchAndCut,yesNo,strength);
    if (yesNo) {
      solver.setSpecialOptions(solver.specialOptions()|1024);
    }
    solver.initialSolve(options);
    lastAlgorithm_ = 2; // say dual
    // If scaled feasible but unscaled infeasible take action
    if (!solver.status()&&cleanupScaling_) {
      solver.cleanup(cleanupScaling_);
    }
    basis_ = getBasis(&solver);
    //basis_.print();
  }
  solver.messageHandler()->setLogLevel(saveMessageLevel);
  solver.returnModel(*modelPtr_);
  if (startFinishOptions) {
    int save = modelPtr_->logLevel();
    if (save<2) modelPtr_->setLogLevel(0);
    modelPtr_->dual(0,startFinishOptions);
    modelPtr_->setLogLevel(save);
  }
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  if (modelPtr_->problemStatus_==3&&lastAlgorithm_==2)
    modelPtr_->computeObjectiveValue();
  // mark so we can pick up objective value quickly
  modelPtr_->upperIn_=0.0;
  time1 = CoinCpuTime()-time1;
  totalTime += time1;
  if (lastAlgorithm_<1||lastAlgorithm_>2)
    lastAlgorithm_=1;
  //std::cout<<time1<<" seconds - total "<<totalTime<<std::endl;
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::resolve()
{
#ifdef COIN_DEVELOP
  {
    int i;
    int n = getNumCols();
    const double *lower = getColLower() ;
    const double *upper = getColUpper() ;
    for (i=0;i<n;i++) {
      assert (lower[i]<1.0e10);
      assert (upper[i]>-1.0e10);
    }
    n = getNumRows();
    lower = getRowLower() ;
    upper = getRowUpper() ;
    for (i=0;i<n;i++) {
      assert (lower[i]<1.0e10);
      assert (upper[i]>-1.0e10);
    }
  }
#endif
  //void pclp(char *);
  //pclp("res");
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = (getHintParam(OsiDoInBranchAndCut,takeHint,strength));
  assert (gotHint);
  // mark so we can pick up objective value quickly
  modelPtr_->upperIn_=0.0;
  if (((specialOptions_>>10)&2)!=0) {
    // Quick check to see if optimal
    modelPtr_->checkSolutionInternal();
    if (modelPtr_->problemStatus()==0) {
      modelPtr_->setNumberIterations(0);
      return;
    }
  }
  // If using Clp initialSolve and primal - just do here
  gotHint = (getHintParam(OsiDoDualInResolve,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&!takeHint&&solveOptions_.getSpecialOption(6)) {
    ClpSolve options=solveOptions_;
    // presolve
    getHintParam(OsiDoPresolveInResolve,takeHint,strength);
    if (strength!=OsiHintIgnore&&!takeHint)
      options.setPresolveType(ClpSolve::presolveOff);
    int saveOptions = modelPtr_->specialOptions();
    getHintParam(OsiDoInBranchAndCut,takeHint,strength);
    if (takeHint) {
      modelPtr_->setSpecialOptions(modelPtr_->specialOptions()|1024);
    }
    setBasis(basis_,modelPtr_);
    modelPtr_->initialSolve(options);
    lastAlgorithm_ = 1; // say primal
    // If scaled feasible but unscaled infeasible take action
    if (!modelPtr_->status()&&cleanupScaling_) {
      modelPtr_->cleanup(cleanupScaling_);
    }
    modelPtr_->setSpecialOptions(saveOptions); // restore
    basis_ = getBasis(modelPtr_);
  }
  int saveSolveType=modelPtr_->solveType();
  bool doingPrimal = modelPtr_->algorithm()>0;
  if (saveSolveType==2) {
    disableSimplexInterface();
  }
  int saveOptions = modelPtr_->specialOptions();
  int startFinishOptions=0;
  if (specialOptions_!=0x80000000) {
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
  int saveMessageLevel=modelPtr_->logLevel();
  int messageLevel=messageHandler()->logLevel();
  bool oldDefault;
  CoinMessageHandler * saveHandler = modelPtr_->pushMessageHandler(handler_,oldDefault);
  //printf("basis before dual\n");
  //basis_.print();
  setBasis(basis_,modelPtr_);
#ifdef SAVE_MODEL
  resolveTry++;
#if SAVE_MODEL > 1
  if (resolveTry>=loResolveTry&&
      resolveTry<=hiResolveTry) {
    char fileName[20];
    sprintf(fileName,"save%d.mod",resolveTry);
    modelPtr_->saveModel(fileName);
  }
#endif
#endif
  // set reasonable defaults
  // Switch off printing if asked to
  gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    if (messageLevel>0)
      messageLevel--;
  }
  if (messageLevel<modelPtr_->messageHandler()->logLevel())
    modelPtr_->messageHandler()->setLogLevel(messageLevel);
  // See if user set factorization frequency
  int userFactorizationFrequency = modelPtr_->factorization()->maximumPivots();
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
    if ((specialOptions_&128)!=0) {
      specialOptions_ &= ~128;
      modelPtr_->deleteAuxiliaryModel();
    }
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
    if (modelPtr_->factorization()->maximumPivots()==200)
      model2->factorization()->maximumPivots(100+model2->numberRows()/50);
    else
      model2->factorization()->maximumPivots(userFactorizationFrequency);
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
      //writeMps("try1");
      int savePerturbation = modelPtr_->perturbation();
      if ((specialOptions_&2)!=0)
	modelPtr_->setPerturbation(100);
      //modelPtr_->messageHandler()->setLogLevel(1);
      //writeMpsNative("bad",NULL,NULL,2,1,1.0);
      OsiClpDisasterHandler handler(this);
      bool inCbcOrOther = (modelPtr_->specialOptions()&0x03000000)!=0;
      if (((modelPtr_->specialOptions()&1024)==0||(specialOptions_ &128)!=0)&&
          modelPtr_->auxiliaryModel_) {
        if ((specialOptions_&128)==0) {
	  handler.setWhereFrom(0); // dual
	  if (inCbcOrOther)
	    modelPtr_->setDisasterHandler(&handler);
	  bool specialScale;
	  if ((specialOptions_&131072)!=0&&!modelPtr_->rowScale_) {
	    modelPtr_->rowScale_ = rowScale_.array();
	    modelPtr_->columnScale_ = columnScale_.array();
	    specialScale=true;
	  } else {
	    specialScale=false;
	  }
          modelPtr_->dual(0,startFinishOptions);
	  if (specialScale) {
	    modelPtr_->rowScale_ = NULL;
	    modelPtr_->columnScale_ = NULL;
	  }
        } else {
          double * rhs = modelPtr_->auxiliaryModel_->lower_;
          int numberTightened = ((ClpSimplexOther *)modelPtr_)->tightenIntegerBounds(rhs);
          if (numberTightened>=0) {
	    handler.setWhereFrom(0); // dual
	    if (inCbcOrOther)
	      modelPtr_->setDisasterHandler(&handler);
	    bool specialScale;
	    if ((specialOptions_&131072)!=0&&!modelPtr_->rowScale_) {
	      modelPtr_->rowScale_ = rowScale_.array();
	      modelPtr_->columnScale_ = columnScale_.array();
	      specialScale=true;
	    } else {
	      specialScale=false;
	    }
	    modelPtr_->dual(0,0);
	    if (specialScale) {
	      modelPtr_->rowScale_ = NULL;
	      modelPtr_->columnScale_ = NULL;
	    }
          } else {
            modelPtr_->setProblemStatus(1);
	  }
        }
      } else {
 	if((specialOptions_&1)==0) {
	  handler.setWhereFrom(0); // dual
	  if (inCbcOrOther)
	    modelPtr_->setDisasterHandler(&handler);
	  bool specialScale;
	  if ((specialOptions_&131072)!=0&&!modelPtr_->rowScale_) {
	    modelPtr_->rowScale_ = rowScale_.array();
	    modelPtr_->columnScale_ = columnScale_.array();
	    specialScale=true;
	  } else {
	    specialScale=false;
	  }
          modelPtr_->dual(0,startFinishOptions);
	  if (specialScale) {
	    modelPtr_->rowScale_ = NULL;
	    modelPtr_->columnScale_ = NULL;
	  }
        } else {
          crunch();
	  // should have already been fixed if problems
	  inCbcOrOther=false;
        }
      }
      if (inCbcOrOther) {
	if(handler.inTrouble()) {
	  // try just going back in
	  handler.setPhase(1);
	  modelPtr_->dual();
	  if (handler.inTrouble()) {
	    // try primal with original basis
	    handler.setPhase(2);
	    setBasis(basis_,modelPtr_);
	    modelPtr_->primal();
	  }
	  if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	    printf("disaster - treat as infeasible\n");
#endif
	    modelPtr_->setProblemStatus(1);
	  }
	}
	// reset
	modelPtr_->setDisasterHandler(NULL);
      }
      if (modelPtr_->problemStatus()==4) {
	// bad bounds?
	modelPtr_->setProblemStatus(1);
      }
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
#ifdef COIN_DEVELOP
	  printf("in trouble - try all slack\n");
#endif
	  CoinWarmStartBasis allSlack;
	  setBasis(allSlack,modelPtr_);
	  modelPtr_->dual();
          if (modelPtr_->status()==3&&!modelPtr_->hitMaximumIterations()) {
	    if (modelPtr_->numberPrimalInfeasibilities()) {
#ifdef COIN_DEVELOP
	      printf("Real real trouble - treat as infeasible\n");
#endif
	      modelPtr_->setProblemStatus(1);
	    } else {
#ifdef COIN_DEVELOP
	      printf("Real real trouble - treat as optimal\n");
#endif
	      modelPtr_->setProblemStatus(0);
	    }
	  }
	}
      }
      assert (modelPtr_->objectiveValue()<1.0e100);
    } else {
      //printf("doing primal\n");
      modelPtr_->primal(1,startFinishOptions);
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
  modelPtr_->popMessageHandler(saveHandler,oldDefault);
  modelPtr_->messageHandler()->setLogLevel(saveMessageLevel);
  if (saveSolveType==2) {
    enableSimplexInterface(doingPrimal);
  }
  //modelPtr_->setSolveType(saveSolveType);
  modelPtr_->setSpecialOptions(saveOptions); // restore
  if (modelPtr_->problemStatus_==3&&lastAlgorithm_==2)
    modelPtr_->computeObjectiveValue();
  if (lastAlgorithm_<1||lastAlgorithm_>2)
    lastAlgorithm_=1;
#ifdef SAVE_MODEL
  if (resolveTry>=loResolveTry&&
      resolveTry<=hiResolveTry) {
    printf("resolve %d took %d iterations - algorithm %d\n",resolveTry,modelPtr_->numberIterations(),lastAlgorithm_);
  }
#endif
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
    getHintParam(OsiDoReducePrint,takeHint,strength);
    int messageLevel=messageHandler()->logLevel();
    if (strength!=OsiHintIgnore&&takeHint) 
      messageLevel--;
    stopPrinting = (messageLevel<=0);
  }
  if (stopPrinting) {
    CoinMessages * messagesPointer = modelPtr_->messagesPointer();
    // won't even build messages 
    messagesPointer->setDetailMessages(100,10000,(int *) NULL);
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
  return modelPtr_->setIntParam((ClpIntParam) key, value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setDblParam(OsiDblParam key, double value)
{
  if (key != OsiLastDblParam ) {
    if (key==OsiDualObjectiveLimit||key==OsiPrimalObjectiveLimit)
      value *= modelPtr_->optimizationDirection();
    return modelPtr_->setDblParam((ClpDblParam) key, value);
  } else {
    return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  assert (key!=OsiSolverName);
  if (key != OsiLastStrParam ) {
    return modelPtr_->setStrParam((ClpStrParam) key, value);
  } else {
    return false;
  }
}


//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getIntParam(OsiIntParam key, int& value) const 
{
  return modelPtr_->getIntParam((ClpIntParam) key, value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  if (key != OsiLastDblParam ) {
    bool condition =  modelPtr_->getDblParam((ClpDblParam) key, value);
    if (key==OsiDualObjectiveLimit||key==OsiPrimalObjectiveLimit)
      value *= modelPtr_->optimizationDirection();
    return condition;
  } else {
    return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  if ( key==OsiSolverName ) {
    value = "clp";
    return true;
  }
  if (key != OsiLastStrParam ) {
    return modelPtr_->getStrParam((ClpStrParam) key, value);
  } else {
    return false;
  }
}


//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiClpSolverInterface::isAbandoned() const
{
  // not sure about -1 (should not happen)
  return (modelPtr_->status()==4||modelPtr_->status()==-1||
	  (modelPtr_->status()==1&&modelPtr_->secondaryStatus()==8));
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
  if (fabs(limit) > 1e30) {
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
  if (fabs(limit) > 1e30) {
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
  modelPtr_->setProblemStatus(0);
#define CLEAN_HOT_START
#ifdef CLEAN_HOT_START
  if ((specialOptions_&65536)!=0) {
    assert ((specialOptions_&128)==0);
    // space for save arrays
    int numberColumns = modelPtr_->numberColumns();
    int numberRows = modelPtr_->numberRows();
    // Get space for strong branching 
    int size = (1+4*(numberRows+numberColumns))*sizeof(double);
    // and for save of original column bounds
    size += 2*numberColumns*sizeof(double);
    size += (1+4*numberRows+2*numberColumns)*sizeof(int);
    size += numberRows+numberColumns;
    assert (spareArrays_==NULL);
    spareArrays_ = new char[size];
    // Setup for strong branching
    assert (factorization_==NULL);
    if ((specialOptions_&131072)!=0) {
      assert (lastNumberRows_>=0);
      if (modelPtr_->rowScale_!=rowScale_.array()) {
	assert(modelPtr_->columnScale_!=columnScale_.array());
	delete [] modelPtr_->rowScale_;
	modelPtr_->rowScale_=NULL;
	delete [] modelPtr_->columnScale_;
	modelPtr_->columnScale_=NULL;
	if (lastNumberRows_==modelPtr_->numberRows()) {
	  // use scaling
	  modelPtr_->rowScale_ = rowScale_.array();
	  modelPtr_->columnScale_ = columnScale_.array();
	} else {
	  specialOptions_ &= ~131072;
	  modelPtr_->specialOptions_ &= ~131072;
	}
      }
      lastNumberRows_ = -1 -lastNumberRows_;
    }
    factorization_ = ((ClpSimplexDual *)modelPtr_)->setupForStrongBranching(spareArrays_,numberRows,
									    numberColumns);
    double * arrayD = (double *) spareArrays_;
    arrayD[0]=modelPtr_->objectiveValue()* modelPtr_->optimizationDirection();
    double * saveSolution = arrayD+1;
    double * saveLower = saveSolution + (numberRows+numberColumns);
    double * saveUpper = saveLower + (numberRows+numberColumns);
    double * saveObjective = saveUpper + (numberRows+numberColumns);
    double * saveLowerOriginal = saveObjective + (numberRows+numberColumns);
    double * saveUpperOriginal = saveLowerOriginal + numberColumns;
    memcpy( saveLowerOriginal, modelPtr_->columnLower(),numberColumns*sizeof(double));
    memcpy( saveUpperOriginal, modelPtr_->columnUpper(),numberColumns*sizeof(double));
#if 0
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
      modelPtr_->gutsOfSolution(NULL,NULL,false);
      // Tell code we can increase costs in some cases
      modelPtr_->setCurrentDualTolerance(0.0);
      ((ClpSimplexOther *) modelPtr_)->dualRanging(numberToDo,which,
						   upRange, whichUp, downRange, whichDown);
      delete [] whichDown;
      delete [] whichUp;
      delete [] which;
      rowActivity_=upRange;
      columnActivity_=downRange;
    }
#endif
    return;
  }
#endif
  if ((specialOptions_&8192)==0) { // ||(specialOptions_&65536)!=0) {
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
#if 0
    int saveLevel = modelPtr_->logLevel();
    modelPtr_->setLogLevel(0);
    //modelPtr_->dual();
    OsiClpSolverInterface::resolve();
    if (modelPtr_->numberIterations()>0)
      printf("**** iterated large %d\n",modelPtr_->numberIterations());
    //else
    //printf("no iterations\n");
    modelPtr_->setLogLevel(saveLevel);
#endif
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
    int nBound=0;
    bool keepModel = modelPtr_->auxiliaryModel_!=NULL;
    ClpSimplex * small;
    if (!keepModel) {
      small = ((ClpSimplexOther *) modelPtr_)->crunch(rhs,whichRow,whichColumn,nBound,true);
      if (small&&(specialOptions_&131072)!=0) {
	assert (lastNumberRows_>=0);
	int numberRows2 = small->numberRows();
	int numberColumns2 = small->numberColumns();
	double * rowScale2 = new double [2*numberRows2];
	const double * rowScale = rowScale_.array();
	double * inverseScale2 = rowScale2+numberRows2;
	const double * inverseScale = rowScale+modelPtr_->numberRows_;
	int i;
	for (i=0;i<numberRows2;i++) {
	  int iRow = whichRow[i];
	  rowScale2[i]=rowScale[iRow];
	  inverseScale2[i]=inverseScale[iRow];
	}
	small->setRowScale(rowScale2);
	double * columnScale2 = new double [2*numberColumns2];
	const double * columnScale = columnScale_.array();
	inverseScale2 = columnScale2+numberColumns2;
	inverseScale = columnScale+modelPtr_->numberColumns_;
	for (i=0;i<numberColumns2;i++) {
	  int iColumn = whichColumn[i];
	  columnScale2[i]=columnScale[iColumn];
	  inverseScale2[i]=inverseScale[iColumn];
	}
	small->setColumnScale(columnScale2);
	small->specialOptions_ |= 131072;
      }
    } else {
      // save stuff
      small=modelPtr_;
      modelPtr_->auxiliaryModel_->numberPrimalInfeasibilities_=modelPtr_->logLevel();
      // make sure auxiliary model won't get deleted
      //modelPtr_->whatsChanged_ |= 511;
      int itlim;
      modelPtr_->getIntParam(ClpMaxNumIteration, itlim);
      modelPtr_->auxiliaryModel_->numberDualInfeasibilities_=itlim;
      CoinIotaN(whichRow,numberRows,0);
      CoinIotaN(whichColumn,numberColumns,0);
      CoinIotaN(whichColumn+numberColumns,numberColumns,0);
    }
    if (!small) {
      // should never be infeasible .... but
      delete [] spareArrays_;
      spareArrays_=NULL;
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
      modelPtr_->setProblemStatus(1);
      return;
    }
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
    assert (smallModel_==NULL);
    if ((specialOptions_&256)!=0||1) {
      // only need to do this on second pass in CbcNode
      if (modelPtr_->logLevel()<2) small->setLogLevel(0);
      small->dual();
      if (keepModel&&!small->auxiliaryModel_) {
        // put back
        synchronizeModel();
	// make sure auxiliary model won't get deleted
	modelPtr_->whatsChanged_ |= 511;
      }
      if (small->numberIterations()>0&&small->logLevel()>2)
	printf("**** iterated small %d\n",small->numberIterations());
      //small->setLogLevel(0);
      // Could be infeasible if forced one way (and other way stopped on iterations)
      if (small->status()==1) {
	if (small!=modelPtr_)
	  delete small;
	delete [] spareArrays_;
	spareArrays_=NULL;
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
	modelPtr_->setProblemStatus(1);
	return;
      }
    }
    smallModel_=small;
    if (modelPtr_->logLevel()<2) smallModel_->setLogLevel(0);
    // Setup for strong branching
    assert (factorization_==NULL);
    factorization_ = ((ClpSimplexDual *)smallModel_)->setupForStrongBranching(spareArrays_,numberRows,
                                                          numberColumns);
    int numberColumns2 = smallModel_->numberColumns();
    memcpy( saveLowerOriginal, smallModel_->columnLower(),numberColumns2*sizeof(double));
    memcpy( saveUpperOriginal, smallModel_->columnUpper(),numberColumns2*sizeof(double));
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
      // Tell code we can increase costs in some cases
      small->setCurrentDualTolerance(0.0);
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
  modelPtr_->getIntParam(ClpMaxNumIteration,itlimOrig_);
  int itlim;
  modelPtr_->getIntParam(ClpMaxNumIterationHotStart, itlim);
#ifdef CLEAN_HOT_START
  if ((specialOptions_&65536)!=0) {
    double * arrayD = (double *) spareArrays_;
    double saveObjectiveValue = arrayD[0];
    double * saveSolution = arrayD+1;
    int number = numberRows+numberColumns;
    memcpy(modelPtr_->solutionRegion(),saveSolution,number*sizeof(double));
    double * saveLower = saveSolution + (numberRows+numberColumns);
    memcpy(modelPtr_->lowerRegion(),saveLower,number*sizeof(double));
    double * saveUpper = saveLower + (numberRows+numberColumns);
    memcpy(modelPtr_->upperRegion(),saveUpper,number*sizeof(double));
    double * saveObjective = saveUpper + (numberRows+numberColumns);
    memcpy(modelPtr_->costRegion(),saveObjective,number*sizeof(double));
    double * saveLowerOriginal = saveObjective + (numberRows+numberColumns);
    double * saveUpperOriginal = saveLowerOriginal + numberColumns;
    arrayD = saveUpperOriginal + numberColumns;
    int * savePivot = (int *) arrayD;
    memcpy(modelPtr_->pivotVariable(),savePivot,numberRows*sizeof(int));
    int * whichRow = savePivot+numberRows;
    int * whichColumn = whichRow + 3*numberRows;
    int * arrayI = whichColumn + 2*numberColumns;
    unsigned char * saveStatus = (unsigned char *) (arrayI+1);
    memcpy(modelPtr_->statusArray(),saveStatus,number);
    modelPtr_->setFactorization(*factorization_);
    double * columnLower = modelPtr_->columnLower();
    double * columnUpper = modelPtr_->columnUpper();
    // make sure whatsChanged_ has 1 set
    modelPtr_->setWhatsChanged(511);
    double * lowerInternal = modelPtr_->lowerRegion();
    double * upperInternal = modelPtr_->upperRegion();
    double rhsScale = modelPtr_->rhsScale();
    const double * columnScale = NULL;
    if (modelPtr_->scalingFlag()>0) 
      columnScale = modelPtr_->columnScale() ;
    // and do bounds in case dual needs them
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (columnLower[iColumn]>saveLowerOriginal[iColumn]) {
	double value = columnLower[iColumn];
	value *= rhsScale;
	if (columnScale)
	  value /= columnScale[iColumn];
	lowerInternal[iColumn]=value;
      }
      if (columnUpper[iColumn]<saveUpperOriginal[iColumn]) {
	double value = columnUpper[iColumn];
	value *= rhsScale;
	if (columnScale)
	  value /= columnScale[iColumn];
	upperInternal[iColumn]=value;
      }
    }
    // Start of fast iterations
    bool alwaysFinish= ((specialOptions_&32)==0) ? true : false;
    //modelPtr_->setLogLevel(1);
    modelPtr_->setIntParam(ClpMaxNumIteration,itlim);
    int status = ((ClpSimplexDual *)modelPtr_)->fastDual(alwaysFinish);
    
    int problemStatus = modelPtr_->problemStatus();
    double objectiveValue =modelPtr_->objectiveValue() * modelPtr_->optimizationDirection();
    CoinAssert (modelPtr_->problemStatus()||modelPtr_->objectiveValue()<1.0e50);
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
      modelPtr_->checkPrimalSolution(modelPtr_->solutionRegion(0),
                                     modelPtr_->solutionRegion(1));
      //modelPtr_->gutsOfSolution(NULL,NULL,0);
      //if (problemStatus==3)
      //modelPtr_->computeObjectiveValue();
      objectiveValue =modelPtr_->objectiveValue() *
	modelPtr_->optimizationDirection();
      obj = CoinMax(objectiveValue,saveObjectiveValue);
      if (!modelPtr_->numberDualInfeasibilities()) { 
	double limit = 0.0;
	modelPtr_->getDblParam(ClpDualObjectiveLimit, limit);
	if (modelPtr_->secondaryStatus()==1&&!problemStatus&&obj<limit) {
	  obj=limit;
	  problemStatus=3;
	}
	if (!modelPtr_->numberPrimalInfeasibilities()&&obj<limit) { 
	  problemStatus=0;
	} else if (problemStatus==10) {
	  problemStatus=3;
	} else if (!modelPtr_->numberPrimalInfeasibilities()) {
	  problemStatus=1; // infeasible
	} 
      } else {
	// can't say much
	//if (problemStatus==3)
	//modelPtr_->computeObjectiveValue();
	lastAlgorithm_=1; // so won't fail on cutoff (in CbcNode)
	problemStatus=3;
      }
    } else if (!problemStatus) {
      if (modelPtr_->isDualObjectiveLimitReached()) 
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
    double * solution = modelPtr_->primalColumnSolution();
    const double * solution2 = modelPtr_->solutionRegion();
    // could just do changed bounds - also try double size scale so can use * not /
    if (!columnScale) {
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	solution[iColumn]= solution2[iColumn];
      }
    } else {
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	solution[iColumn]= solution2[iColumn]*columnScale[iColumn];
      }
    }
    memcpy(columnLower,saveLowerOriginal,numberColumns*sizeof(double));
    memcpy(columnUpper,saveUpperOriginal,numberColumns*sizeof(double));
#if 0
    // could combine with loop above
    if (modelPtr_==modelPtr_)
      modelPtr_->computeObjectiveValue();
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
    // and back bounds
    memcpy(modelPtr_->lowerRegion(),saveLower,number*sizeof(double));
    memcpy(modelPtr_->upperRegion(),saveUpper,number*sizeof(double));
    modelPtr_->setIntParam(ClpMaxNumIteration,itlimOrig_);
    return;
  }
#endif
  if (smallModel_==NULL) {
    setWarmStart(ws_);
    memcpy(modelPtr_->primalRowSolution(),
           rowActivity_,numberRows*sizeof(double));
    memcpy(modelPtr_->primalColumnSolution(),columnActivity_,
           numberColumns*sizeof(double));
    modelPtr_->setIntParam(ClpMaxNumIteration,CoinMin(itlim,9999));
    resolve();
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
    smallModel_->setWhatsChanged(511);
    double * lowerSmall = smallModel_->lowerRegion();
    double * upperSmall = smallModel_->upperRegion();
    double * lowerSmallReal = smallModel_->columnLower();
    double * upperSmallReal = smallModel_->columnUpper();
    int i;
    double rhsScale = smallModel_->rhsScale();
    const double * columnScale = NULL;
    if (smallModel_->scalingFlag()>0) 
      columnScale = (smallModel_->auxiliaryModel_==NULL) ? smallModel_->columnScale() 
        : smallModel_->auxiliaryModel_->columnScale();
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
    bool alwaysFinish= ((specialOptions_&32)==0) ? true : false;
    //smallModel_->setLogLevel(1);
    smallModel_->setIntParam(ClpMaxNumIteration,itlim);
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
      //if (problemStatus==3)
      //smallModel_->computeObjectiveValue();
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
        //if (problemStatus==3)
        //smallModel_->computeObjectiveValue();
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
    if (!columnScale) {
      for (i=0;i<numberColumns2;i++) {
        int iColumn = whichColumn[i];
        solution[iColumn]= solution2[i];
        lowerSmallReal[i]=saveLowerOriginal[i];
        upperSmallReal[i]=saveUpperOriginal[i];
      }
    } else {
      for (i=0;i<numberColumns2;i++) {
        int iColumn = whichColumn[i];
        solution[iColumn]= solution2[i]*columnScale[i];
        lowerSmallReal[i]=saveLowerOriginal[i];
        upperSmallReal[i]=saveUpperOriginal[i];
      }
    }
    // could combine with loop above
    if (modelPtr_==smallModel_)
      modelPtr_->computeObjectiveValue();
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
  modelPtr_->setIntParam(ClpMaxNumIteration,itlimOrig_);
}

void OsiClpSolverInterface::unmarkHotStart()
{
#ifdef CLEAN_HOT_START
  if ((specialOptions_&65536)!=0) {
    modelPtr_->deleteRim(0);
    if (lastNumberRows_<0) {
      specialOptions_ |= 131072;
      modelPtr_->specialOptions_ |= 131072;
      lastNumberRows_ = -1 -lastNumberRows_;
      if (modelPtr_->rowScale_) {
	if (modelPtr_->rowScale_!=rowScale_.array()) {
	  delete [] modelPtr_->rowScale_;
	  delete [] modelPtr_->columnScale_;
	}
	modelPtr_->rowScale_=NULL;
	modelPtr_->columnScale_=NULL;
      }
    }
    delete factorization_;
    delete [] spareArrays_;
    smallModel_=NULL;
    spareArrays_=NULL;
    factorization_=NULL;
    delete [] rowActivity_;
    delete [] columnActivity_;
    rowActivity_=NULL;
    columnActivity_=NULL;
    return;
  }
#endif
  if (smallModel_==NULL) {
    setWarmStart(ws_);
    int numberRows = modelPtr_->numberRows();
    int numberColumns = modelPtr_->numberColumns();
    memcpy(modelPtr_->primalRowSolution(),
           rowActivity_,numberRows*sizeof(double));
    memcpy(modelPtr_->primalColumnSolution(),columnActivity_,
           numberColumns*sizeof(double));
    delete ws_;
    ws_ = NULL;
  } else {
    if (!modelPtr_->auxiliaryModel_) {
      if (smallModel_!=modelPtr_)
	delete smallModel_;
    } else {
      modelPtr_->deleteRim(0);
      //modelPtr_->setLogLevel(modelPtr_->auxiliaryModel_->numberPrimalInfeasibilities_);
      //modelPtr_->setIntParam(ClpMaxNumIteration,modelPtr_->auxiliaryModel_->numberDualInfeasibilities_);
    }
    delete factorization_;
    delete [] spareArrays_;
    smallModel_=NULL;
    spareArrays_=NULL;
    factorization_=NULL;
  }
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
  
// Get pointer to mutable column-wise copy of matrix (returns NULL if not meaningful)
CoinPackedMatrix * 
OsiClpSolverInterface::getMutableMatrixByCol() const 
{
  ClpPackedMatrix * matrix = dynamic_cast<ClpPackedMatrix *>(modelPtr_->matrix_) ;
  if (matrix)
    return matrix->getPackedMatrix();
  else
    return NULL;
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
/* Set the objective coefficients for all columns
    array [getNumCols()] is an array of values for the objective.
    This defaults to a series of set operations and is here for speed.
*/
void 
OsiClpSolverInterface::setObjective(const double * array)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  CoinMemcpyN(array,modelPtr_->numberColumns(),
		    modelPtr_->objective());
}
/* Set the lower bounds for all columns
    array [getNumCols()] is an array of values for the objective.
    This defaults to a series of set operations and is here for speed.
*/
void 
OsiClpSolverInterface::setColLower(const double * array)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  CoinMemcpyN(array,modelPtr_->numberColumns(),
		    modelPtr_->columnLower());
}
/* Set the upper bounds for all columns
    array [getNumCols()] is an array of values for the objective.
    This defaults to a series of set operations and is here for speed.
*/
void 
OsiClpSolverInterface::setColUpper(const double * array)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  CoinMemcpyN(array,modelPtr_->numberColumns(),
		    modelPtr_->columnUpper());
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setColSolution(const double * cs) 
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
  if (collb) {
    for (iCol = 0; iCol < numcols; iCol++) {
      lower[iCol]= forceIntoRange(collb[iCol], -OsiClpInfinity, OsiClpInfinity);
      if (lower[iCol]<-1.0e27)
	lower[iCol]=-COIN_DBL_MAX;
    }
  } else {
    CoinFillN ( lower, numcols,0.0);
  }
  if (colub) {
    for (iCol = 0; iCol < numcols; iCol++) {
      upper[iCol]= forceIntoRange(colub[iCol], -OsiClpInfinity, OsiClpInfinity);
      if (upper[iCol]>1.0e27)
	upper[iCol]=COIN_DBL_MAX;
    }
  } else {
    CoinFillN ( upper, numcols,COIN_DBL_MAX);
  }
  if (obj) {
    for (iCol = 0; iCol < numcols; iCol++) {
      objective[iCol] = obj[iCol];
    }
  } else {
    CoinFillN ( objective, numcols,0.0);
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
void 
OsiClpSolverInterface::addCols(const int numcols,
			       const int * columnStarts, const int * rows, const double * elements,
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
  if (collb) {
    for (iCol = 0; iCol < numcols; iCol++) {
      lower[iCol]= forceIntoRange(collb[iCol], -OsiClpInfinity, OsiClpInfinity);
      if (lower[iCol]<-1.0e27)
	lower[iCol]=-COIN_DBL_MAX;
    }
  } else {
    CoinFillN ( lower, numcols,0.0);
  }
  if (colub) {
    for (iCol = 0; iCol < numcols; iCol++) {
      upper[iCol]= forceIntoRange(colub[iCol], -OsiClpInfinity, OsiClpInfinity);
      if (upper[iCol]>1.0e27)
	upper[iCol]=COIN_DBL_MAX;
    }
  } else {
    CoinFillN ( upper, numcols,COIN_DBL_MAX);
  }
  if (obj) {
    for (iCol = 0; iCol < numcols; iCol++) {
      objective[iCol] = obj[iCol];
    }
  } else {
    CoinFillN ( objective, numcols,0.0);
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendCols(numcols,columnStarts,rows,elements);
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
  findIntegers(false);
  deleteBranchingInfo(num,columnIndices);
  modelPtr_->deleteColumns(num,columnIndices);
  int nameDiscipline;
  getIntParam(OsiNameDiscipline,nameDiscipline) ;
  if (num&&nameDiscipline) {
    // Very clumsy (and inefficient) - need to sort and then go backwards in ? chunks
    int * indices = CoinCopyOfArray(columnIndices,num);
    std::sort(indices,indices+num);
    int num2=num;
    while(num2) {
      int next = indices[num2-1];
      int firstDelete = num2-1;
      int i;
      for (i=num2-2;i>=0;i--) {
	if (indices[i]+1==next) {
	  next --;
	  firstDelete=i;
	} else {
	  break;
	}
      }
      OsiSolverInterface::deleteColNames(firstDelete,num2-firstDelete);
      num2 = firstDelete;
      assert (num2>=0);
    }
    delete [] indices;
  }
  // synchronize integers (again)
  if (integerInformation_) {
    int numberColumns = modelPtr_->numberColumns();
    for (int i=0;i<numberColumns;i++) {
      if (modelPtr_->isInteger(i))
	integerInformation_[i]=1;
      else
	integerInformation_[i]=0;
    }
  }
  basis_.deleteColumns(num,columnIndices);
  linearObjective_ = modelPtr_->objective();
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  setRowBounds(numberRows,rowlb,rowub);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults1();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  double rowlb = 0, rowub = 0;
  convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
  setRowBounds(numberRows,rowlb,rowub);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults1();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(int numberElements, const int * columns, const double * elements,
			   const double rowlb, const double rowub) 
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  setRowBounds(numberRows,rowlb,rowub);
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRow(numberElements, columns, elements);
  CoinBigIndex starts[2];
  starts[0]=0;
  starts[1]=numberElements;
  redoScaleFactors( 1,starts, columns, elements);
  freeCachedResults1();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    if (rowlb) 
      lower[iRow]= forceIntoRange(rowlb[iRow], -OsiClpInfinity, OsiClpInfinity);
    else 
      lower[iRow]=-OsiClpInfinity;
    if (rowub) 
      upper[iRow]= forceIntoRange(rowub[iRow], -OsiClpInfinity, OsiClpInfinity);
    else 
      upper[iRow]=OsiClpInfinity;
    if (lower[iRow]<-1.0e27)
      lower[iRow]=-COIN_DBL_MAX;
    if (upper[iRow]>1.0e27)
      upper[iRow]=COIN_DBL_MAX;
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRows(numrows,rows);
  freeCachedResults1();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    double rowlb = 0, rowub = 0;
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
  freeCachedResults1();
}
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const int * rowStarts, const int * columns, const double * element,
			       const double* rowlb, const double* rowub)
{
  freeCachedResults0();
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    if (rowlb) 
      lower[iRow]= forceIntoRange(rowlb[iRow], -OsiClpInfinity, OsiClpInfinity);
    else 
      lower[iRow]=-OsiClpInfinity;
    if (rowub) 
      upper[iRow]= forceIntoRange(rowub[iRow], -OsiClpInfinity, OsiClpInfinity);
    else 
      upper[iRow]=OsiClpInfinity;
    if (lower[iRow]<-1.0e27)
      lower[iRow]=-COIN_DBL_MAX;
    if (upper[iRow]>1.0e27)
      upper[iRow]=COIN_DBL_MAX;
  }
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  modelPtr_->matrix()->appendRows(numrows,rowStarts,columns,element);
  redoScaleFactors( numrows,rowStarts, columns, element);
  freeCachedResults1();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  // will still be optimal if all rows basic
  bool allBasic=true;
  int numBasis = basis_.getNumArtificial();
  for (int i=0;i<num;i++) {
    int iRow = rowIndices[i];
    if (iRow<numBasis) {
      if (basis_.getArtifStatus(iRow)!=CoinWarmStartBasis::basic) {
	allBasic=false;
	break;
      }
    }
  }
  int saveAlgorithm = allBasic ? lastAlgorithm_ : 999;
  modelPtr_->deleteRows(num,rowIndices);
  int nameDiscipline;
  getIntParam(OsiNameDiscipline,nameDiscipline) ;
  if (num&&nameDiscipline) {
    // Very clumsy (and inefficient) - need to sort and then go backwards in ? chunks
    int * indices = CoinCopyOfArray(rowIndices,num);
    std::sort(indices,indices+num);
    int num2=num;
    while(num2) {
      int next = indices[num2-1];
      int firstDelete = num2-1;
      int i;
      for (i=num2-2;i>=0;i--) {
	if (indices[i]+1==next) {
	  next --;
	  firstDelete=i;
	} else {
	  break;
	}
      }
      OsiSolverInterface::deleteRowNames(firstDelete,num2-firstDelete);
      num2 = firstDelete;
      assert (num2>=0);
    }
    delete [] indices;
  }
  basis_.deleteRows(num,rowIndices);
  freeCachedResults();
  lastAlgorithm_ = saveAlgorithm;
  if ((specialOptions_&131072)!=0) 
    lastNumberRows_=modelPtr_->numberRows();
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
  // assert( rowsen != NULL );
  // assert( rowrhs != NULL );
  // If any of Rhs NULLs then create arrays
  int numrows = matrix.getNumRows();
  const char * rowsenUse = rowsen;
  if (!rowsen) {
    char * rowsen = new char [numrows];
    for (int i=0;i<numrows;i++)
      rowsen[i]='G';
    rowsenUse = rowsen;
  } 
  const double * rowrhsUse = rowrhs;
  if (!rowrhs) {
    double * rowrhs = new double [numrows];
    for (int i=0;i<numrows;i++)
      rowrhs[i]=0.0;
    rowrhsUse = rowrhs;
  }
  const double * rowrngUse = rowrng;
  if (!rowrng) {
    double * rowrng = new double [numrows];
    for (int i=0;i<numrows;i++)
      rowrng[i]=0.0;
    rowrngUse = rowrng;
  }
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsenUse[i],rowrhsUse[i],rowrngUse[i],rowlb[i],rowub[i]);
  }
  if (rowsen!=rowsenUse)
    delete [] rowsenUse;
  if (rowrhs!=rowrhsUse)
    delete [] rowrhsUse;
  if (rowrng!=rowrngUse)
    delete [] rowrngUse;
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
  // If any of Rhs NULLs then create arrays
  const char * rowsenUse = rowsen;
  if (!rowsen) {
    char * rowsen = new char [numrows];
    for (int i=0;i<numrows;i++)
      rowsen[i]='G';
    rowsenUse = rowsen;
  } 
  const double * rowrhsUse = rowrhs;
  if (!rowrhs) {
    double * rowrhs = new double [numrows];
    for (int i=0;i<numrows;i++)
      rowrhs[i]=0.0;
    rowrhsUse = rowrhs;
  }
  const double * rowrngUse = rowrng;
  if (!rowrng) {
    double * rowrng = new double [numrows];
    for (int i=0;i<numrows;i++)
      rowrng[i]=0.0;
    rowrngUse = rowrng;
  }
  double * rowlb = new double[numrows];
  double * rowub = new double[numrows];
  for (int i = numrows-1; i >= 0; --i) {   
    convertSenseToBound(rowsenUse[i],rowrhsUse[i],rowrngUse[i],rowlb[i],rowub[i]);
  }
  if (rowsen!=rowsenUse)
    delete [] rowsenUse;
  if (rowrhs!=rowrhsUse)
    delete [] rowrhsUse;
  if (rowrng!=rowrngUse)
    delete [] rowrngUse;
  loadProblem(numcols, numrows, start,  index, value, collb, colub, obj,
	      rowlb,  rowub);
  delete[] rowlb;
  delete[] rowub;
}
// This loads a model from a coinModel object - returns number of errors
int 
OsiClpSolverInterface::loadFromCoinModel (  CoinModel & modelObject, bool keepSolution)
{
  int numberErrors = 0;
  // Set arrays for normal use
  double * rowLower = modelObject.rowLowerArray();
  double * rowUpper = modelObject.rowUpperArray();
  double * columnLower = modelObject.columnLowerArray();
  double * columnUpper = modelObject.columnUpperArray();
  double * objective = modelObject.objectiveArray();
  int * integerType = modelObject.integerTypeArray();
  double * associated = modelObject.associatedArray();
  // If strings then do copies
  if (modelObject.stringsExist()) {
    numberErrors = modelObject.createArrays(rowLower, rowUpper, columnLower, columnUpper,
                                            objective, integerType,associated);
  }
  CoinPackedMatrix matrix;
  modelObject.createPackedMatrix(matrix,associated);
  int numberRows = modelObject.numberRows();
  int numberColumns = modelObject.numberColumns();
  CoinWarmStart * ws = getWarmStart();
  bool restoreBasis = keepSolution && numberRows&&numberRows==getNumRows()&&
    numberColumns==getNumCols();
  loadProblem(matrix, 
              columnLower, columnUpper, objective, rowLower, rowUpper);
  if (restoreBasis)
    setWarmStart(ws);
  delete ws;
  // Do names if wanted
  int numberItems;
  numberItems = modelObject.rowNames()->numberItems();
  if (numberItems) {
    const char *const * rowNames=modelObject.rowNames()->names();
    modelPtr_->copyRowNames(rowNames,0,numberItems);
  }
  numberItems = modelObject.columnNames()->numberItems();
  if (numberItems) {
    const char *const * columnNames=modelObject.columnNames()->names();
    modelPtr_->copyColumnNames(columnNames,0,numberItems);
  }
  // Do integers if wanted
  assert(integerType);
  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
    if (integerType[iColumn])
      setInteger(iColumn);
  }
  if (rowLower!=modelObject.rowLowerArray()||
      columnLower!=modelObject.columnLowerArray()) {
    delete [] rowLower;
    delete [] rowUpper;
    delete [] columnLower;
    delete [] columnUpper;
    delete [] objective;
    delete [] integerType;
    delete [] associated;
    //if (numberErrors)
    //  handler_->message(CLP_BAD_STRING_VALUES,messages_)
    //    <<numberErrors
    //    <<CoinMessageEol;
  }
  modelPtr_->optimizationDirection_ = modelObject.optimizationDirection();  
  return numberErrors;
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
  // get names
  const char * const * const rowNames = modelPtr_->rowNamesAsChar();
  const char * const * const columnNames = modelPtr_->columnNamesAsChar();
  // Fall back on Osi version - possibly with names
  OsiSolverInterface::writeMpsNative(fullname.c_str(), 
				     const_cast<const char **>(rowNames),
                                     const_cast<const char **>(columnNames),0,2,objSense,
				     numberSOS_,setInfo_);
  if (rowNames) {
    modelPtr_->deleteNamesAsChar(rowNames, modelPtr_->numberRows_+1);
    modelPtr_->deleteNamesAsChar(columnNames, modelPtr_->numberColumns_);
  }
}

int 
OsiClpSolverInterface::writeMpsNative(const char *filename, 
		  const char ** rowNames, const char ** columnNames,
		  int formatType,int numberAcross,double objSense) const 
{
  return OsiSolverInterface::writeMpsNative(filename, rowNames, columnNames,
			       formatType, numberAcross,objSense,
					    numberSOS_,setInfo_);
}

//#############################################################################
// CLP specific public interfaces
//#############################################################################

ClpSimplex * OsiClpSolverInterface::getModelPtr() const
{
  int saveAlgorithm = lastAlgorithm_;
  freeCachedResults();
  lastAlgorithm_ = saveAlgorithm;
  //bool inCbcOrOther = (modelPtr_->specialOptions()&0x03000000)!=0;
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
numberSOS_(0),
setInfo_(NULL),
smallModel_(NULL),
factorization_(NULL),
smallestElementInCut_(1.0e-15),
smallestChangeInCut_(1.0e-10),
spareArrays_(NULL),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL),
cleanupScaling_(0),
specialOptions_(0x80000000),
baseModel_(NULL),
lastNumberRows_(0)
{
  //printf("in default %x\n",this);
  modelPtr_=NULL;
  notOwned_=false;
  reset();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * OsiClpSolverInterface::clone(bool CopyData) const
{
  //printf("in clone %x\n",this);
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
: OsiSolverInterface(rhs),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
rowActivity_(NULL),
columnActivity_(NULL),
numberSOS_(rhs.numberSOS_),
setInfo_(NULL),
smallModel_(NULL),
factorization_(NULL),
smallestElementInCut_(rhs.smallestElementInCut_),
smallestChangeInCut_(rhs.smallestChangeInCut_),
spareArrays_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL)
{
  //printf("in copy %x - > %x\n",&rhs,this);
  if ( rhs.modelPtr_  ) 
    modelPtr_ = new ClpSimplex(*rhs.modelPtr_);
  else
    modelPtr_ = new ClpSimplex();
  if ( rhs.baseModel_  ) 
    baseModel_ = new ClpSimplex(*rhs.baseModel_);
  else
    baseModel_ = NULL;
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
  solveOptions_  = rhs.solveOptions_;
  cleanupScaling_ = rhs.cleanupScaling_;
  specialOptions_ = rhs.specialOptions_;
  lastNumberRows_ = rhs.lastNumberRows_;
  rowScale_ = rhs.rowScale_;
  columnScale_ = rhs.columnScale_;
  fillParamMaps();
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  if (numberSOS_) {
    setInfo_ = new CoinSet[numberSOS_];
    for (int i=0;i<numberSOS_;i++)
      setInfo_[i]=rhs.setInfo_[i];
  }
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
numberSOS_(0),
setInfo_(NULL),
smallModel_(NULL),
factorization_(NULL),
smallestElementInCut_(1.0e-15),
smallestChangeInCut_(1.0e-10),
spareArrays_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
notOwned_(false),
matrixByRow_(NULL),
integerInformation_(NULL),
whichRange_(NULL),
cleanupScaling_(0),
specialOptions_(0x80000000),
baseModel_(NULL),
lastNumberRows_(0)
{
  //printf("in borrow %x - > %x\n",&rhs,this);
  modelPtr_ = rhs;
  basis_.resize(modelPtr_->numberRows(),modelPtr_->numberColumns());
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
  //printf("in destructor %x\n",this);
  freeCachedResults();
  if (!notOwned_)
    delete modelPtr_;
  delete baseModel_;
  delete ws_;
  delete [] rowActivity_;
  delete [] columnActivity_;
  delete [] setInfo_;
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
    //printf("in = %x - > %x\n",&rhs,this);
    OsiSolverInterface::operator=(rhs);
    freeCachedResults();
    if (!notOwned_)
      delete modelPtr_;
    delete ws_;
    if ( rhs.modelPtr_  ) 
      modelPtr_ = new ClpSimplex(*rhs.modelPtr_);
    delete baseModel_;
    if ( rhs.baseModel_  ) 
      baseModel_ = new ClpSimplex(*rhs.baseModel_);
    else
      baseModel_ = NULL;
    notOwned_=false;
    linearObjective_ = modelPtr_->objective();
    saveData_ = rhs.saveData_;
    solveOptions_  = rhs.solveOptions_;
    cleanupScaling_ = rhs.cleanupScaling_;
    specialOptions_ = rhs.specialOptions_;
    lastNumberRows_ = rhs.lastNumberRows_;
    rowScale_ = rhs.rowScale_;
    columnScale_ = rhs.columnScale_;
    basis_ = rhs.basis_;
    if (rhs.integerInformation_) {
      int numberColumns = modelPtr_->numberColumns();
      integerInformation_ = new char[numberColumns];
      memcpy(integerInformation_,rhs.integerInformation_,
	     numberColumns*sizeof(char));
    }
    if ( rhs.ws_ ) 
      ws_ = new CoinWarmStartBasis(*rhs.ws_);
    else
      ws_=NULL;
    delete [] rowActivity_;
    delete [] columnActivity_;
    rowActivity_=NULL;
    columnActivity_=NULL;
    delete [] setInfo_;
    numberSOS_ = rhs.numberSOS_;
    setInfo_=NULL;
    if (numberSOS_) {
      setInfo_ = new CoinSet[numberSOS_];
      for (int i=0;i<numberSOS_;i++)
	setInfo_[i]=rhs.setInfo_[i];
    }
    assert(smallModel_==NULL);
    assert(factorization_==NULL);
    smallestElementInCut_ = rhs.smallestElementInCut_;
    smallestChangeInCut_ = rhs.smallestChangeInCut_;
    assert(spareArrays_==NULL);
    basis_ = rhs.basis_;
    fillParamMaps();
    messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiClpSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  applyRowCuts(1, &rowCut);
}
/* Apply a collection of row cuts which are all effective.
   applyCuts seems to do one at a time which seems inefficient.
*/
void 
OsiClpSolverInterface::applyRowCuts(int numberCuts, const OsiRowCut * cuts)
{
  if (numberCuts) {
    // Say can't gurantee optimal basis etc
    lastAlgorithm_=999;

    // Thanks to js
    const OsiRowCut * * cutsp = new const OsiRowCut * [numberCuts];
    for (int i=0;i<numberCuts;i++) 
      cutsp[i] = &cuts[i];
    
    applyRowCuts(numberCuts, cutsp);
    
    delete [] cutsp;
  }
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
#if 0 // was #ifndef NDEBUG
  int nameDiscipline;
  getIntParam(OsiNameDiscipline,nameDiscipline) ;
  assert (!nameDiscipline);
#endif
  freeCachedResults0();
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numberCuts,modelPtr_->numberColumns());
  basis_.resize(numberRows+numberCuts,modelPtr_->numberColumns());
  // redo as relaxed - use modelPtr_-> addRows with starts etc
  int size = 0;
  for (i=0;i<numberCuts;i++) 
    size += cuts[i]->row().getNumElements();
  CoinBigIndex * starts = new CoinBigIndex [numberCuts+1];
  int * indices = new int[size];
  double * elements = new double[size];
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  const double * columnLower = modelPtr_->columnLower();
  const double * columnUpper = modelPtr_->columnUpper();
  size=0;
  for (i=0;i<numberCuts;i++) {
    double rowLb = cuts[i]->lb();
    double rowUb = cuts[i]->ub();
    int n=cuts[i]->row().getNumElements();
    const int * index = cuts[i]->row().getIndices();
    const double * elem = cuts[i]->row().getElements();
    starts[i]=size;
    for (int j=0;j<n;j++) {
      double value = elem[j];
      int column = index[j];
      if (fabs(value)>=smallestChangeInCut_) {
        // always take
        indices[size]=column;
        elements[size++]=value;
      } else if (fabs(value)>=smallestElementInCut_) {
        double lowerValue = columnLower[column];
        double upperValue = columnUpper[column];
        double difference = upperValue-lowerValue;
        if (difference<1.0e20&&difference*fabs(value)<smallestChangeInCut_&&
            (rowLb<-1.0e20||rowUb>1.0e20)) {
          // Take out and adjust to relax
          //printf("small el %g adjusted\n",value);
          if (rowLb>-1.0e20) {
            // just lower bound on row
            if (value>0.0) {
              // pretend at upper
              rowLb -= value*upperValue;
            } else {
              // pretend at lower
              rowLb -= value*lowerValue;
            }
          } else {
            // just upper bound on row
            if (value>0.0) {
              // pretend at lower
              rowUb -= value*lowerValue;
            } else {
              // pretend at upper
              rowUb -= value*upperValue;
            }
          }
        } else {
          // take (unwillingly)
          indices[size]=column;
          elements[size++]=value;
        }
      } else {
        //printf("small el %g ignored\n",value);
      }
    }
    lower[i]= forceIntoRange(rowLb, -OsiClpInfinity, OsiClpInfinity);
    upper[i]= forceIntoRange(rowUb, -OsiClpInfinity, OsiClpInfinity);
    if (lower[i]<-1.0e27)
      lower[i]=-COIN_DBL_MAX;
    if (upper[i]>1.0e27)
      upper[i]=COIN_DBL_MAX;
  }
  starts[numberCuts]=size;
  if (!modelPtr_->clpMatrix())
    modelPtr_->createEmptyMatrix();
  //modelPtr_->matrix()->appendRows(numberCuts,rows);
  modelPtr_->clpMatrix()->appendMatrix(numberCuts,0,starts,indices,elements);
  freeCachedResults1();
  redoScaleFactors( numberCuts,starts, indices, elements);
  delete [] starts;
  delete [] indices;
  delete [] elements;

}
// Extend scale factors
void 
OsiClpSolverInterface::redoScaleFactors(int numberAdd,const CoinBigIndex * starts,
					const int * indices, const double * elements)
{
  if ((specialOptions_&131072)!=0) {
    int numberRows = modelPtr_->numberRows()-numberAdd;
    assert (lastNumberRows_==numberRows); // ???
    int iRow;
    int newNumberRows = numberRows + numberAdd;
    rowScale_.extend(2*newNumberRows*sizeof(double));
    double * rowScale = rowScale_.array();
    double * oldInverseScale = rowScale + lastNumberRows_;
    double * inverseRowScale = rowScale + newNumberRows;
    for (iRow=lastNumberRows_-1;iRow>=0;iRow--)
      inverseRowScale[iRow] = oldInverseScale[iRow] ;
    //int numberColumns = baseModel_->numberColumns();
    const double * columnScale = columnScale_.array();
    //const double * inverseColumnScale = columnScale + numberColumns;
    // Geometric mean on row scales
    // adjust arrays
    rowScale += lastNumberRows_;
    inverseRowScale += lastNumberRows_;
    for (iRow=0;iRow<numberAdd;iRow++) {
      CoinBigIndex j;
      double largest=1.0e-20;
      double smallest=1.0e50;
      for (j=starts[iRow];j<starts[iRow+1];j++) {
	int iColumn = indices[j];
	double value = fabs(elements[j]);
	// Don't bother with tiny elements
	if (value>1.0e-20) {
	  value *= columnScale[iColumn];
	  largest = CoinMax(largest,value);
	  smallest = CoinMin(smallest,value);
	}
      }
      double scale=sqrt(smallest*largest);
      scale=CoinMax(1.0e-10,CoinMin(1.0e10,scale));
      inverseRowScale[iRow]=scale;
      rowScale[iRow]=1.0/scale;
    }
    lastNumberRows_=newNumberRows;
  }
}
// Delete all scale factor stuff and reset option
void OsiClpSolverInterface::deleteScaleFactors()
{
  delete baseModel_;
  baseModel_=NULL;
  lastNumberRows_=0;
  specialOptions_ &= ~131072;
}
//-----------------------------------------------------------------------------

void OsiClpSolverInterface::applyColCut( const OsiColCut & cc )
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  delete [] rowsense_;
  delete [] rhs_;
  delete [] rowrange_;
  delete matrixByRow_;
  //delete ws_;
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
  matrixByRow_=NULL;
  //ws_ = NULL;
  if (modelPtr_&&modelPtr_->clpMatrix()) {
    modelPtr_->clpMatrix()->refresh(modelPtr_); // make sure all clean
#ifndef NDEBUG
    ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (modelPtr_->clpMatrix());
    if (clpMatrix) {
      assert (clpMatrix->getNumRows()==modelPtr_->getNumRows());
      assert (clpMatrix->getNumCols()==modelPtr_->getNumCols());
    }
#endif
  }
}

//------------------------------------------------------------------- 

void OsiClpSolverInterface::freeCachedResults0() const
{  
  delete [] rowsense_;
  delete [] rhs_;
  delete [] rowrange_;
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
}

//------------------------------------------------------------------- 

void OsiClpSolverInterface::freeCachedResults1() const
{  
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  delete matrixByRow_;
  matrixByRow_=NULL;
  //ws_ = NULL;
  if (modelPtr_&&modelPtr_->clpMatrix()) {
    modelPtr_->clpMatrix()->refresh(modelPtr_); // make sure all clean
#ifndef NDEBUG
    ClpPackedMatrix * clpMatrix = dynamic_cast<ClpPackedMatrix *> (modelPtr_->clpMatrix());
    if (clpMatrix) {
      assert (clpMatrix->getNumRows()==modelPtr_->getNumRows());
      assert (clpMatrix->getNumCols()==modelPtr_->getNumCols());
    }
#endif
  }
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
   assert ((int) OsiMaxNumIteration==        (int)ClpMaxNumIteration);
   assert ((int) OsiMaxNumIterationHotStart==(int)ClpMaxNumIterationHotStart);
   //assert ((int) OsiLastIntParam==           (int)ClpLastIntParam);

   assert ((int) OsiDualObjectiveLimit==  (int)ClpDualObjectiveLimit);
   assert ((int) OsiPrimalObjectiveLimit==(int)ClpPrimalObjectiveLimit);
   assert ((int) OsiDualTolerance==       (int)ClpDualTolerance);
   assert ((int) OsiPrimalTolerance==     (int)ClpPrimalTolerance);
   assert ((int) OsiObjOffset==           (int)ClpObjOffset);
   //assert ((int) OsiLastDblParam==        (int)ClpLastDblParam);

   assert ((int) OsiProbName==    (int) ClpProbName);
   //strParamMap_[OsiLastStrParam] = ClpLastStrParam;
}
// Sets up basis
void 
OsiClpSolverInterface::setBasis ( const CoinWarmStartBasis & basis)
{
  setBasis(basis,modelPtr_);
  setWarmStart(&basis); 
}
//#define NEW_STATUS
#ifdef NEW_STATUS
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
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
    int stat = basis2.getArtifStatus(iRow);
    model->setRowStatus(iRow, (ClpSimplex::Status) stat);
  }
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    model->setColumnStatus(iColumn,
		    (ClpSimplex::Status) basis2.getStructStatus(iColumn));
  }
}
#else
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
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
#endif
/* Read an mps file from the given filename (defaults to Osi reader) - returns
   number of errors (see OsiMpsReader class) */
int 
OsiClpSolverInterface::readMps(const char *filename,
			       const char *extension ) 
{
  // Get rid of integer stuff
  delete [] integerInformation_;
  integerInformation_=NULL;
  freeCachedResults();
  
  CoinMpsIO m;
  m.setInfinity(getInfinity());
  m.passInMessageHandler(modelPtr_->messageHandler());
  *m.messagesPointer()=modelPtr_->coinMessages();

  delete [] setInfo_;
  setInfo_=NULL;
  numberSOS_=0;
  CoinSet ** sets=NULL;
  int numberErrors = m.readMps(filename,extension,numberSOS_,sets);
  if (numberSOS_) {
    setInfo_ = new CoinSet[numberSOS_];
    for (int i=0;i<numberSOS_;i++) {
      setInfo_[i]=*sets[i];
      delete sets[i];
    }
    delete [] sets;
  }
  handler_->message(COIN_SOLVER_MPS,messages_)
    <<m.getProblemName()<< numberErrors <<CoinMessageEol;
  if (!numberErrors) {

    // set objective function offest
    setDblParam(OsiObjOffset,m.objectiveOffset());

    // set problem name
    setStrParam(OsiProbName,m.getProblemName());
    
    // no errors
    loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		m.getObjCoefficients(),m.getRowSense(),m.getRightHandSide(),
		m.getRowRange());
    const char * integer = m.integerColumns();
    int nCols=m.getNumCols();
    int nRows=m.getNumRows();
    if (integer) {
      int i,n=0;
      int * index = new int [nCols];
      for (i=0;i<nCols;i++) {
	if (integer[i]) {
	  index[n++]=i;
        }
      }
      setInteger(index,n);
      delete [] index;
      if (n) 
        modelPtr_->copyInIntegerInformation(integer);
    }

    // set objective name
    setObjName(m.getObjectiveName());

    // Always keep names
    int nameDiscipline;
    getIntParam(OsiNameDiscipline,nameDiscipline) ;
    int iRow;
    std::vector<std::string> rowNames = std::vector<std::string> ();
    std::vector<std::string> columnNames = std::vector<std::string> ();
    rowNames.reserve(nRows);
    for (iRow=0;iRow<nRows;iRow++) {
      const char * name = m.rowName(iRow);
      rowNames.push_back(name);
      if (nameDiscipline) 
	OsiSolverInterface::setRowName(iRow,name) ;
    }
    
    int iColumn;
    columnNames.reserve(nCols);
    for (iColumn=0;iColumn<nCols;iColumn++) {
      const char * name = m.columnName(iColumn);
      columnNames.push_back(name);
      if (nameDiscipline) 
	OsiSolverInterface::setColName(iColumn,name) ;
    }
    modelPtr_->copyNames(rowNames,columnNames);
  }
  return numberErrors;
}
/* Read an mps file from the given filename returns
   number of errors (see OsiMpsReader class) */
int 
OsiClpSolverInterface::readMps(const char *filename,bool keepNames,bool allowErrors)
{
  // Get rid of integer stuff
  delete [] integerInformation_;
  integerInformation_=NULL;
  freeCachedResults();
  
  CoinMpsIO m;
  m.setInfinity(getInfinity());
  m.passInMessageHandler(modelPtr_->messageHandler());
  *m.messagesPointer()=modelPtr_->coinMessages();

  delete [] setInfo_;
  setInfo_=NULL;
  numberSOS_=0;
  CoinSet ** sets=NULL;
  int numberErrors = m.readMps(filename,"",numberSOS_,sets);
  if (numberSOS_) {
    setInfo_ = new CoinSet[numberSOS_];
    for (int i=0;i<numberSOS_;i++) {
      setInfo_[i]=*sets[i];
      delete sets[i];
    }
    delete [] sets;
  }
  handler_->message(COIN_SOLVER_MPS,messages_)
    <<m.getProblemName()<< numberErrors <<CoinMessageEol;
  if (!numberErrors||((numberErrors>0&&numberErrors<100000)&&allowErrors)) {

    // set objective function offest
    setDblParam(OsiObjOffset,m.objectiveOffset());

    // set problem name
    setStrParam(OsiProbName,m.getProblemName());

    // set objective name
    setObjName(m.getObjectiveName());

    // no errors
    loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		m.getObjCoefficients(),m.getRowSense(),m.getRightHandSide(),
		m.getRowRange());
    int nCols=m.getNumCols();
    // get quadratic part
    if (m.reader()->whichSection (  ) == COIN_QUAD_SECTION ) {
      int * start=NULL;
      int * column = NULL;
      double * element = NULL;
      int status=m.readQuadraticMps(NULL,start,column,element,2);
      if (!status) 
	modelPtr_->loadQuadraticObjective(nCols,start,column,element);
      delete [] start;
      delete [] column;
      delete [] element;
    }
    const char * integer = m.integerColumns();
    int nRows=m.getNumRows();
    if (integer) {
      int i,n=0;
      int * index = new int [nCols];
      for (i=0;i<nCols;i++) {
	if (integer[i]) {
	  index[n++]=i;
        }
      }
      setInteger(index,n);
      delete [] index;
      if (n) 
        modelPtr_->copyInIntegerInformation(integer);
    }
    if (keepNames) {
      // keep names
      int nameDiscipline;
      getIntParam(OsiNameDiscipline,nameDiscipline) ;
      int iRow;
      std::vector<std::string> rowNames = std::vector<std::string> ();
      std::vector<std::string> columnNames = std::vector<std::string> ();
      rowNames.reserve(nRows);
      for (iRow=0;iRow<nRows;iRow++) {
	const char * name = m.rowName(iRow);
	rowNames.push_back(name);
	if (nameDiscipline) 
	  OsiSolverInterface::setRowName(iRow,name) ;
      }
      
      int iColumn;
      columnNames.reserve(nCols);
      for (iColumn=0;iColumn<nCols;iColumn++) {
	const char * name = m.columnName(iColumn);
	columnNames.push_back(name);
	if (nameDiscipline) 
	  OsiSolverInterface::setColName(iColumn,name) ;
      }
      modelPtr_->copyNames(rowNames,columnNames);
    }
  }
  return numberErrors;
}
// Read file in LP format (with names)
int 
OsiClpSolverInterface::readLp(const char *filename, const double epsilon )
{
  CoinLpIO m;
  m.readLp(filename, epsilon);
  freeCachedResults();

  // set objective function offest
  setDblParam(OsiObjOffset, 0);

  // set problem name
  setStrParam(OsiProbName, m.getProblemName());

  // set objective name
  setObjName(m.getObjName());

  // no errors
  loadProblem(*m.getMatrixByRow(), m.getColLower(), m.getColUpper(),
	      m.getObjCoefficients(), m.getRowLower(), m.getRowUpper());

  const char *integer = m.integerColumns();
  int nCols = m.getNumCols();
  int nRows = m.getNumRows();
  if (integer) {
    int i, n = 0;
    int *index = new int [nCols];
    for (i=0; i<nCols; i++) {
      if (integer[i]) {
	index[n++] = i;
      }
    }
    setInteger(index,n);
    delete [] index;
  }
  // Always keep names
  int nameDiscipline;
  getIntParam(OsiNameDiscipline,nameDiscipline) ;
  int iRow;
  std::vector<std::string> rowNames = std::vector<std::string> ();
  std::vector<std::string> columnNames = std::vector<std::string> ();
  rowNames.reserve(nRows);
  for (iRow=0;iRow<nRows;iRow++) {
    const char * name = m.rowName(iRow);
    rowNames.push_back(name);
    if (nameDiscipline) 
      OsiSolverInterface::setRowName(iRow,name) ;
  }
  
  int iColumn;
  columnNames.reserve(nCols);
  for (iColumn=0;iColumn<nCols;iColumn++) {
    const char * name = m.columnName(iColumn);
    columnNames.push_back(name);
    if (nameDiscipline) 
      OsiSolverInterface::setColName(iColumn,name) ;
  }
  modelPtr_->copyNames(rowNames,columnNames);
  return(0);
}
/* Write the problem into an Lp file of the given filename.
   If objSense is non zero then -1.0 forces the code to write a
   maximization objective and +1.0 to write a minimization one.
   If 0.0 then solver can do what it wants.
   This version calls writeLpNative with names */
void 
OsiClpSolverInterface::writeLp(const char *filename,
                               const char *extension ,
                               double epsilon ,
                               int numberAcross ,
                               int decimals ,
                               double objSense ,
                               bool changeNameOnRange) const
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
  // get names
  const char * const * const rowNames = modelPtr_->rowNamesAsChar();
  const char * const * const columnNames = modelPtr_->columnNamesAsChar();
  // Fall back on Osi version - possibly with names
  OsiSolverInterface::writeLpNative(fullname.c_str(), 
				    rowNames,columnNames, epsilon, numberAcross,
				    decimals, objSense,changeNameOnRange);
  if (rowNames) {
    modelPtr_->deleteNamesAsChar(rowNames, modelPtr_->numberRows_+1);
    modelPtr_->deleteNamesAsChar(columnNames, modelPtr_->numberColumns_);
  }
}
void 
OsiClpSolverInterface::writeLp(FILE * fp,
                               double epsilon ,
                               int numberAcross ,
                               int decimals ,
                               double objSense ,
                               bool changeNameOnRange) const
{
  // get names
  const char * const * const rowNames = modelPtr_->rowNamesAsChar();
  const char * const * const columnNames = modelPtr_->columnNamesAsChar();
  // Fall back on Osi version - possibly with names
  OsiSolverInterface::writeLpNative(fp,
				    rowNames,columnNames, epsilon, numberAcross,
				    decimals, objSense,changeNameOnRange);
  if (rowNames) {
    modelPtr_->deleteNamesAsChar(rowNames, modelPtr_->numberRows_+1);
    modelPtr_->deleteNamesAsChar(columnNames, modelPtr_->numberColumns_);
  }
}
/*
  I (JJF) am getting incredibly annoyed because I can't just replace a matrix.
  The default behavior of this is do nothing so only use where that would not matter
  e.g. strengthening a matrix for MIP
*/
void 
OsiClpSolverInterface::replaceMatrixOptional(const CoinPackedMatrix & matrix)
{
  replaceMatrix(matrix);
}
// And if it does matter (not used at present)
void 
OsiClpSolverInterface::replaceMatrix(const CoinPackedMatrix & matrix)
{
  delete modelPtr_->matrix_;
  delete modelPtr_->rowCopy_;
  modelPtr_->rowCopy_=NULL;
  if (matrix.isColOrdered()) {
    modelPtr_->matrix_=new ClpPackedMatrix(matrix);
  } else {
    CoinPackedMatrix matrix2;
    matrix2.reverseOrderedCopyOf(matrix);
    modelPtr_->matrix_=new ClpPackedMatrix(matrix2);
  }    
  modelPtr_->matrix_->setDimensions(modelPtr_->numberRows_,modelPtr_->numberColumns_);
  freeCachedResults();
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
    //printf("obj a %g %g\n",modelPtr_->objectiveValue(),
    //     OsiSolverInterface::getObjValue());
    return modelPtr_->objectiveValue();
  } else {
    return OsiSolverInterface::getObjValue();
  }
}

/* Set an objective function coefficient */
void 
OsiClpSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
OsiClpSolverInterface::setColLower( int index, double elementValue )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (index<0||index>=n) {
    indexError(index,"setColLower");
  }
#endif
  double currentValue = modelPtr_->columnActivity_[index];
  bool changed=(currentValue<elementValue-modelPtr_->primalTolerance()||
                index>=basis_.getNumStructural()||
                basis_.getStructStatus(index)==CoinWarmStartBasis::atLowerBound);
  // Say can't gurantee optimal basis etc
  if (changed)
    lastAlgorithm_=999;
  if (!modelPtr_->lower_)
    modelPtr_->whatsChanged_=0; // switch off
  modelPtr_->setColumnLower(index,elementValue);
}
      
/* Set a single column upper bound<br>
   Use DBL_MAX for infinity. */
void 
OsiClpSolverInterface::setColUpper( int index, double elementValue )
{
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (index<0||index>=n) {
    indexError(index,"setColUpper");
  }
#endif
  double currentValue = modelPtr_->columnActivity_[index];
  bool changed=(currentValue>elementValue+modelPtr_->primalTolerance()||
                index>=basis_.getNumStructural()||
                basis_.getStructStatus(index)==CoinWarmStartBasis::atUpperBound);
  // Say can't gurantee optimal basis etc
  if (changed)
    lastAlgorithm_=999;
  if (!modelPtr_->upper_)
    modelPtr_->whatsChanged_=0; // switch off
  modelPtr_->setColumnUpper(index,elementValue);
}

/* Set a single column lower and upper bound */
void 
OsiClpSolverInterface::setColBounds( int elementIndex,
				     double lower, double upper )
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
#ifndef NDEBUG
  int n = modelPtr_->numberColumns();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setColBounds");
  }
#endif
  if (!modelPtr_->lower_)
    modelPtr_->whatsChanged_=0; // switch off
  modelPtr_->setColumnBounds(elementIndex,lower,upper);
}
void OsiClpSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowLower");
  }
#endif
  modelPtr_->setRowLower(elementIndex , elementValue);
  if (rowsense_!=NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(modelPtr_->rowLower_[elementIndex], 
			modelPtr_->rowUpper_[elementIndex],
			rowsense_[elementIndex], rhs_[elementIndex], rowrange_[elementIndex]);
  }
}
      
/* Set a single row upper bound<br>
   Use DBL_MAX for infinity. */
void 
OsiClpSolverInterface::setRowUpper( int elementIndex, double elementValue ) {
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowUpper");
  }
#endif
  modelPtr_->setRowUpper(elementIndex , elementValue);
  if (rowsense_!=NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(modelPtr_->rowLower_[elementIndex], 
			modelPtr_->rowUpper_[elementIndex],
			rowsense_[elementIndex], rhs_[elementIndex], rowrange_[elementIndex]);
  }
}
    
/* Set a single row lower and upper bound */
void 
OsiClpSolverInterface::setRowBounds( int elementIndex,
	      double lower, double upper ) {
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (elementIndex<0||elementIndex>=n) {
    indexError(elementIndex,"setRowBounds");
  }
#endif
  modelPtr_->setRowBounds(elementIndex,lower,upper);
  if (rowsense_!=NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(modelPtr_->rowLower_[elementIndex], 
			modelPtr_->rowUpper_[elementIndex],
			rowsense_[elementIndex], rhs_[elementIndex], rowrange_[elementIndex]);
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (i<0||i>=n) {
    indexError(i,"setRowType");
  }
#endif
  double lower = 0, upper = 0;
  convertSenseToBound(sense, rightHandSide, range, lower, upper);
  setRowBounds(i, lower, upper);
  // If user is using sense then set
  if (rowsense_) {
    rowsense_[i] = sense;
    rhs_[i] = rightHandSide;
    rowrange_[i] = range;
  }
}
// Set name of row
void 
//OsiClpSolverInterface::setRowName(int rowIndex, std::string & name) 
OsiClpSolverInterface::setRowName(int rowIndex, std::string name) 
{
  if (rowIndex>=0&&rowIndex<modelPtr_->numberRows()) {
    int nameDiscipline;
    getIntParam(OsiNameDiscipline,nameDiscipline) ;
    if (nameDiscipline) {
      modelPtr_->setRowName(rowIndex,name);
      OsiSolverInterface::setRowName(rowIndex,name) ;
    }
  }
}
// Return name of row if one exists or Rnnnnnnn
// we ignore maxLen
std::string 
OsiClpSolverInterface::getRowName(int rowIndex, unsigned maxLen) const
{ 
	if (rowIndex == getNumRows())
		return getObjName();
  return modelPtr_->getRowName(rowIndex);
}
    
// Set name of col
void 
//OsiClpSolverInterface::setColName(int colIndex, std::string & name) 
OsiClpSolverInterface::setColName(int colIndex, std::string name) 
{
  if (colIndex>=0&&colIndex<modelPtr_->numberColumns()) {
    int nameDiscipline;
    getIntParam(OsiNameDiscipline,nameDiscipline) ;
    if (nameDiscipline) {
      modelPtr_->setColumnName(colIndex,name);
      OsiSolverInterface::setColName(colIndex,name) ;
    }
  }
}
// Return name of col if one exists or Rnnnnnnn
std::string 
OsiClpSolverInterface::getColName(int colIndex, unsigned maxLen) const
{
  return modelPtr_->getColumnName(colIndex);
}
    
    
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
    double lowerValue = 0;
    double upperValue = 0;
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
  int saveIts = modelPtr_->numberIterations_;
  modelPtr_->setSolveType(2);
  if (doingPrimal)
    modelPtr_->setAlgorithm(1);
  else
    modelPtr_->setAlgorithm(-1);
  // Do initialization
  saveData_ = modelPtr_->saveData();
  saveData_.scalingFlag_=modelPtr_->scalingFlag();
  modelPtr_->scaling(0);
  specialOptions_ = 0x80000000;
  // set infeasibility cost up
  modelPtr_->setInfeasibilityCost(1.0e12);
  ClpDualRowDantzig dantzig;
  modelPtr_->setDualRowPivotAlgorithm(dantzig);
  ClpPrimalColumnDantzig dantzigP;
  dantzigP.saveWeights(modelPtr_,0); // set modelPtr 
  modelPtr_->setPrimalColumnPivotAlgorithm(dantzigP);
#ifdef NDEBUG
  modelPtr_->startup(0);
#else
  int returnCode=modelPtr_->startup(0);
  assert (!returnCode);
#endif
  modelPtr_->numberIterations_=saveIts;
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
void 
OsiClpSolverInterface::enableFactorization() const
{
  saveData_.scalingFlag_=specialOptions_;
  int saveStatus = modelPtr_->problemStatus_;
  if ((specialOptions_&(1+8))!=1+8)
    setSpecialOptionsMutable((1+8)|specialOptions_);
#ifdef NDEBUG
  modelPtr_->startup(0);
#else
  int returnCode=modelPtr_->startup(0);
  assert (!returnCode);
#endif
  modelPtr_->problemStatus_=saveStatus;
}

//Undo whatever setting changes the above method had to make
void 
OsiClpSolverInterface::disableFactorization() const
{
  specialOptions_=saveData_.scalingFlag_;
  // declare optimality anyway  (for message handler)
  modelPtr_->setProblemStatus(0);
  // message will not appear anyway
  int saveMessageLevel=modelPtr_->messageHandler()->logLevel();
  modelPtr_->messageHandler()->setLogLevel(0);
  // Should re-do - for moment save arrays 
  double * sol = CoinCopyOfArray(modelPtr_->columnActivity_,modelPtr_->numberColumns_);
  double * dj = CoinCopyOfArray(modelPtr_->reducedCost_,modelPtr_->numberColumns_);
  double * rsol = CoinCopyOfArray(modelPtr_->rowActivity_,modelPtr_->numberRows_);
  double * dual = CoinCopyOfArray(modelPtr_->dual_,modelPtr_->numberRows_);
  modelPtr_->finish();
  memcpy(modelPtr_->columnActivity_,sol,modelPtr_->numberColumns_*sizeof(double));
  memcpy(modelPtr_->reducedCost_,dj,modelPtr_->numberColumns_*sizeof(double));
  memcpy(modelPtr_->rowActivity_,rsol,modelPtr_->numberRows_*sizeof(double));
  memcpy(modelPtr_->dual_,dual,modelPtr_->numberRows_*sizeof(double));
  delete [] sol;
  delete [] dj;
  delete [] rsol;
  delete [] dual;
  modelPtr_->messageHandler()->setLogLevel(saveMessageLevel);
}
#ifdef NEW_STATUS
/* The following two methods may be replaced by the
   methods of OsiSolverInterface using OsiWarmStartBasis if:
   1. OsiWarmStartBasis resize operation is implemented
   more efficiently and
   2. It is ensured that effects on the solver are the same
   
   Returns a basis status of the structural/artificial variables 
*/
void 
OsiClpSolverInterface::getBasisStatus(int* cstat, int* rstat) const
{
  int iRow,iColumn;
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const double * pi = modelPtr_->dualRowSolution();
  const double * dj = modelPtr_->dualColumnSolution();
  double multiplier = modelPtr_->optimizationDirection();

  int lookup[]={0,1,2,3,0,3};
  for (iRow=0;iRow<numberRows;iRow++) {
    int iStatus = modelPtr_->getRowStatus(iRow);
    if (iStatus==5) {
      // Fixed - look at reduced cost
      if (pi[iRow]*multiplier<-1.0e-7)
        iStatus = 2;
    }
    iStatus = lookup[iStatus];
    rstat[iRow]=iStatus;
  }
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int iStatus = modelPtr_->getColumnStatus(iColumn);
    if (iStatus==5) {
      // Fixed - look at reduced cost
      if (dj[iColumn]*multiplier<-1.0e-7)
        iStatus = 2;
    }
    iStatus = lookup[iStatus];
    cstat[iColumn]=iStatus;
  }
}

//Set the status of structural/artificial variables 
int 
OsiClpSolverInterface::setBasisStatus(const int* cstat, const int* rstat)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
  modelPtr_->createStatus();
  int i, n;
  double * lower, * upper, * solution;
  n=modelPtr_->numberRows();
  lower = modelPtr_->rowLower();
  upper = modelPtr_->rowUpper();
  solution = modelPtr_->primalRowSolution();
  // For rows lower and upper are flipped
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
  // Save 
  basis_ = getBasis(modelPtr_);
  return 0;
}
#else
/* The following two methods may be replaced by the
   methods of OsiSolverInterface using OsiWarmStartBasis if:
   1. OsiWarmStartBasis resize operation is implemented
   more efficiently and
   2. It is ensured that effects on the solver are the same
   
   Returns a basis status of the structural/artificial variables 
*/
void 
OsiClpSolverInterface::getBasisStatus(int* cstat, int* rstat) const
{
  int iRow,iColumn;
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const double * pi = modelPtr_->dualRowSolution();
  const double * dj = modelPtr_->dualColumnSolution();
  double multiplier = modelPtr_->optimizationDirection();
  // Flip slacks
  int lookupA[]={0,1,3,2,0,3};
  for (iRow=0;iRow<numberRows;iRow++) {
    int iStatus = modelPtr_->getRowStatus(iRow);
    if (iStatus==5) {
      // Fixed - look at reduced cost
      if (pi[iRow]*multiplier>1.0e-7)
        iStatus = 3;
    }
    iStatus = lookupA[iStatus];
    rstat[iRow]=iStatus;
  }
  int lookupS[]={0,1,2,3,0,3};
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int iStatus = modelPtr_->getColumnStatus(iColumn);
    if (iStatus==5) {
      // Fixed - look at reduced cost
      if (dj[iColumn]*multiplier<-1.0e-7)
        iStatus = 2;
    }
    iStatus = lookupS[iStatus];
    cstat[iColumn]=iStatus;
  }
}

//Set the status of structural/artificial variables 
//Returns 0 if OK, 1 if problem is bad e.g. duplicate elements, too large ...
int 
OsiClpSolverInterface::setBasisStatus(const int* cstat, const int* rstat)
{
  // Say can't gurantee optimal basis etc
  lastAlgorithm_=999;
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
  // May be bad model
  if (modelPtr_->status()==4)
    return 1;
  // Save 
  basis_ = getBasis(modelPtr_);
  return 0;
}
#endif

/* Perform a pivot by substituting a colIn for colOut in the basis. 
   The status of the leaving variable is given in statOut. Where
   1 is to upper bound, -1 to lower bound
   Return code is 0 for okay,
   1 if inaccuracy forced re-factorization (should be okay) and
   -1 for singular factorization
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
  return modelPtr_->pivot();
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
   Primal and dual solutions are updated
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
OsiClpSolverInterface::getBInvARow(int row, double* z, double * slack) const
{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvARow");
  }
#endif
  //assert (modelPtr_->solveType()==2||(specialOptions_&1));
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
  // If user is sophisticated then let her/him do work
  if ((specialOptions_&512)==0) {
    // otherwise copy and clear
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
    columnArray0->clear();
    rowArray1->clear();
  }
  // don't need to clear everything always, but doesn't cost
  rowArray0->clear();
  columnArray1->clear();
}
//Get a row of the tableau (slack part in slack if not NULL)
void 
OsiClpSolverInterface::getBInvARow(int row, CoinIndexedVector * columnArray0, CoinIndexedVector * slack,
				   bool keepScaled) const
{
#ifndef NDEBUG
  int nx = modelPtr_->numberRows();
  if (row<0||row>=nx) {
    indexError(row,"getBInvARow");
  }
#endif
  //assert (modelPtr_->solveType()==2||(specialOptions_&1));
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = slack ? slack : modelPtr_->rowArray(1);
  CoinIndexedVector * columnArray1 = modelPtr_->columnArray(1);
  rowArray0->clear();
  rowArray1->clear();
  columnArray0->clear();
  columnArray1->clear();
  //int numberRows = modelPtr_->numberRows();
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
  int n;
  const int * which;
  double * array;
  // deal with scaling etc
  if (rowScale&&!keepScaled) {
    int j;
    // First columns
    n = columnArray0->getNumElements();
    which = columnArray0->getIndices();
    array = columnArray0->denseVector();
    for (j=0; j < n; j++) {
      int k=which[j];
      array[k] /= columnScale[k];
    }
    if (slack) {
      n = slack->getNumElements();
      which = slack->getIndices();
      array = slack->denseVector();
      for(j=0; j < n; j++) {
	int k=which[j];
	array[k] *= rowScale[k];
      }
    }
  }
  if (!slack)
    rowArray1->clear();
}

//Get a row of the basis inverse
void 
OsiClpSolverInterface::getBInvRow(int row, double* z) const

{
#ifndef NDEBUG
  int n = modelPtr_->numberRows();
  if (row<0||row>=n) {
    indexError(row,"getBInvRow");
  }
#endif
  //assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  ClpFactorization * factorization = modelPtr_->factorization();
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  CoinIndexedVector * rowArray1 = modelPtr_->rowArray(1);
  rowArray0->clear();
  rowArray1->clear();
  // put +1 in row
  // But swap if pivot variable was slack as clp stores slack as -1.0
  double value = (modelPtr_->pivotVariable()[row]<modelPtr_->numberColumns()) ? 1.0 : -1.0;
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const double * rowScale = modelPtr_->rowScale();
  const double * columnScale = modelPtr_->columnScale();
  const int * pivotVariable = modelPtr_->pivotVariable();
  // but scale
  if (rowScale) {
    int pivot = pivotVariable[row];
    if (pivot<numberColumns) 
      value *= columnScale[pivot];
    else
      value /= rowScale[pivot-numberColumns];
  }
  rowArray1->insert(row,value);
  factorization->updateColumnTranspose(rowArray0,rowArray1);
  // If user is sophisticated then let her/him do work
  if ((specialOptions_&512)==0) {
    // otherwise copy and clear
    if (!rowScale) {
      memcpy(z,rowArray1->denseVector(),modelPtr_->numberRows()*sizeof(double));
    } else {
      double * array = rowArray1->denseVector();
      for (int i=0;i<numberRows;i++) {
        z[i] = array[i] * rowScale[i];
      }
    }
    rowArray1->clear();
  }
}

//Get a column of the tableau
void 
OsiClpSolverInterface::getBInvACol(int col, double* vec) const
{
  //assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
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
  // If user is sophisticated then let her/him do work
  if ((specialOptions_&512)==0) {
    // otherwise copy and clear
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
}
//Get a column of the tableau
void 
OsiClpSolverInterface::getBInvACol(int col, CoinIndexedVector * rowArray1) const
{
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  rowArray0->clear();
  rowArray1->clear();
  // get column of matrix
#ifndef NDEBUG
  int nx = modelPtr_->numberColumns()+modelPtr_->numberRows();
  if (col<0||col>=nx) {
    indexError(col,"getBInvACol");
  }
#endif
  //int numberRows = modelPtr_->numberRows();
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
  // Deal with stuff
  int n = rowArray1->getNumElements();
  const int * which = rowArray1->getIndices();
  double * array = rowArray1->denseVector();
  for(int j=0; j < n; j++){
    int k=which[j];
    // need to know pivot variable for +1/-1 (slack) and row/column scaling
    int pivot = pivotVariable[k];
    if (pivot<numberColumns) {
      if (columnScale) 
	array[k] *= columnScale[pivot];
    } else {
      if (!rowScale) {
	array[k] = -array[k];
      } else {
	array[k] = -array[k]/rowScale[pivot-numberColumns];
      }
    }
  }
}

//Get an updated column
void 
OsiClpSolverInterface::getBInvACol(CoinIndexedVector * rowArray1) const
{
  CoinIndexedVector * rowArray0 = modelPtr_->rowArray(0);
  rowArray0->clear();
  // get column of matrix
  //int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const int * pivotVariable = modelPtr_->pivotVariable();
  const double * rowScale = modelPtr_->rowScale();
  const double * columnScale = modelPtr_->columnScale();
  // rowArray1 is not a column - so column scale can't be applied before
  modelPtr_->factorization()->updateColumn(rowArray0,rowArray1,false);
  // Deal with stuff
  int n = rowArray1->getNumElements();
  const int * which = rowArray1->getIndices();
  double * array = rowArray1->denseVector();
  for(int j=0; j < n; j++){
    int k=which[j];
    // need to know pivot variable for +1/-1 (slack) and row/column scaling
    int pivot = pivotVariable[k];
    if (pivot<numberColumns) {
      if (columnScale) 
	array[k] *= columnScale[pivot];
    } else {
      if (!rowScale) {
	array[k] = -array[k];
      } else {
	array[k] = -array[k]/rowScale[pivot-numberColumns];
      }
    }
  }
}

//Get a column of the basis inverse
void 
OsiClpSolverInterface::getBInvCol(int col, double* vec) const
{
  //assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
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
  int numberRows = modelPtr_->numberRows();
  int numberColumns = modelPtr_->numberColumns();
  const double * rowScale = modelPtr_->rowScale();
  const double * columnScale = modelPtr_->columnScale();
  const int * pivotVariable = modelPtr_->pivotVariable();
  // but scale
  double value;
  if (!rowScale) {
    value=1.0;
  } else {
    value = rowScale[col];
  }
  rowArray1->insert(col,value);
  factorization->updateColumn(rowArray0,rowArray1,false);
  // If user is sophisticated then let her/him do work
  if ((specialOptions_&512)==0) {
    // otherwise copy and clear
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
        double value = array[i];
        if (pivot<numberColumns) 
          vec[i] = value * columnScale[pivot];
        else
          vec[i] = - value / rowScale[pivot-numberColumns];
      }
    }
    rowArray1->clear();
  }
}

/* Get basic indices (order of indices corresponds to the
   order of elements in a vector returned by getBInvACol() and
   getBInvCol()).
*/
void 
OsiClpSolverInterface::getBasics(int* index) const
{
  //assert (modelPtr_->solveType()==2||(specialOptions_&1)!=0);
  assert (index);
  assert (modelPtr_->pivotVariable());
  memcpy(index,modelPtr_->pivotVariable(),
	 modelPtr_->numberRows()*sizeof(int));
}
//Returns true if a basis is available and optimal
bool 
OsiClpSolverInterface::basisIsAvailable() const 
{
  return (lastAlgorithm_==1||lastAlgorithm_==2)&&(!modelPtr_->problemStatus_);
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
  ws_ = NULL;
  delete [] rowActivity_;
  delete [] columnActivity_;
  assert(smallModel_==NULL);
  assert(factorization_==NULL);
  smallestElementInCut_ = 1.0e-15;
  smallestChangeInCut_ = 1.0e-10;
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
      if ( specialOptions_==0x80000000) {
        setupForRepeatedUse(0,0);
        specialOptions_=0;
      }
      // set normal
      specialOptions_ &= (1023+3*8192+3*65536);
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
  bool tightenBounds = ((specialOptions_&64)==0) ? false : true; 
  ClpSimplex * small = ((ClpSimplexOther *) modelPtr_)->crunch(rhs,whichRow,whichColumn,
                                                               nBound,false,tightenBounds);
  if (small) {
    if ((specialOptions_&131072)!=0) {
      assert (lastNumberRows_>=0);
      int numberRows2 = small->numberRows();
      int numberColumns2 = small->numberColumns();
      double * rowScale2 = new double [2*numberRows2];
      assert (rowScale_.getSize()>=2*numberRows);
      const double * rowScale = rowScale_.array();
      double * inverseScale2 = rowScale2+numberRows2;
      const double * inverseScale = rowScale+modelPtr_->numberRows_;
      int i;
      for (i=0;i<numberRows2;i++) {
	int iRow = whichRow[i];
	//assert (iRow<numberRows);
	rowScale2[i]=rowScale[iRow];
	inverseScale2[i]=inverseScale[iRow];
      }
      small->setRowScale(rowScale2);
      double * columnScale2 = new double [2*numberColumns2];
      assert (columnScale_.getSize()>=2*numberColumns);
      const double * columnScale = columnScale_.array();
      inverseScale2 = columnScale2+numberColumns2;
      inverseScale = columnScale+modelPtr_->numberColumns_;
      for (i=0;i<numberColumns2;i++) {
	int iColumn = whichColumn[i];
	//assert (iColumn<numberColumns);
	columnScale2[i]=columnScale[iColumn];
	inverseScale2[i]=inverseScale[iColumn];
      }
      small->setColumnScale(columnScale2);
      small->specialOptions_ |= 131072;
    }
    OsiClpDisasterHandler handler(this);
    bool inCbcOrOther = (modelPtr_->specialOptions()&0x03000000)!=0;
    if (inCbcOrOther) {
      handler.setSimplex(small);
      handler.setWhereFrom(1); // crunch
      small->setDisasterHandler(&handler);
    }
    small->dual();
    if (small->problemStatus()==0) {
      modelPtr_->setProblemStatus(0);
      ((ClpSimplexOther *) modelPtr_)->afterCrunch(*small,whichRow,whichColumn,nBound);
    } else if (small->problemStatus()!=3) {
      modelPtr_->setProblemStatus(1);
    } else {
      if (small->problemStatus_==3) {
	// may be problems
	if (inCbcOrOther&&handler.inTrouble()) {
	  // in case scaling bad
	  small->setRowScale(NULL);
	  small->setColumnScale(NULL);
    	  // try just going back in
	  handler.setPhase(1);
	  small->dual();
	  if (handler.inTrouble()) {
	    // try primal on original model
	    handler.setPhase(2);
	    handler.setOsiModel(this);
	    modelPtr_->setDisasterHandler(&handler);
	    modelPtr_->primal();
	    if(handler.inTrouble()) {
#ifdef COIN_DEVELOP
	      printf("disaster crunch - treat as infeasible\n");
#endif
	      modelPtr_->setProblemStatus(1);
	    }
	    // give up for now
	    modelPtr_->setDisasterHandler(NULL);
	  }
	} else {
	  small->computeObjectiveValue();
	  modelPtr_->setObjectiveValue(small->objectiveValue());
	  modelPtr_->setProblemStatus(3);
	}
      } else {
	modelPtr_->setProblemStatus(3);
      }
    }
    delete small;
  } else {
    modelPtr_->setProblemStatus(1);
  }
  delete [] whichRow;
  delete [] whichColumn;
}
// Synchronize model 
void 
OsiClpSolverInterface::synchronizeModel() 
{
  if ((specialOptions_ &128)!=0) {
    if (!modelPtr_->rowScale_&&(specialOptions_&131072)!=0) {
      assert (lastNumberRows_==modelPtr_->numberRows_);
      int numberRows = modelPtr_->numberRows();
      int numberColumns = modelPtr_->numberColumns();
      double * rowScale = CoinCopyOfArray(rowScale_.array(),2*numberRows);
      modelPtr_->setRowScale(rowScale);
      double * columnScale = CoinCopyOfArray(columnScale_.array(),2*numberColumns);
      modelPtr_->setColumnScale(columnScale);
      modelPtr_->auxiliaryModel(63-2);
      modelPtr_->setRowScale(NULL);
      modelPtr_->setColumnScale(NULL);
    } else {
      modelPtr_->auxiliaryModel(63-2);
    }
  }
}
// Returns true if has OsiSimplex methods
/* Returns 1 if can just do getBInv etc
   2 if has all OsiSimplex methods
   and 0 if it has none */
int
OsiClpSolverInterface::canDoSimplexInterface() const
{
  return 2;
}
// Pass in sos stuff from AMPl
void 
OsiClpSolverInterface::setSOSData(int numberSOS,const char * type,
				  const int * start,const int * indices, const double * weights)
{
  delete [] setInfo_;
  setInfo_=NULL;
  numberSOS_=numberSOS;
  if (numberSOS_) {
    setInfo_ = new CoinSet[numberSOS_];
    for (int i=0;i<numberSOS_;i++) {
      int iStart = start[i];
      setInfo_[i]=CoinSosSet(start[i+1]-iStart,indices+iStart,weights ? weights+iStart : NULL,
			     type[i]);
    }
  }
}
/* Identify integer variables and SOS and create corresponding objects.
  
      Record integer variables and create an OsiSimpleInteger object for each
      one.  All existing OsiSimpleInteger objects will be destroyed.
      If the solver supports SOS then do the same for SOS.

      If justCount then no objects created and we just store numberIntegers_
      Returns number of SOS
*/
int 
OsiClpSolverInterface::findIntegersAndSOS(bool justCount)
{
  findIntegers(justCount);
  int nObjects=0;
  OsiObject ** oldObject = object_;
  int iObject;
  int numberSOS=0;
  for (iObject = 0;iObject<numberObjects_;iObject++) {
    OsiSOS * obj =
      dynamic_cast <OsiSOS *>(oldObject[iObject]) ;
    if (obj) 
      numberSOS++;
  }
  if (numberSOS_&&!numberSOS) {
    // make a large enough array for new objects
    nObjects = numberObjects_;
    numberObjects_=numberSOS_+nObjects;
    if (numberObjects_)
      object_ = new OsiObject * [numberObjects_];
    else
      object_=NULL;
    // copy
    memcpy(object_,oldObject,nObjects*sizeof(OsiObject *));
    // Delete old array (just array)
    delete [] oldObject;
    
    for (int i=0;i<numberSOS_;i++) {
      CoinSet * set =  setInfo_+i;
      object_[nObjects++] =
	new OsiSOS(this,set->numberEntries(),set->which(),set->weights(),
		   set->setType());
    }
  } else if (!numberSOS_&&numberSOS) {
    // create Coin sets
    assert (!setInfo_);
    setInfo_ = new CoinSet[numberSOS];
    for (iObject = 0;iObject<numberObjects_;iObject++) {
      OsiSOS * obj =
	dynamic_cast <OsiSOS *>(oldObject[iObject]) ;
      if (obj) 
	setInfo_[numberSOS_++]=CoinSosSet(obj->numberMembers(),obj->members(),obj->weights(),obj->sosType());
    }
  } else if (numberSOS!=numberSOS_) {
    printf("mismatch on SOS\n");
  }
  return numberSOS_;
}
// below needed for pathetic branch and bound code
#include <vector>
#include <map>

// Trivial class for Branch and Bound

class OsiNodeSimple  {
  
public:
    
  // Default Constructor 
  OsiNodeSimple ();

  // Constructor from current state (and list of integers)
  // Also chooses branching variable (if none set to -1)
  OsiNodeSimple (OsiSolverInterface &model,
                 int numberIntegers, int * integer,
                 CoinWarmStart * basis);
  
  // Copy constructor 
  OsiNodeSimple ( const OsiNodeSimple &);
   
  // Assignment operator 
  OsiNodeSimple & operator=( const OsiNodeSimple& rhs);

  // Destructor 
  ~OsiNodeSimple ();
  
  // Public data
  // Basis (should use tree, but not as wasteful as bounds!)
  CoinWarmStart * basis_;
  // Objective value
  double objectiveValue_;
  // Branching variable (0 is first integer)
  int variable_;
  // Way to branch - -1 down (first), 1 up, -2 down (second), 2 up (second)
  int way_;
  // Number of integers (for length of arrays)
  int numberIntegers_;
  // Current value
  double value_;
  // Now I must use tree
  // Bounds stored in full (for integers)
  int * lower_;
  int * upper_;
};


OsiNodeSimple::OsiNodeSimple() :
  basis_(NULL),
  objectiveValue_(1.0e100),
  variable_(-100),
  way_(-1),
  numberIntegers_(0),
  value_(0.5),
  lower_(NULL),
  upper_(NULL)
{
}
OsiNodeSimple::OsiNodeSimple(OsiSolverInterface & model,
		 int numberIntegers, int * integer,CoinWarmStart * basis)
{
  basis_ = basis;
  variable_=-1;
  way_=-1;
  numberIntegers_=numberIntegers;
  value_=0.0;
  if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
    objectiveValue_ = model.getObjSense()*model.getObjValue();
  } else {
    objectiveValue_ = 1.0e100;
    lower_ = NULL;
    upper_ = NULL;
    return; // node cutoff
  }
  lower_ = new int [numberIntegers_];
  upper_ = new int [numberIntegers_];
  assert (upper_!=NULL);
  const double * lower = model.getColLower();
  const double * upper = model.getColUpper();
  const double * solution = model.getColSolution();
  int i;
  // Hard coded integer tolerance
#define INTEGER_TOLERANCE 1.0e-6
  // Number of strong branching candidates
#define STRONG_BRANCHING 5
#ifdef STRONG_BRANCHING
  double upMovement[STRONG_BRANCHING];
  double downMovement[STRONG_BRANCHING];
  double solutionValue[STRONG_BRANCHING];
  int chosen[STRONG_BRANCHING];
  int iSmallest=0;
  // initialize distance from integer
  for (i=0;i<STRONG_BRANCHING;i++) {
    upMovement[i]=0.0;
    chosen[i]=-1;
  }
#endif
  variable_=-1;
  // This has hard coded integer tolerance
  double mostAway=INTEGER_TOLERANCE;
  int numberAway=0;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    lower_[i]=(int)lower[iColumn];
    upper_[i]=(int)upper[iColumn];
    double value = solution[iColumn];
    value = max(value,(double) lower_[i]);
    value = min(value,(double) upper_[i]);
    double nearest = floor(value+0.5);
    if (fabs(value-nearest)>INTEGER_TOLERANCE)
      numberAway++;
    if (fabs(value-nearest)>mostAway) {
#ifdef STRONG_BRANCHING
      double away = fabs(value-nearest);
      if (away>upMovement[iSmallest]) {
	//add to list
	upMovement[iSmallest]=away;
	solutionValue[iSmallest]=value;
	chosen[iSmallest]=i;
	int j;
	iSmallest=-1;
	double smallest = 1.0;
	for (j=0;j<STRONG_BRANCHING;j++) {
	  if (upMovement[j]<smallest) {
	    smallest=upMovement[j];
	    iSmallest=j;
	  }
	}
      }
#else
      mostAway=fabs(value-nearest);
      variable_=i;
      value_=value;
      if (value<=nearest)
	way_=1; // up
      else
	way_=-1; // down
#endif
    }
  }
#ifdef STRONG_BRANCHING
  int numberStrong=0;
  for (i=0;i<STRONG_BRANCHING;i++) {
    if (chosen[i]>=0) { 
      numberStrong ++;
      variable_ = chosen[i];
    }
  }
  // out strong branching if bit set
  OsiClpSolverInterface* clp =
    dynamic_cast<OsiClpSolverInterface*>(&model);
  if (clp&&(clp->specialOptions()&16)!=0&&numberStrong>1) {
    int j;
    int iBest=-1;
    double best = 0.0;
    for (j=0;j<STRONG_BRANCHING;j++) {
      if (upMovement[j]>best) {
        best=upMovement[j];
        iBest=j;
      }
    }
    numberStrong=1;
    variable_=chosen[iBest];
  }
  if (numberStrong==1) {
    // just one - makes it easy
    int iColumn = integer[variable_];
    double value = solution[iColumn];
    value = max(value,(double) lower_[variable_]);
    value = min(value,(double) upper_[variable_]);
    double nearest = floor(value+0.5);
    value_=value;
    if (value<=nearest)
      way_=1; // up
    else
      way_=-1; // down
  } else if (numberStrong) {
    // more than one - choose
    bool chooseOne=true;
    model.markHotStart();
    for (i=0;i<STRONG_BRANCHING;i++) {
      int iInt = chosen[i];
      if (iInt>=0) {
	int iColumn = integer[iInt];
	double value = solutionValue[i]; // value of variable in original
	double objectiveChange;
	value = max(value,(double) lower_[iInt]);
	value = min(value,(double) upper_[iInt]);

	// try down

	model.setColUpper(iColumn,floor(value));
	model.solveFromHotStart();
	model.setColUpper(iColumn,upper_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	downMovement[i]=objectiveChange;

	// try up

	model.setColLower(iColumn,ceil(value));
	model.solveFromHotStart();
	model.setColLower(iColumn,lower_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	upMovement[i]=objectiveChange;
	
	/* Possibilities are:
	   Both sides feasible - store
	   Neither side feasible - set objective high and exit
	   One side feasible - change bounds and resolve
	*/
	bool solveAgain=false;
	if (upMovement[i]<1.0e100) {
	  if(downMovement[i]<1.0e100) {
	    // feasible - no action
	  } else {
	    // up feasible, down infeasible
	    solveAgain = true;
	    model.setColLower(iColumn,ceil(value));
	  }
	} else {
	  if(downMovement[i]<1.0e100) {
	    // down feasible, up infeasible
	    solveAgain = true;
	    model.setColUpper(iColumn,floor(value));
	  } else {
	    // neither side feasible
	    objectiveValue_=1.0e100;
	    chooseOne=false;
	    break;
	  }
	}
	if (solveAgain) {
	  // need to solve problem again - signal this
	  variable_ = numberIntegers;
	  chooseOne=false;
	  break;
	}
      }
    }
    if (chooseOne) {
      // choose the one that makes most difference both ways
      double best = -1.0;
      double best2 = -1.0;
      for (i=0;i<STRONG_BRANCHING;i++) {
	int iInt = chosen[i];
	if (iInt>=0) {
	  //std::cout<<"Strong branching on "
          //   <<i<<""<<iInt<<" down "<<downMovement[i]
          //   <<" up "<<upMovement[i]
          //   <<" value "<<solutionValue[i]
          //   <<std::endl;
	  bool better = false;
	  if (min(upMovement[i],downMovement[i])>best) {
	    // smaller is better
	    better=true;
	  } else if (min(upMovement[i],downMovement[i])>best-1.0e-5) {
	    if (max(upMovement[i],downMovement[i])>best2+1.0e-5) {
	      // smaller is about same, but larger is better
	      better=true;
	    }
	  }
	  if (better) {
	    best = min(upMovement[i],downMovement[i]);
	    best2 = max(upMovement[i],downMovement[i]);
	    variable_ = iInt;
	    double value = solutionValue[i];
	    value = max(value,(double) lower_[variable_]);
	    value = min(value,(double) upper_[variable_]);
	    value_=value;
	    if (upMovement[i]<=downMovement[i])
	      way_=1; // up
	    else
	      way_=-1; // down
	  }
	}
      }
    }
    // Delete the snapshot
    model.unmarkHotStart();
  }
#endif
}

OsiNodeSimple::OsiNodeSimple(const OsiNodeSimple & rhs) 
{  
  basis_=rhs.basis_->clone();
  objectiveValue_=rhs.objectiveValue_;
  variable_=rhs.variable_;
  way_=rhs.way_;
  numberIntegers_=rhs.numberIntegers_;
  value_=rhs.value_;
  lower_=NULL;
  upper_=NULL;
  if (rhs.lower_!=NULL) {
    lower_ = new int [numberIntegers_];
    upper_ = new int [numberIntegers_];
    assert (upper_!=NULL);
    memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
    memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
  }
}

OsiNodeSimple &
OsiNodeSimple::operator=(const OsiNodeSimple & rhs)
{
  if (this != &rhs) {
    delete basis_;
    basis_=rhs.basis_->clone();
    objectiveValue_=rhs.objectiveValue_;
    variable_=rhs.variable_;
    way_=rhs.way_;
    numberIntegers_=rhs.numberIntegers_;
    value_=rhs.value_;
    delete [] lower_;
    delete [] upper_;
    lower_=NULL;
    upper_=NULL;
    if (rhs.lower_!=NULL) {
      lower_ = new int [numberIntegers_];
      upper_ = new int [numberIntegers_];
      assert (upper_!=NULL);
      memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
      memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
    }
  }
  return *this;
}


OsiNodeSimple::~OsiNodeSimple ()
{
  delete [] lower_;
  delete [] upper_;
  delete basis_;
}

#include <vector>

// Vector of OsiNodeSimples 
typedef std::vector<OsiNodeSimple>    OsiVectorNode;

// Invoke solver's built-in enumeration algorithm
void 
OsiClpSolverInterface::branchAndBound() {
  // solve LP
  initialSolve();

  if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
    int numberIntegers=0;
    int numberColumns = getNumCols();
    int iColumn;
    int i;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn))
        numberIntegers++;
    }
    if (!numberIntegers) {
      std::cout<<"No integer variables"
               <<std::endl;
      return;
    }
    int * which = new int[numberIntegers]; // which variables are integer
    numberIntegers=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn))
        which[numberIntegers++]=iColumn;
    }
    double direction = getObjSense();
    // empty tree
    OsiVectorNode branchingTree;
    
    // Add continuous to it;
    OsiNodeSimple rootNode(*this,numberIntegers,which,getWarmStart());
    // something extra may have been fixed by strong branching
    // if so go round again
    while (rootNode.variable_==numberIntegers) {
      resolve();
      rootNode = OsiNodeSimple(*this,numberIntegers,which,getWarmStart());
    }
    if (rootNode.objectiveValue_<1.0e100) {
      // push on stack
      branchingTree.push_back(rootNode);
    }
    
    // For printing totals
    int numberIterations=0;
    int numberNodes =0;
    
    OsiNodeSimple bestNode;
    // while until nothing on stack
    while (branchingTree.size()) {
      // last node
      OsiNodeSimple node = branchingTree.back();
      branchingTree.pop_back();
      numberNodes++;
      if (node.variable_>=0) {
        // branch - do bounds
        for (i=0;i<numberIntegers;i++) {
          iColumn=which[i];
          setColBounds( iColumn,node.lower_[i],node.upper_[i]);
        }
        // move basis
        setWarmStart(node.basis_);
        // do branching variable
        if (node.way_<0) {
          setColUpper(which[node.variable_],floor(node.value_));
          // now push back node if more to come
          if (node.way_==-1) { 
            node.way_=+2;	  // Swap direction
            branchingTree.push_back(node);
          }
        } else {
          setColLower(which[node.variable_],ceil(node.value_));
          // now push back node if more to come
          if (node.way_==1) { 
            node.way_=-2;	  // Swap direction
            branchingTree.push_back(node);
          }
        }
        // solve
        resolve();
        CoinWarmStart * ws = getWarmStart();
        const CoinWarmStartBasis* wsb =
          dynamic_cast<const CoinWarmStartBasis*>(ws);
        assert (wsb!=NULL); // make sure not volume
        numberIterations += getIterationCount();
        // fix on djs
        int nFixed0=0,nFixed1=0;
        double cutoff;
        getDblParam(OsiDualObjectiveLimit,cutoff);
        double gap=(cutoff-modelPtr_->objectiveValue())*direction+1.0e-4;
        if (gap<1.0e10&&isProvenOptimal()&&!isDualObjectiveLimitReached()) {
          const double * dj = getReducedCost();
          const double * lower = getColLower();
          const double * upper = getColUpper();
          for (i=0;i<numberIntegers;i++) {
            iColumn=which[i];
            if (upper[iColumn]>lower[iColumn]) {
              double djValue = dj[iColumn]*direction;
              if (wsb->getStructStatus(iColumn)==CoinWarmStartBasis::atLowerBound&&
                  djValue>gap) {
                nFixed0++;
                setColUpper(iColumn,lower[iColumn]);
              } else if (wsb->getStructStatus(iColumn)==CoinWarmStartBasis::atUpperBound&&
                         -djValue>gap) {
                nFixed1++;
                setColLower(iColumn,upper[iColumn]);
              }
            }
          }
          //if (nFixed0+nFixed1)
          //printf("%d fixed to lower, %d fixed to upper\n",nFixed0,nFixed1);
        }
        if (!isIterationLimitReached()) {
          OsiNodeSimple newNode(*this,numberIntegers,which,ws);
          // something extra may have been fixed by strong branching
          // if so go round again
          while (newNode.variable_==numberIntegers) {
            resolve();
            newNode = OsiNodeSimple(*this,numberIntegers,which,getWarmStart());
          }
          if (newNode.objectiveValue_<1.0e100) {
            // push on stack
            branchingTree.push_back(newNode);
          }
        } else {
          // maximum iterations - exit
          std::cout<<"Exiting on maximum iterations"
                   <<std::endl;
	  break;
        }
      } else {
        // integer solution - save
        bestNode = node;
        // set cutoff (hard coded tolerance)
        setDblParam(OsiDualObjectiveLimit,(bestNode.objectiveValue_-1.0e-5)*direction);
        std::cout<<"Integer solution of "
                 <<bestNode.objectiveValue_
                 <<" found after "<<numberIterations
                 <<" iterations and "<<numberNodes<<" nodes"
                 <<std::endl;
      }
    }
    std::cout<<"Search took "
             <<numberIterations
             <<" iterations and "<<numberNodes<<" nodes"
             <<std::endl;
    if (bestNode.numberIntegers_) {
      // we have a solution restore
      // do bounds
      for (i=0;i<numberIntegers;i++) {
        iColumn=which[i];
        setColBounds( iColumn,bestNode.lower_[i],bestNode.upper_[i]);
      }
      // move basis
      setWarmStart(bestNode.basis_);
      // set cutoff so will be good (hard coded tolerance)
      setDblParam(OsiDualObjectiveLimit,(bestNode.objectiveValue_+1.0e-5)*direction);
      resolve();
    } else {
      modelPtr_->setProblemStatus(1);
    }
    delete [] which;
  } else {
    std::cout<<"The LP relaxation is infeasible"
             <<std::endl;
    modelPtr_->setProblemStatus(1);
    //throw CoinError("The LP relaxation is infeasible or too expensive",
    //"branchAndBound", "OsiClpSolverInterface");
  }
}
void 
OsiClpSolverInterface::setSpecialOptions(unsigned int value)
{ 
  if ((value&131072)!=0&&(specialOptions_&131072)==0) {
    // Try and keep scaling factors around
    delete baseModel_;
    baseModel_ = new ClpSimplex(*modelPtr_);
    ClpPackedMatrix * clpMatrix = 
      dynamic_cast< ClpPackedMatrix*>(baseModel_->matrix_);
    if (!clpMatrix||clpMatrix->scale(baseModel_)) {
      // switch off again
      delete baseModel_;
      baseModel_=NULL;
      value &= ~131072;
    } else {
      // Off current scaling
      modelPtr_->setRowScale(NULL);
      modelPtr_->setColumnScale(NULL);
      lastNumberRows_=baseModel_->numberRows();
      rowScale_ = CoinDoubleArrayWithLength(2*lastNumberRows_,0);
      int i;
      double * scale;
      double * inverseScale;
      scale = rowScale_.array();
      inverseScale = scale + lastNumberRows_;
      const double * rowScale = baseModel_->rowScale_;
      for (i=0;i<lastNumberRows_;i++) {
	scale[i] = rowScale[i];
	inverseScale[i] = 1.0/scale[i];
      }
      int numberColumns = baseModel_->numberColumns();
      columnScale_ = CoinDoubleArrayWithLength(2*numberColumns,0);
      scale = columnScale_.array();
      inverseScale = scale + numberColumns;
      const double * columnScale = baseModel_->columnScale_;
      for (i=0;i<numberColumns;i++) {
	scale[i] = columnScale[i];
	inverseScale[i] = 1.0/scale[i];
      }
    }
  }
  specialOptions_=value;
  if ((specialOptions_&0x80000000)!=0) {
    // unset top bit if anything set
    if (specialOptions_!=0x80000000) 
      specialOptions_ &= 0x7fffffff;
  }
}
void 
OsiClpSolverInterface::setSpecialOptionsMutable(unsigned int value) const
{ 
  specialOptions_=value;
  if ((specialOptions_&0x80000000)!=0) {
    // unset top bit if anything set
    if (specialOptions_!=0x80000000) 
      specialOptions_ &= 0x7fffffff;
  }
}
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiClpDisasterHandler::OsiClpDisasterHandler (OsiClpSolverInterface * model) 
  : ClpDisasterHandler(),
    osiModel_(model),
    whereFrom_(0),
    phase_(0),
    inTrouble_(false)
{
  if (model)
    setSimplex(model->getModelPtr());
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiClpDisasterHandler::OsiClpDisasterHandler (const OsiClpDisasterHandler & rhs) 
  : ClpDisasterHandler(rhs),
    osiModel_(rhs.osiModel_),
    whereFrom_(rhs.whereFrom_),
    phase_(rhs.phase_),
    inTrouble_(rhs.inTrouble_)
{  
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiClpDisasterHandler::~OsiClpDisasterHandler ()
{
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiClpDisasterHandler &
OsiClpDisasterHandler::operator=(const OsiClpDisasterHandler& rhs)
{
  if (this != &rhs) {
    ClpDisasterHandler::operator=(rhs);
    osiModel_ = rhs.osiModel_;
    whereFrom_ = rhs.whereFrom_;
    phase_ = rhs.phase_;
    inTrouble_ = rhs.inTrouble_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpDisasterHandler * OsiClpDisasterHandler::clone() const
{
  return new OsiClpDisasterHandler(*this);
}

void
OsiClpDisasterHandler::intoSimplex()
{
  inTrouble_=false;
}
bool
OsiClpDisasterHandler::check() const
{
  // Exit if really large number of iterations
  if (model_->numberIterations()> 100000+100*(model_->numberRows()+model_->numberColumns()))
    return true;
  if ((whereFrom_&2)==0||!model_->nonLinearCost()) {
    // dual
    if (model_->numberIterations()<model_->numberRows()+1000) {
      return false;
    } else if (phase_<2) {
      if (model_->numberIterations()> 2*model_->numberRows()+2000||
	  model_->numberDualInfeasibilitiesWithoutFree()||
	  model_->largestDualError()>=1.0e-1) {
#ifdef COIN_DEVELOP
	printf("trouble in phase %d\n",phase_);
#endif
	return true;
      } else {
	return false;
      }
    } else {
      assert (phase_==2);
      if (model_->numberIterations()> 3*model_->numberRows()+2000||
	  model_->largestPrimalError()>=1.0e3) {
#ifdef COIN_DEVELOP
	printf("trouble in phase %d\n",phase_);
#endif
	return true;
      } else {
	return false;
      }
    }
  } else {
    // primal
    if (model_->numberIterations()<model_->numberRows()+4000) {
      return false;
    } else if (phase_<2) {
      if (model_->numberIterations()> 2*model_->numberRows()+2000+
	  model_->numberColumns()/2&&
	  model_->numberDualInfeasibilitiesWithoutFree()>0&&
	  model_->numberPrimalInfeasibilities()>0&&
	  model_->nonLinearCost()->changeInCost()>1.0e8) {
#ifdef COIN_DEVELOP
	printf("trouble in phase %d\n",phase_);
#endif
	return true;
      } else {
	return false;
      }
    } else {
      assert (phase_==2);
      if (model_->numberIterations()> 3*model_->numberRows()+2000||
	  model_->largestPrimalError()>=1.0e3) {
#ifdef COIN_DEVELOP
	printf("trouble in phase %d\n",phase_);
#endif
	return true;
      } else {
	return false;
      }
    }
  }
}
void
OsiClpDisasterHandler::saveInfo()
{
  inTrouble_=true;
}
/* set model. */
void 
OsiClpDisasterHandler::setOsiModel(OsiClpSolverInterface * model)
{
  osiModel_=model;
  model_=model->getModelPtr();
}
// Create C++ lines to get to current state
void 
OsiClpSolverInterface::generateCpp( FILE * fp)
{
  modelPtr_->generateCpp(fp,true);
  // Stuff that can't be done easily
  // setupForRepeatedUse here
  if (!messageHandler()->prefix())
    fprintf(fp,"3  clpModel->messageHandler()->setPrefix(false);\n");
  OsiClpSolverInterface defaultModel;
  OsiClpSolverInterface * other = &defaultModel;
  int iValue1, iValue2;
  double dValue1, dValue2;
  bool takeHint1,takeHint2;
  int add;
  OsiHintStrength strength1,strength2;
  std::string strengthName[] = {"OsiHintIgnore","OsiHintTry","OsiHintDo",
				"OsiForceDo"};
  iValue1 = this->specialOptions();
  iValue2 = other->specialOptions();
  fprintf(fp,"%d  int save_specialOptions = osiclpModel->specialOptions();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setSpecialOptions(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  osiclpModel->setSpecialOptions(save_specialOptions);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->messageHandler()->logLevel();
  iValue2 = other->messageHandler()->logLevel();
  fprintf(fp,"%d  int save_messageHandler = osiclpModel->messageHandler()->logLevel();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->messageHandler()->setLogLevel(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  osiclpModel->messageHandler()->setLogLevel(save_messageHandler);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->cleanupScaling();
  iValue2 = other->cleanupScaling();
  fprintf(fp,"%d  int save_cleanupScaling = osiclpModel->cleanupScaling();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setCleanupScaling(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  osiclpModel->setCleanupScaling(save_cleanupScaling);\n",iValue1==iValue2 ? 7 : 6);
  dValue1 = this->smallestElementInCut();
  dValue2 = other->smallestElementInCut();
  fprintf(fp,"%d  double save_smallestElementInCut = osiclpModel->smallestElementInCut();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setSmallestElementInCut(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  osiclpModel->setSmallestElementInCut(save_smallestElementInCut);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->smallestChangeInCut();
  dValue2 = other->smallestChangeInCut();
  fprintf(fp,"%d  double save_smallestChangeInCut = osiclpModel->smallestChangeInCut();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setSmallestChangeInCut(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  osiclpModel->setSmallestChangeInCut(save_smallestChangeInCut);\n",dValue1==dValue2 ? 7 : 6);
  this->getIntParam(OsiMaxNumIterationHotStart,iValue1);
  other->getIntParam(OsiMaxNumIterationHotStart,iValue2);
  fprintf(fp,"%d  int save_OsiMaxNumIterationHotStart;\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->getIntParam(OsiMaxNumIterationHotStart,save_OsiMaxNumIterationHotStart);\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setIntParam(OsiMaxNumIterationHotStart,%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  osiclpModel->setIntParam(OsiMaxNumIterationHotStart,save_OsiMaxNumIterationHotStart);\n",iValue1==iValue2 ? 7 : 6);
  this->getDblParam(OsiDualObjectiveLimit,dValue1);
  other->getDblParam(OsiDualObjectiveLimit,dValue2);
  fprintf(fp,"%d  double save_OsiDualObjectiveLimit;\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->getDblParam(OsiDualObjectiveLimit,save_OsiDualObjectiveLimit);\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setDblParam(OsiDualObjectiveLimit,%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  osiclpModel->setDblParam(OsiDualObjectiveLimit,save_OsiDualObjectiveLimit);\n",dValue1==dValue2 ? 7 : 6);
  this->getDblParam(OsiPrimalObjectiveLimit,dValue1);
  other->getDblParam(OsiPrimalObjectiveLimit,dValue2);
  fprintf(fp,"%d  double save_OsiPrimalObjectiveLimit;\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->getDblParam(OsiPrimalObjectiveLimit,save_OsiPrimalObjectiveLimit);\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  osiclpModel->setDblParam(OsiPrimalObjectiveLimit,%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  osiclpModel->setDblParam(OsiPrimalObjectiveLimit,save_OsiPrimalObjectiveLimit);\n",dValue1==dValue2 ? 7 : 6);
  this->getHintParam(OsiDoPresolveInInitial,takeHint1,strength1);
  other->getHintParam(OsiDoPresolveInInitial,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoPresolveInInitial;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoPresolveInInitial;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoPresolveInInitial,saveHint_OsiDoPresolveInInitial,saveStrength_OsiDoPresolveInInitial);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoPresolveInInitial,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoPresolveInInitial,saveHint_OsiDoPresolveInInitial,saveStrength_OsiDoPresolveInInitial);\n",add+6);
  this->getHintParam(OsiDoDualInInitial,takeHint1,strength1);
  other->getHintParam(OsiDoDualInInitial,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoDualInInitial;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoDualInInitial;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoDualInInitial,saveHint_OsiDoDualInInitial,saveStrength_OsiDoDualInInitial);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoDualInInitial,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoDualInInitial,saveHint_OsiDoDualInInitial,saveStrength_OsiDoDualInInitial);\n",add+6);
  this->getHintParam(OsiDoPresolveInResolve,takeHint1,strength1);
  other->getHintParam(OsiDoPresolveInResolve,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoPresolveInResolve;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoPresolveInResolve;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoPresolveInResolve,saveHint_OsiDoPresolveInResolve,saveStrength_OsiDoPresolveInResolve);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoPresolveInResolve,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoPresolveInResolve,saveHint_OsiDoPresolveInResolve,saveStrength_OsiDoPresolveInResolve);\n",add+6);
  this->getHintParam(OsiDoDualInResolve,takeHint1,strength1);
  other->getHintParam(OsiDoDualInResolve,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoDualInResolve;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoDualInResolve;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoDualInResolve,saveHint_OsiDoDualInResolve,saveStrength_OsiDoDualInResolve);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoDualInResolve,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoDualInResolve,saveHint_OsiDoDualInResolve,saveStrength_OsiDoDualInResolve);\n",add+6);
  this->getHintParam(OsiDoScale,takeHint1,strength1);
  other->getHintParam(OsiDoScale,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoScale;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoScale;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoScale,saveHint_OsiDoScale,saveStrength_OsiDoScale);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoScale,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoScale,saveHint_OsiDoScale,saveStrength_OsiDoScale);\n",add+6);
  this->getHintParam(OsiDoCrash,takeHint1,strength1);
  other->getHintParam(OsiDoCrash,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoCrash;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoCrash;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoCrash,saveHint_OsiDoCrash,saveStrength_OsiDoCrash);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoCrash,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoCrash,saveHint_OsiDoCrash,saveStrength_OsiDoCrash);\n",add+6);
  this->getHintParam(OsiDoReducePrint,takeHint1,strength1);
  other->getHintParam(OsiDoReducePrint,takeHint2,strength2);
  add = ((takeHint1==takeHint2)&&(strength1==strength2)) ? 1 : 0;
  fprintf(fp,"%d  bool saveHint_OsiDoReducePrint;\n",add+1);
  fprintf(fp,"%d  OsiHintStrength saveStrength_OsiDoReducePrint;\n",add+1);
  fprintf(fp,"%d  osiclpModel->getHintParam(OsiDoReducePrint,saveHint_OsiDoReducePrint,saveStrength_OsiDoReducePrint);\n",add+1);
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoReducePrint,%s,%s);\n",add+3,takeHint1 ? "true" : "false",strengthName[strength1].c_str());
  fprintf(fp,"%d  osiclpModel->setHintParam(OsiDoReducePrint,saveHint_OsiDoReducePrint,saveStrength_OsiDoReducePrint);\n",add+6);
}
