// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
 
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <cstdio>

#include "OsiSolverInterface.hpp"
#include "OsiSimplexInterface.hpp"
//#############################################################################

void
OsiSimplexInterfaceCommonUnitTest(const OsiSolverInterface* emptySi,
				 const std::string & mpsDir)
{
  OsiSolverInterface * si = emptySi->clone();
  std::string solverName;
  si->getStrParam(OsiSolverName,solverName);
  OsiSimplexInterface * simplexSi = dynamic_cast<OsiSimplexInterface *> (si);
  
  if (simplexSi) {
    // solve an lp by hand
    
    std::string fn = mpsDir+"p0033";
    si->readMps(fn.c_str(),"mps");
    si->setObjSense(-1.0);
    si->initialSolve();
    si->setObjSense(1.0);
    // enable special mode
    simplexSi->enableSimplexInterface(true);
    // we happen to know that variables are 0-1 and rows are L
    int numberIterations=0;
    int numberColumns = si->getNumCols();
    int numberRows = si->getNumRows();
    double * fakeCost = new double[numberColumns];
    double * duals = new double [numberRows];
    double * djs = new double [numberColumns];
    const double * solution = si->getColSolution();
    memcpy(fakeCost,si->getObjCoefficients(),numberColumns*sizeof(double));
    while (1) {
      const double * dj;
      const double * dual;
      if ((numberIterations&1)==0) {
	// use given ones
	dj = si->getReducedCost();
	dual = si->getRowPrice();
      } else {
	// create
	dj = djs;
	dual = duals;
	simplexSi->getReducedGradient(djs,duals,fakeCost);
      }
      int i;
      int colIn=9999;
      int direction=1;
      double best=1.0e-6;
      // find most negative reduced cost
      // Should check basic - but should be okay on this problem
      for (i=0;i<numberRows;i++) {
	double value=dual[i];
	if (value>best) {
	  direction=-1;
	  best=value;
	  colIn=-i-1;
	}
      }
      for (i=0;i<numberColumns;i++) {
	double value=dj[i];
	if (value<-best&&solution[i]<1.0e-6) {
	  direction=1;
	  best=-value;
	  colIn=i;
	} else if (value>best&&solution[i]>1.0-1.0e-6) {
	  direction=-1;
	  best=value;
	  colIn=i;
	}
      }
      if (colIn==9999)
	break; // should be optimal
      int colOut;
      int outStatus;
      double theta;
      assert(!simplexSi->primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL));
      printf("out %d, direction %d theta %g\n",
	     colOut,outStatus,theta);
      numberIterations++;
    }
    delete [] fakeCost;
    delete [] duals;
    delete [] djs;
    // exit special mode
    simplexSi->disableSimplexInterface();
    si->resolve();
    assert (!si->getIterationCount());
    si->setObjSense(-1.0);
    si->initialSolve();
    std::cout<<solverName<<" passed OsiSimplexInterface test"<<std::endl;
  } else {
    std::cout<<solverName<<" has no OsiSimplexInterface"<<std::endl;
  }
  delete si;
}

