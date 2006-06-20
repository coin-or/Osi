/*! \legal
  Copyright (C) 2002, 2003, 2004.
  Lou Hafer, Stephen Tse, International Business Machines Corporation and
  others. All Rights Reserved.

  This file is a portion of the COIN/OSI interface for dylp.
*/

#include "OsiConfig.h"

#ifdef COIN_HAS_DYLP

#if defined(_MSC_VER)

/* Turn off compiler warning about long names */
#  pragma warning(disable:4786)

#endif

/*
  Unset NDEBUG, if it happens to be set. This code needs to be rewritten if
  it's going to run without assert.
*/
#ifdef NDEBUG
#undef NDEBUG
#endif

/* 
  Rudimentary tests for the Dylp OSI layer. The assumption is that
  OsiSolverInterfaceCommonUnitTest checks basic functionality. The routine
  test_starts does a cursory check of cold/warm/hot starts.

  These tests need to be sharpened considerably, if they are to have any
  teeth.
*/

#include <iostream>
#include <assert.h>
#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpWarmStartBasis.hpp"
#include "OsiDylpMessages.hpp"

namespace {
  char sccsid[] UNUSED = "@(#)OsiDylpSolverInterfaceTest.cpp	1.11	09/25/04" ;
  char cvsid[] UNUSED = "$Id$" ;
}


void test_starts (const std::string& mpsDir)
/*
  This routine makes a number of checks for warm and hot start capabilities.
    * Create and attempt to set an empty warm start object.
    * Create an ODSI object and solve the exmip1 sample MPS file with
      initialSolve.
    * Get a warm start object, then destroy the ODSI object. Create a new ODSI
      object, clone the saved warm start, install it, and resolve. Check that
      the objective is the same and that we did not pivot.
    * Change the objective sense and resolve from hot start.
*/
{ OsiDylpSolverInterface *osi = new OsiDylpSolverInterface ;
  OsiHintStrength strength ;
  bool sense ;
  void *p_info ;

/*
  Read in exmip1 and solve it.
*/
  std::cout << "Boosting verbosity.\n" ;
  osi->setHintParam(OsiDoReducePrint,false) ;

  std::string exmpsfile = mpsDir+"exmip1" ;
  std::string probname ;
  std::cout << "Reading mps file \"" << exmpsfile << "\"\n" ;
  osi->readMps(exmpsfile.c_str(), "mps") ;
  assert(osi->getStrParam(OsiProbName,probname)) ;
  std::cout << "Solving " << probname << " ... \n" ;
  osi->initialSolve() ;
  double val = osi->getObjValue() ;
  std::cout << "And the answer is " << val << ".\n" ;
  assert(fabs(val - 3.23) < 0.01) ;
/*
  Grab a warm start object for later use.
*/
  std::cout << "Getting a warm start object ... \n" ;
  CoinWarmStart *ws = osi->getWarmStart() ;
  assert(ws) ;
/*
  Brief interruption for an idiot check: are the signs of the reduced costs
  correct in the solution, given minimisation? Easy to test with status info
  from the warm start object.
*/
  { const double *cbar = osi->getReducedCost() ;
    std::cout << "Performing sanity test on reduced costs.\n" ;
    const CoinWarmStartBasis *wsb =
	dynamic_cast<CoinWarmStartBasis *>(ws) ;
    for (int j = 0 ; j < osi->getNumCols() ; j++)
    { switch (wsb->getStructStatus(j))
      { case CoinWarmStartBasis::atUpperBound:
	{ if (cbar[j] > 0)
	  { std::cout << "Sign error! " ;
	    std::cout << "Var " << j
		      << " at upper bound, cbar = " << cbar[j] << ".\n" ; }
	  break ; }
        case CoinWarmStartBasis::atLowerBound:
	{ if (cbar[j] < 0)
	  { std::cout << "Sign error! " ;
	    std::cout << "Var " << j
		      << " at lower bound, cbar = " << cbar[j] << ".\n" ; }
	  break ; }
        case CoinWarmStartBasis::basic:
	{ if (fabs(cbar[j]) > 1.0e-5)
	  { std::cout << "Value error! " ;
	    std::cout << "Var " << j << " is basic, cbar = " << cbar[j] ;
	    std::cout << ", should be zero" ;
	    std::cout << ".\n" ; }
	  break ; }
	default:
	{ break ; } } } }
/*
  Back to our regular programming. Create an empty warm start object and
  set it as the warm start. Then call resolve(). The call to setWarmStart
  should remove the warm start information in the solver, and the call to
  resolve() should throw.
*/
  { std::cout << "Checking behaviour for empty warm start object.\n" ;
    std::cout << "Acquiring ... " ;
    CoinWarmStart *ws = osi->getEmptyWarmStart() ;
    assert(ws) ;
    std::cout << "setting ... " ;
    assert(osi->setWarmStart(ws)) ;
    std::cout << "calling resolve (throw expected) ... " ;
    bool throwSeen = false ;
    try
    { osi->resolve() ; }
    catch (CoinError &ce)
    { std::cout << "\n" << ce.methodName() << ":" << ce.message() ;
      throwSeen = true ; }
    if (throwSeen)
    { std::cout << "\n caught ... success!\n" ; }
    else
    { std::cout << " no throw! ... FAILURE!" ; }
    delete ws ; }
/*
  Back to our regular programming.
*/
  std::cout << "Discarding current ODSI object ... \n" ;
  delete osi ;
/*
  We've discarded the first solver. Create a second solver, and read in the
  problem. Clone the original warm start object and destroy the original.
  Install the clone in the new solver, and resolve. Check that we did not
  pivot and that the objective hasn't changed. Set the print level quite high.
*/
  std::cout << "Creating new ODSI object ... \n" ;
  osi = new OsiDylpSolverInterface ;
  assert(osi) ;
  std::cout << "Testing anonymous clone for warm start ... \n" ;
  CoinWarmStart *ws_clone = ws->clone() ;
  assert(ws_clone) ;
  delete ws ;
  ws = ws_clone ;
  ws_clone = 0 ;

  int level = 5 ;
  level |= 0x10 ;
  osi->setHintParam(OsiDoReducePrint,true,OsiForceDo,&level) ;
  osi->getHintParam(OsiDoReducePrint,sense,strength,p_info) ;
  std::cout << "Verbosity now maxed at "
	    << *reinterpret_cast<int *>(p_info) << ".\n" ;

  osi->readMps(exmpsfile.c_str(), "mps") ;

  std::cout << "Installing cloned warm start object ... \n" ;
  assert(osi->setWarmStart(ws)) ;
  std::cout << "Resolving the lp ... \n" ;

  osi->resolve() ;
  val = osi->getObjValue() ;
  int pivots = osi->getIterationCount() ;
  std::cout << "\nAnd the answer is " << val << " after " <<
	       pivots << " pivots.\n" ;
  std::cout << "(Expecting 3.23684 with 0 <= pivots <= 1.)\n" ;
  delete ws ;
  assert(fabs(val - 3.23) < 0.01 && pivots <= 1) ;

  osi->setHintParam(OsiDoReducePrint,true,OsiForceDo) ;
  std::cout << "Reducing verbosity.\n" ;
  std::cout << "Changing objective sense ..." ;
  osi->setObjSense(-1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->markHotStart() ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  std::cout << "\nAnd the answer is " << val
	    << " (expecting " << 4.5 << ").\n" ;
/*
  Brief interruption for an idiot check again: are the signs of the reduced
  costs correct in the solution, given maximisation?
*/
  { const double *cbar = osi->getReducedCost() ;
    ws = osi->getWarmStart() ;
    const OsiDylpWarmStartBasis *odsi_wsb =
	dynamic_cast<OsiDylpWarmStartBasis *>(ws) ;
    std::cout << "Performing sanity test on reduced costs.\n" ;
    for (int j = 0 ; j < osi->getNumCols() ; j++)
    { switch (odsi_wsb->getStructStatus(j))
      { case CoinWarmStartBasis::atUpperBound:
	{ if (cbar[j] < 0)
	  { std::cout << "Sign error! " ;
	    std::cout << "Var " << j
		      << " at upper bound, cbar = " << cbar[j] << ".\n" ; }
	  break ; }
        case CoinWarmStartBasis::atLowerBound:
	{ if (cbar[j] > 0)
	  { std::cout << "Sign error! " ;
	    std::cout << "Var " << j
		      << " at lower bound, cbar = " << cbar[j] << ".\n" ; }
	  break ; }
        case CoinWarmStartBasis::basic:
	{ if (fabs(cbar[j]) > 1.0e-5)
	  { std::cout << "Value error! " ;
	    std::cout << "Var " << j << " is basic, cbar = " << cbar[j] ;
	    std::cout << ", should be zero" ;
	    std::cout << ".\n" ; }
	  break ; }
	default:
	{ break ; } } }
    delete ws ; }
/*
  Turn off printing, to make sure we can get dylp to shut up.
*/
  level = 0 ;
  osi->setHintParam(OsiDoReducePrint,true,OsiForceDo,&level) ;
  osi->getHintParam(OsiDoReducePrint,sense,strength,p_info) ;
  std::cout << "Verbosity now at "
	    << *reinterpret_cast<int *>(p_info) << ".\n" ;

  std::cout << "And back ..." ;
  osi->setObjSense(1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  std::cout << "\nAnd the answer is " << val
	    << " (expecting " << 3.23684 << ").\n" ;
  assert(fabs(val - 3.23) < 0.01) ;

  delete osi ;

  return ; }


/*!
  This is the unit test routine for OsiDylpSolverInterface. It tests for
  problems that have been uncovered and fixed already. If it fails, you've
  probably tickled a new bug. Please file a bug report.
*/

void OsiDylpSolverInterfaceUnitTest (const std::string &mpsDir,
				     const std::string &netLibDir)
/*
  Dylp unit test driver.
*/
{ std::cout << "Starting dylp OSI interface tests ...\n" ;
  OsiDylpSolverInterface* osi = new OsiDylpSolverInterface ;
  osi->handler_->setLogLevel(3) ;
  osi->handler_->message(ODSI_TEST_MSG,osi->messages_) ;
  osi->newLanguage(CoinMessages::uk_en) ;
  osi->handler_->message(ODSI_TEST_MSG,osi->messages_) ;
  osi->handler_->finish() ;
  std::cout <<
    "Calling OsiSolverInterfaceCommonUnitTest for basic tests ...\n" ;
  OsiSolverInterfaceCommonUnitTest(osi,mpsDir,netLibDir) ;
/*
  Test the reset function.
*/
  std::cout << "Testing reset ...\n" ;
  OsiDylpSolverInterface* osi2 = new OsiDylpSolverInterface ;
  osi->reset() ;
  osi->assert_same(*osi,*osi2,true) ;
  delete osi ;
  delete osi2 ;
  std::cout <<
    "Testing cold/warm/hot start ...\n" ;
  test_starts(mpsDir) ;

  std::cout << "\n dylp tests completed.\n\n" ;

  return ; }

#endif // COIN_HAS_DYLP

