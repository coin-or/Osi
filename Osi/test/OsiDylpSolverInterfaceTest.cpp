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
#include <iomanip>
#include <assert.h>
#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpWarmStartBasis.hpp"
#include "OsiDylpMessages.hpp"

namespace {
  char sccsid[] UNUSED = "@(#)OsiDylpSolverInterfaceTest.cpp	1.11	09/25/04" ;
  char cvsid[] UNUSED = "$Id$" ;
}


int test_starts (const std::string& mpsDir)

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

  int retval, errCnt = 0 ;

  CoinRelFltEq eq ;

  double exmip1MinObj = 3.23684210526 ;
  double exmip1MaxObj = 4.5 ;
  std::streamsize old_prec = std::cout.precision() ;

  if (!osi)
  { std::cout
      << "Failed to create first ODSI object." << std::endl ;
    return (++errCnt) ; }
/*
  Read in exmip1 and solve it.
*/
  std::cout << "Boosting verbosity." << std::endl ;
  osi->setHintParam(OsiDoReducePrint,false) ;

  std::string exmpsfile = mpsDir+"exmip1" ;
  std::string probname ;
  std::cout << "Reading mps file \"" << exmpsfile << "\"" << std::endl ;
  retval = osi->readMps(exmpsfile.c_str(), "mps") ;
  if (retval)
  { std::cout
      << "Encountered " << retval
      << " errors reading mps file (1)." << std::endl ;
    return (retval) ; }
  if (osi->getStrParam(OsiProbName,probname) == false)
  { std::cout
      << "Failed to read back problem name." ;
    errCnt++ ; }
  std::cout << "Solving " << probname << " ... " << std::endl ;
  osi->initialSolve() ;
  double val = osi->getObjValue() ;
  if (!eq(val,exmip1MinObj))
  { std::cout
      << "Incorrect objective " << std::setprecision(12) << val
      << "; expecting " << exmip1MinObj
      << ", diff " << val-exmip1MinObj << "." << std::endl ;
    std::cout.precision(old_prec) ;
    errCnt++ ; }
  else
  { std::cout << "And the answer is " << val << "." << std::endl ; }
/*
  Grab a warm start object for later use. No point in continuing if this
  fails.
*/
  std::cout << "Getting a warm start object ... " << std::endl ;
  CoinWarmStart *ws = osi->getWarmStart() ;
  if (!ws)
  { std::cout
      << "Failed to acquire a warm start." << std::endl ;
    return (++errCnt) ; }
/*
  Brief interruption for an idiot check: are the signs of the reduced costs
  correct in the solution, given minimisation? Easy to test with status info
  from the warm start object.
*/
  { const double *cbar = osi->getReducedCost() ;
    std::cout << "Performing sanity test on reduced costs." << std::endl ;
    const CoinWarmStartBasis *wsb =
	dynamic_cast<CoinWarmStartBasis *>(ws) ;
    for (int j = 0 ; j < osi->getNumCols() ; j++)
    { switch (wsb->getStructStatus(j))
      { case CoinWarmStartBasis::atUpperBound:
	{ if (cbar[j] > 0)
	  { std::cout
	      << "Sign error! " << "Var " << j
	      << " at upper bound, cbar = " << cbar[j] << "." << std::endl ;
	    errCnt++ ; }
	  break ; }
        case CoinWarmStartBasis::atLowerBound:
	{ if (cbar[j] < 0)
	  { std::cout
	      << "Sign error! " << "Var " << j
	      << " at lower bound, cbar = " << cbar[j] << "." << std::endl ;
	    errCnt++ ; }
	  break ; }
        case CoinWarmStartBasis::basic:
	{ if (fabs(cbar[j]) > 1.0e-5)
	  { std::cout
	      << "Value error! " << "Var " << j
	      << " is basic, cbar = " << cbar[j]
	      << ", should be zero" << "." << std::endl ;
	    errCnt++ ; }
	  break ; }
	default:
	{ break ; } } } }
/*
  Back to our regular programming. Create an empty warm start object and
  set it as the warm start. Then call resolve(). The call to setWarmStart
  should remove the warm start information in the solver, and the call to
  resolve() should throw.
*/
  { std::cout
      << "Checking behaviour for empty warm start object." << std::endl ;
    std::cout << "Acquiring ... " ;
    CoinWarmStart *emptyWS = osi->getEmptyWarmStart() ;
    if (!emptyWS)
    { std::cout
	<< "Failed to acquire empty warm start." << std::endl ;
      return (++errCnt) ; }
    std::cout << "setting ... " ;
    if (osi->setWarmStart(emptyWS) == false)
    { std::cout
	<< "Failed to install empty warm start." << std::endl ;
      return (++errCnt) ; }
    std::cout << "calling resolve (throw expected) ... " ;
    bool throwSeen = false ;
    try
    { osi->resolve() ; }
    catch (CoinError &ce)
    { std::cout << std::endl << ce.methodName() << ":" << ce.message() ;
      throwSeen = true ; }
    if (throwSeen)
    { std::cout << std::endl << " caught ... success!" << std::endl ; }
    else
    { std::cout << " no throw! ... FAILURE!" ;
      errCnt++ ; }
    delete emptyWS ; }
/*
  Make sure that the warm start information is sufficient (and persistent) by
  discarding the current ODSI object and then installing the warm start
  information in a new ODSI object.
*/
  std::cout << "Discarding current ODSI object ... " << std::endl ;
  delete osi ;
  osi = 0 ;
/*
  We've discarded the first solver. Clone the original warm start object and
  destroy the original.
*/
  std::cout << "Cloning warm start ... " << std::endl ;
  CoinWarmStart *ws_clone = ws->clone() ;
  if (!ws_clone)
  { std::cout
      << "Failed to clone warm start." << std::endl ;
    return (++errCnt) ; }
  delete ws ;
  ws = ws_clone ;
  ws_clone = 0 ;
/*
  Create a second solver, and read in exmip1.  Install the cloned warm start
  in the new solver.
*/
  int level = 5 ;
  level |= 0x10 ;
  std::cout << "Creating new ODSI object ... " << std::endl ;
  osi = new OsiDylpSolverInterface ;
  if (!osi)
  { std::cout
      << "Failed to create second ODSI object." << std::endl ;
    return (++errCnt) ; }

  osi->setHintParam(OsiDoReducePrint,false,OsiForceDo,&level) ;
  osi->getHintParam(OsiDoReducePrint,sense,strength,p_info) ;
  std::cout << "Verbosity now maxed at "
	    << *reinterpret_cast<int *>(p_info) << "." << std::endl ;

  retval = osi->readMps(exmpsfile.c_str(), "mps") ;
  if (retval)
  { std::cout
      << "Encountered " << retval
      << " errors reading mps file (2)." << std::endl ;
    return (retval) ; }
  std::cout << "Installing cloned warm start object ... " << std::endl ;
  if (osi->setWarmStart(ws) == false)
  { std::cout
      << "Failed to install valid warm start after deleting original solver."
      << std::endl ;
    return (++errCnt) ; }
/*
  Resolve. Check that we did not pivot (much) and that the objective hasn't
  changed. Set the print level quite high (we need to do this at some
  point).
*/
  std::cout << "Resolving the lp ... " << std::endl ;
  osi->resolve() ;
  val = osi->getObjValue() ;
  int pivots = osi->getIterationCount() ;
  if (!eq(val,exmip1MinObj))
  { std::cout
      << "Incorrect objective " << std::setprecision(12) << val
      << "; expecting " << exmip1MinObj
      << ", diff " << val-exmip1MinObj << "." << std::endl ;
    std::cout.precision(old_prec) ;
    errCnt++ ; }
  else
  if (pivots > 1)
  { std::cout
      << "Excessive pivots; counted "
      << pivots << ", expected <= 1." << std::endl ;
    errCnt++ ; }
  else
  { std::cout
      << std::endl << "And the answer is " << val << " after "
      << pivots << " pivots." << std::endl ; }
  delete ws ;
  ws = 0 ;
/*
  Flip the objective and do a hot start.
*/
  osi->setHintParam(OsiDoReducePrint,true,OsiForceDo) ;
  std::cout << "Reducing verbosity." << std::endl ;
  std::cout << "Changing objective sense to maximisation ..." ;
  osi->setObjSense(-1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->markHotStart() ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  if (!eq(val,exmip1MaxObj))
  { std::cout
      << "Incorrect objective " << std::setprecision(12) << val
      << "; expecting " << exmip1MaxObj
      << ", diff " << val-exmip1MaxObj << "." << std::endl ;
    std::cout.precision(old_prec) ;
    errCnt++ ; }
  else
  { std::cout
      << std::endl << "And the answer is " << val << "." << std::endl ; }
/*
  Another brief interruption for an idiot check: are the signs of the reduced
  costs correct in the solution, given maximisation?
*/
  { const double *cbar = osi->getReducedCost() ;
    ws = osi->getWarmStart() ;
    const OsiDylpWarmStartBasis *odsi_wsb =
	dynamic_cast<OsiDylpWarmStartBasis *>(ws) ;
    std::cout
      << "Performing sanity test on reduced costs." << std::endl ;
    for (int j = 0 ; j < osi->getNumCols() ; j++)
    { switch (odsi_wsb->getStructStatus(j))
      { case CoinWarmStartBasis::atUpperBound:
	{ if (cbar[j] < 0)
	  { std::cout
	      << "Sign error! " << "Var " << j
	      << " at upper bound, cbar = " << cbar[j] << "." << std::endl ;
	    errCnt++ ; }
	  break ; }
        case CoinWarmStartBasis::atLowerBound:
	{ if (cbar[j] > 0)
	  { std::cout
	      << "Sign error! " << "Var " << j
	      << " at lower bound, cbar = " << cbar[j] << "." << std::endl ;
	    errCnt++ ; }
	  break ; }
        case CoinWarmStartBasis::basic:
	{ if (fabs(cbar[j]) > 1.0e-5)
	  { std::cout
	      << "Value error! " << "Var " << j
	      << " is basic, cbar = " << cbar[j] << ", should be zero"
	      << "." << std::endl ;
	    errCnt++ ; }
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
	    << *reinterpret_cast<int *>(p_info) << "." << std::endl ;
/*
  And return to minimisation.
*/
  std::cout << "And back to minimisation ..." ;
  osi->setObjSense(1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  if (!eq(val,exmip1MinObj))
  { std::cout
      << "Incorrect objective " << std::setprecision(12) << val
      << "; expecting " << exmip1MinObj
      << ", diff " << val-exmip1MinObj << "." << std::endl ;
    std::cout.precision(old_prec) ;
    errCnt++ ; }
  else
  { std::cout
      << std::endl << "And the answer is " << val << "." << std::endl ; }

  delete osi ;

  return (errCnt) ; }


/*! OsiDylp unit test driver

  This is the unit test routine for OsiDylpSolverInterface. It tests for
  problems that have been uncovered and fixed already. If it fails, you've
  probably tickled a new bug. Please file a bug report.
*/

int OsiDylpSolverInterfaceUnitTest (const std::string &mpsDir,
				     const std::string &netLibDir)

{ int errCnt = 0 ; 

  std::cout << "Starting dylp OSI interface tests ... " << std::endl ;
  OsiDylpSolverInterface* osi = new OsiDylpSolverInterface ;
  osi->handler_->setLogLevel(3) ;
  osi->handler_->message(ODSI_TEST_MSG,osi->messages_) ;
  osi->newLanguage(CoinMessages::uk_en) ;
  osi->handler_->message(ODSI_TEST_MSG,osi->messages_) ;
  osi->handler_->finish() ;
  std::cout <<
    "Calling OsiSolverInterfaceCommonUnitTest for basic tests ... "
    << std::endl ;
  errCnt += OsiSolverInterfaceCommonUnitTest(osi,mpsDir,netLibDir) ;
  if (errCnt != 0)
  { std::cout
      << "ODSIUnitTest: " << errCnt << " errors after common unit test."
      << std::endl ; }
/*
  Test the reset function.
*/
  std::cout << "Testing reset ... " << std::endl ;
  OsiDylpSolverInterface* osi2 = new OsiDylpSolverInterface ;
  osi->reset() ;
# ifndef _MSC_VER
  osi->assert_same(*osi,*osi2,true) ;
# endif
  delete osi ;
  delete osi2 ;
  std::cout <<
    "Testing cold/warm/hot start ... " << std::endl ;
  errCnt += test_starts(mpsDir) ;

  std::cout
    << std::endl << " dylp tests completed, "
    << errCnt << " errors." << std::endl << std::endl ;

  return (errCnt) ; }

#endif // COIN_HAS_DYLP

