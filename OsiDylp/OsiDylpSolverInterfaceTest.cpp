//-----------------------------------------------------------------------------
// Copyright (C) 2002, Lou Hafer, Stephen Tse, International Business Machines
// Corporation and others.  All Rights Reserved.
//-----------------------------------------------------------------------------


#ifdef COIN_USE_DYLP

#if defined(_MSC_VER)

/* Turn off compiler warning about long names */
#  pragma warning(disable:4786)

#endif

/* 
  Rudimentary tests for Dylp OSI layer. The assumption is that
  OsiSolverInterfaceCommonUnitTest checks basic functionality. The routine
  test_starts does a cursory check of cold/warm/hot starts.

  These tests need to be sharpened considerably, if they are to have any
  teeth.
*/

#include <iostream>
#include <assert.h>
#include <OsiWarmStartBasis.hpp>
#include "OsiDylpSolverInterface.hpp"


static char sccsid[] = "%W%	%G%" ;


void test_starts (const std::string& mpsDir)
/*
  Some basic tests: Solve the exmip1 sample MPS file with the initialSolve
  routine, then check warm and hot starts. Return the ODSI object for
  further tests.
*/
{ OsiDylpSolverInterface *osi = new OsiDylpSolverInterface ;
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

  std::cout << "Getting a warm start ... \n" ;
  OsiWarmStart *wsb = osi->getWarmStart() ;
  std::cout << "Setting a warm start ... \n" ;
  osi->setWarmStart(wsb) ;
  std::cout << "Resolving the lp ... \n" ;
  osi->resolve() ;
  val = osi->getObjValue() ;
  std::cout << "And the answer is " << val << ".\n" ;
  delete wsb ;
  assert(fabs(val - 3.23) < 0.01) ;

  std::cout << "Changing objective sense ..." ;
  osi->setObjSense(-1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->markHotStart() ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  std::cout << "And the answer is " << val << ".\n" ;

  std::cout << "And back ..." ;
  osi->setObjSense(1.0) ;
  std::cout << "Attempting hot start ..." ;
  osi->solveFromHotStart() ;
  val = osi->getObjValue() ;
  std::cout << "And the answer is " << val << ".\n" ;
  assert(fabs(val - 3.23) < 0.01) ;

  delete osi ;

  return ; } ;


void OsiDylpSolverInterfaceUnitTest (const std::string& mpsDir)
/*
  Dylp unit test driver.
*/
{ std::cout << "Starting dylp OSI interface tests ...\n" ;
  std::cout <<
    "Calling OsiSolverInterfaceCommonUnitTest for basic tests ...\n" ;
  OsiDylpSolverInterface* osi = new OsiDylpSolverInterface ;
  OsiSolverInterfaceCommonUnitTest(osi,mpsDir) ;
  delete osi ;
  std::cout <<
    "Testing cold/warm/hot start ...\n" ;
  test_starts(mpsDir) ;
  std::cout << "\n dylp tests completed." ;

  return ; }

#endif // COIN_USE_DYLP

