// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// Test individual classes or groups of classes
// This file is licensed under the terms of Eclipse Public License (EPL).

#include "CoinPragma.hpp"

#include "OsiConfig.h"

#include <iostream>
#include <cstdio>

#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiUnitTests.hpp"
#include "OsiTestSolverInterface.hpp"

using namespace OsiUnitTest;

/*
  Currently the Osi unit test is configured to exercise only the external
  solvers. The Osi interfaces for Coin solvers have been moved out to the
  project's repository and each has its own private Osi unit test.

  This unit test will include as many external solvers as are available. If
  none of them are available, the OsiTestSolver (currently a clone of Vol)
  will be used. If any other solver is available, its presence will disable
  use of the test solver. You can disable it manually by undefining
  USETESTSOLVER.

  You may want to use the Osi unit test to compare two or more Coin solvers.
  In particular, OsiSolverInterfaceMpsUnitTest, which runs the Netlib problem
  set, is set up for exactly this sort of comparison.  To take advantage of
  it, you'll need to edit this file and Makefile in order to get it to work.
*/
#define USETESTSOLVER

/*
  Some convenient undef's, to make it easy to isolate a particular solver.
  Uncomment to disable a solver that's included in the build. Leave them
  commented if you're happy with running the unitTest for all solvers in
  the build.
*/
// #undef COIN_HAS_XPR
// #undef COIN_HAS_CPX
// #undef COIN_HAS_GLPK
// #undef COIN_HAS_MSK
// #undef COIN_HAS_GRB
// #undef COIN_HAS_SOPLEX

#ifdef COIN_HAS_XPR
#include "OsiXprSolverInterface.hpp"
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#undef USETESTSOLVER
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef COIN_HAS_GLPK
#include "OsiGlpkSolverInterface.hpp"
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef COIN_HAS_MSK
#include "OsiMskSolverInterface.hpp"
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef COIN_HAS_GRB
#include "OsiGrbSolverInterface.hpp"
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef COIN_HAS_SOPLEX
#include "OsiSpxSolverInterface.hpp"
#ifdef USETESTSOLVER
#undef USETESTSOLVER
#endif
#endif

#ifdef USETESTSOLVER
#include "OsiTestSolverInterface.hpp"
#endif


//----------------------------------------------------------------
// to see parameter list, call unitTest -usage
//----------------------------------------------------------------

int main (int argc, const char *argv[])

{ int totalErrCnt = 0;
  outcomes.clear();

/*
  Start off with various bits of initialisation that don't really belong
  anywhere else.

  First off, synchronise C++ stream i/o with C stdio. This makes debugging
  output a bit more comprehensible. It still suffers from interleave of cout
  (stdout) and cerr (stderr), but -nobuf deals with that.
*/
  std::ios::sync_with_stdio() ;
/*
  Suppress an popup window that Windows shows in response to a crash. See
  note at head of file.
*/
  WindowsErrorPopupBlocker();
/*
  Might as well make use of this convenient Xpress feature.
*/
#ifdef COIN_HAS_XPR
  OsiXprSolverInterface::setLogFileName("xprCallTrace.txt");
#endif

/*
  Process command line parameters.
*/
  std::map<std::string,std::string> parms ;
  if (processParameters(argc,argv,parms) == false)
  { return (1) ; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;

try {
/*
  Test Osi{Row,Col}Cut routines.
*/
#ifdef COIN_HAS_XPR  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiRowCut with OsiXprSolverInterface\n" );
    OsiRowCutUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiColCut with OsiXprSolverInterface\n" );
    OsiColCutUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiXprSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&xprSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_CPX
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiRowCut with OsiCpxSolverInterface\n" );
    OsiRowCutUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiColCut with OsiCpxSolverInterface\n" );
    OsiColCutUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiCpxSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&cpxSi,mpsDir);
  }
#endif

#ifdef USETESTSOLVER
  {
    OsiTestSolverInterface testSi;
    testingMessage( "Testing OsiRowCut with OsiTestSolverInterface\n" );
    OsiRowCutUnitTest(&testSi,mpsDir);
  }
  {
    OsiTestSolverInterface testSi;
    testingMessage( "Testing OsiColCut with OsiTestSolverInterface\n" );
    OsiColCutUnitTest(&testSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_GLPK
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiRowCut with OsiGlpkSolverInterface\n" );
    OsiRowCutUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiColCut with OsiGlpkSolverInterface\n" );
    OsiColCutUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiGlpkSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&glpkSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_MSK  
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiRowCut with OsiMskSolverInterface\n" );
    OsiRowCutUnitTest(&MskSi,mpsDir);
  }
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiColCut with OsiMskSolverInterface\n" );
    OsiColCutUnitTest(&MskSi,mpsDir);
  }
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiMskSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&MskSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_GRB
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiRowCut with OsiGrbSolverInterface\n" );
    OsiRowCutUnitTest(&grbSi,mpsDir);
  }
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiColCut with OsiGrbSolverInterface\n" );
    OsiColCutUnitTest(&grbSi,mpsDir);
  }
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiGrbSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&grbSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_SOPLEX
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiRowCut with OsiSpxSolverInterface\n" );
    OsiRowCutUnitTest(&spxSi,mpsDir);
  }
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiColCut with OsiSpxSolverInterface\n" );
    OsiColCutUnitTest(&spxSi,mpsDir);
  }
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiSpxSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&spxSi,mpsDir);
  }
#endif

  testingMessage( "Testing OsiCuts\n" );
  OsiCutsUnitTest();

/*
  Testing OsiCuts only? A useful option when doing memory access and leak
  checks. Keeps the run time to something reasonable.
*/
  if (parms.find("-cutsOnly") != parms.end()) {
    testingMessage( "Stopping after OsiCuts tests.\n" );
    return 0;
  }

/*
  Run the OsiXXX class test for each solver. It's up to the solver implementor
  to decide whether or not to run OsiSolverInterfaceCommonUnitTest. Arguably
  this should be required.
*/
#ifdef COIN_HAS_XPR
  testingMessage( "Testing OsiXprSolverInterface\n" );
  OsiXprSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_CPX
  testingMessage( "Testing OsiCpxSolverInterface\n" );
  OsiCpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef USETESTSOLVER
  testingMessage( "Testing OsiTestSolverInterface\n" );
  totalErrCnt += OsiTestSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_HAS_GLPK
  testingMessage( "Testing OsiGlpkSolverInterface\n" );
  totalErrCnt += OsiGlpkSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_HAS_MSK
  testingMessage( "Testing OsiMskSolverInterface\n" );
  OsiMskSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_GRB
  testingMessage( "Testing OsiGrbSolverInterface\n" );
  OsiGrbSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_SOPLEX
  testingMessage( "Testing OsiSpxSolverInterface\n" );
  OsiSpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

/*
  Each solver has run its specialised unit test. Check now to see if we need to
  run through the Netlib problems.
*/
  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi;
#   if COIN_HAS_XPR
    OsiSolverInterface * xprSi = new OsiXprSolverInterface;
    vecSi.push_back(xprSi);
#   endif
#   if COIN_HAS_CPX
    OsiSolverInterface * cpxSi = new OsiCpxSolverInterface;
    vecSi.push_back(cpxSi);
#   endif
#   if COIN_HAS_GLPK
    OsiSolverInterface * glpkSi = new OsiGlpkSolverInterface;
    glpkSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry) ;
    glpkSi->setHintParam(OsiDoReducePrint,true,OsiHintDo) ;
    vecSi.push_back(glpkSi);
#   endif
#   if COIN_HAS_MSK
    OsiSolverInterface * MskSi = new OsiMskSolverInterface;
    vecSi.push_back(MskSi);
#   endif
#   if COIN_HAS_GRB
    OsiSolverInterface * grbSi = new OsiGrbSolverInterface;
    vecSi.push_back(grbSi);
#   endif
#   if COIN_HAS_SOPLEX
    OsiSolverInterface * spxSi = new OsiSpxSolverInterface;
    vecSi.push_back(spxSi);
#   endif
#   ifdef USETESTSOLVER
/*
  The test solver is normally Vol, which can't do any of the netlib problems.
  So let's not try.
*/
    {
      std::string solverName ;
      OsiSolverInterface * testSi = new OsiTestSolverInterface;
      testSi->getStrParam(OsiSolverName,solverName) ;
      if (solverName != "vol")
      { vecSi.push_back(testSi); }
      else
      { testingMessage("Test solver vol cannot do Netlib. Skipping.\n") ; }
    }
#   endif

    if (vecSi.size() > 0)
    { testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
      totalErrCnt += OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir); }

    unsigned int i;
    for (i=0; i<vecSi.size(); i++)
      delete vecSi[i];
  }
  else {
    testingMessage( "***Skipped Testing of OsiSolverInterface on Netlib problems***\n" );
    testingMessage( "***use -testOsiSolverInterface to run them.***\n" );
  }
} catch (CoinError& error) {
  std::cout.flush();
  std::cerr << "Caught CoinError exception: ";
  error.print(true);
  totalErrCnt += 1;
} catch (...) {
  std::cout.flush() ;
  std::cerr << "Caught unknown exception." ;
  totalErrCnt += 1 ;
}

/*
  We're done. Report on the results.
*/
  std::cout.flush();
  outcomes.print();

  int nerrors;
  int nerrors_expected;
  outcomes.getCountBySeverity(TestOutcome::ERROR, nerrors, nerrors_expected);

  if (nerrors > nerrors_expected)
    std::cerr << "Tests completed with " << nerrors - nerrors_expected << " unexpected errors." << std::endl ;
  else
    std::cerr << "All tests completed successfully\n";

  return nerrors - nerrors_expected;
}
