// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "OsiTestSolverInterface.hpp"
#include "OsiUnitTests.hpp"

//#############################################################################

//--------------------------------------------------------------------------
// test EKKsolution methods.
int
OsiTestSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

  // Do common solverInterface testing
  {
    OsiTestSolverInterface m;
    return OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }

}
