// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "OsiVolSolverInterface.hpp"

//#############################################################################

//--------------------------------------------------------------------------
// test EKKsolution methods.
int
OsiVolSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

  // Do common solverInterface testing
  {
    OsiVolSolverInterface m;
    return OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }

}
