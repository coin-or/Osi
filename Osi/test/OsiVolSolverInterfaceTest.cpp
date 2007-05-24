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
void
OsiVolSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

  // Do common solverInterface testing
  {
    OsiVolSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }

}
