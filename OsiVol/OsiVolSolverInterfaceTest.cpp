// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifdef COIN_USE_VOL

#include <cassert>

#include "OsiVolSolverInterface.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif

//--------------------------------------------------------------------------
// test EKKsolution methods.
void
OsiVolSolverInterfaceUnitTest(const std::string & mpsDir)
{

  // Do common solverInterface testing
  {
    OsiVolSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir);
  }

}
#endif
