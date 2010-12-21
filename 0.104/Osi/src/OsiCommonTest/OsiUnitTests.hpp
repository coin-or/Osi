// Copyright (C) 2010
// All Rights Reserved.

#ifndef OSISOLVERINTERFACETEST_HPP_
#define OSISOLVERINTERFACETEST_HPP_

#include <string>
#include <vector>

class OsiSolverInterface;
class CoinPackedVectorBase;

/** \brief Run solvers on NetLib problems.

  The routine creates a vector of NetLib problems (problem name, objective,
  various other characteristics), and a vector of solvers to be tested.

  Each solver is run on each problem. The run is deemed successful if the
  solver reports the correct problem size after loading and returns the
  correct objective value after optimization.

  If multiple solvers are available, the results are compared pairwise against
  the results reported by adjacent solvers in the solver vector. Due to
  limitations of the volume solver, it must be the last solver in vecEmptySiP.
*/
int OsiSolverInterfaceMpsUnitTest
  (const std::vector<OsiSolverInterface*> & vecEmptySiP,
   const std::string& mpsDir);

int OsiSolverInterfaceCommonUnitTest
  (const OsiSolverInterface* emptySi,
   const std::string& mpsDir,
   const std::string& netlibDir);

void OsiCutsUnitTest();

void OsiColCutUnitTest(const OsiSolverInterface * baseSiP, const std::string & mpsDir);

void OsiRowCutUnitTest(const OsiSolverInterface * baseSiP, const std::string & mpsDir);

void OsiRowCutDebuggerUnitTest(const OsiSolverInterface * baseSiP, const std::string & mpsDir);

int testSimplexAPI(const OsiSolverInterface* emptySi,
		   const std::string& mpsDir) ;

namespace OsiUnitTest {

void failureMessage(const std::string &solverName,
		    const std::string &message) ;
void failureMessage(const OsiSolverInterface &si,
		    const std::string &message) ;
void testingMessage(const char *const msg) ;

bool equivalentVectors(const OsiSolverInterface * si1,
		       const OsiSolverInterface * si2,
		       double tol,
		       const double * v1,
		       const double * v2,
		       int size) ;

bool compareProblems(OsiSolverInterface *osi1, OsiSolverInterface *osi2) ;

bool isEquivalent(const CoinPackedVectorBase &pv, int n, const double *fv) ;

} // end namespace OsiUnitTest

#endif /*OSISOLVERINTERFACETEST_HPP_*/
