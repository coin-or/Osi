// Example of using COIN-OR OSI
// including accessing solver-specific functions

#include <iostream>
#include "OsiClpSolverInterface.hpp"

int
main(void)
{
   // Create a problem pointer.  We use the base class here.
   OsiSolverInterface *si;

   // When we instantiate the object, we need a specific derived class.
   si = new OsiClpSolverInterface;  

   // The next few lines are solver-dependent!
   ClpSimplex * clpPointer;
   clpPointer = (dynamic_cast<OsiClpSolverInterface *>(si))->getModelPtr();

   clpPointer->setLogLevel(0);
   //clpPointer->setMaximumIterations(10);
   // Could tell Clp many other things

   // Read in an mps file.  This one's from the MIPLIB library.
   si->readMps("../../Data/Sample/p0033");

   // Solve the (relaxation of the) problem
   si->initialSolve();

   // Check the solution
   if ( si->isProvenOptimal() ) { 
      std::cout << "Found optimal solution!" << std::endl; 
      std::cout << "Objective value is " << si->getObjValue() << std::endl;

      int n = si->getNumCols();
      const double *solution;
      solution = si->getColSolution();
      // We could then print the solution or examine it.
   } else {
      std::cout << "Didn't find optimal solution." << std::endl;
      // Could then check other status functions.
   }

  return 0;
}
