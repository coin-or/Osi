/*
  Copyright (C) 2000 -- 2010, Lou Hafer, International Business Machines,
  and others.  All Rights Reserved.
*/

#include "CoinPragma.hpp"

#include "OsiUnitTests.hpp"

#include "OsiConfig.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
/*
  #include "CoinTime.hpp"
  #include <cstdlib>
  #include <cassert>
  #include <vector>
  #include <iostream>
  #include <iomanip>
  #include <sstream>
  #include <cstdio>
*/
/*
  A utility definition which allows for easy suppression of unused variable
  warnings from GCC. Handy in this environment, where we're constantly def'ing
  things in and out.
*/
#ifndef UNUSED
# if defined(__GNUC__)
#   define UNUSED __attribute__((unused))
# else
#   define UNUSED
# endif
#endif

#include "OsiSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
/*
  #include "CoinFloatEqual.hpp"
  #include "CoinPackedVector.hpp"
  #include "CoinWarmStartBasis.hpp"
  #include "OsiRowCut.hpp"
  #include "OsiCuts.hpp"
  #include "OsiPresolve.hpp"
*/


namespace OsiUnitTest {

//#############################################################################
// Helper routines for messages.
//#############################################################################

// A helper function to write out a message about a test failure
void failureMessage( const std::string & solverName,
		     const std::string & message )
{
  std::string messageText;
  messageText = "*** ";
  messageText += solverName + "SolverInterface testing issue: ";
  messageText += message;
  // flush stdout so that error messages are properly interleaved.
  std::cout.flush() ;
  std::cerr << messageText.c_str() << std::endl;
}

void failureMessage( const OsiSolverInterface & si,
		     const std::string & message )
{
  std::string solverName;
  si.getStrParam(OsiSolverName,solverName);
  failureMessage(solverName,message);
}

/*
  Display message on stderr. Flush cout buffer before printing the message,
  so that output comes out in order in spite of buffered cout.
*/
void testingMessage( const char * const msg )
{
  std::cout.flush() ;
  std::cerr << msg;
}
//#############################################################################
// Vector comparison utility.
//#############################################################################

// A helper function to compare the equivalence of two vectors
bool equivalentVectors (const OsiSolverInterface * si1,
			const OsiSolverInterface * si2,
			double tol,
			const double * v1,
			const double * v2,
			int size)
{
  bool retVal = true;
  CoinRelFltEq eq(tol) ;
  int i;
  for ( i=0; i<size; i++ ) {
    if ( !eq(v1[i],v2[i]) ) {
      std::cout.flush() ;
      std::cerr <<"eq " <<i <<" " <<v1[i] <<" " <<v2[i] <<std::endl;
      retVal = false;
      break;
    }
  }
  return retVal;
}

/*
  Check a packed vector for equality with a full vector. The algorithm is to
  first confirm that the elements of the packed vector are present in the full
  vector, then scan the full vector to make sure there are no additional
  nonzeros.
*/
bool isEquivalent (const CoinPackedVectorBase &pv, int n, const double *fv)
{ 
  int pvCnt = pv.getNumElements() ;
  const int *indices = pv.getIndices() ;
  const double *elems = pv.getElements() ;

  bool retval = true ;
  CoinRelFltEq eq ;

  for (int v = 0 ; v < pvCnt ; v++)
  { int k = indices[v] ;
    if (!eq(elems[v],fv[k]))
    { retval = false ;
      break ; } }
  if (retval == true)
  { int fvCnt = 0 ;
    for (int k = 0 ; k < n ; k++)
    { if (!eq(fv[k],0.0))
        fvCnt++ ; }
    if (fvCnt != pvCnt)
      retval = false ; }

  return (retval) ;
}


/*
  Method to compare the problem representation held by a pair of solver
  interfaces.
*/
bool compareProblems (OsiSolverInterface *osi1, OsiSolverInterface *osi2) {

  bool areEquiv = true ;
  std::string si1Name, si2Name ;
  osi1->getStrParam(OsiSolverName,si1Name) ;
  osi2->getStrParam(OsiSolverName,si2Name) ;

  // Compare row and column counts
  int colCnt = 0 ;
  if (osi1->getNumCols() != osi2->getNumCols())
  { std::cerr
      << "  Unequal column count, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }
  else
  { colCnt = osi1->getNumCols() ; }

  int rowCnt = 0 ;
  if (osi1->getNumRows() != osi2->getNumRows())
  { std::cerr
      << "  Unequal row count, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }
  else
  { rowCnt = osi1->getNumRows() ; }

  // Compare column bounds
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
  		osi1->getColLower(),osi2->getColLower(),colCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal column lower bounds, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
  		osi1->getColUpper(),osi2->getColUpper(),colCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal column upper bounds, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare row bounds
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
		osi1->getRowLower(),osi2->getRowLower(),rowCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal row lower bounds, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
		osi1->getRowUpper(),osi2->getRowUpper(),rowCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal row lower bounds, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare row sense
  { const char *rowSense1 = osi1->getRowSense() ;
    const char *rowSense2 = osi2->getRowSense() ;
    areEquiv = true ;
    for (int r = 0 ; r < rowCnt && areEquiv == true ; r++)
    { if (rowSense1[r] != rowSense2[r])
      { areEquiv = false ; } }
    if (areEquiv == false)
    { std::cerr
	<< "  Unequal row sense, "
	<< si1Name << " vs. " << si2Name << std::endl ;
      return (false) ; } }

  // Compare row rhs
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
  		osi1->getRightHandSide(),osi2->getRightHandSide(),rowCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal right-hand-side, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare range
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
  		osi1->getRowRange(),osi2->getRowRange(),rowCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal row range, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare objective sense
  if (osi1->getObjSense() != osi2->getObjSense())
  { std::cerr
      << "  Unequal objective sense, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare objective coefficients
  areEquiv = equivalentVectors(osi1,osi2,1.e-10,
  		osi1->getObjCoefficients(),osi2->getObjCoefficients(),colCnt) ;
  if (areEquiv == false)
  { std::cerr
      << "  Unequal objective coefficients, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare number of elements
  if (osi1->getNumElements() != osi2->getNumElements())
  { std::cerr
      << "  Unequal number of constraint matrix coefficients, "
      << si1Name << " vs. " << si2Name << std::endl ;
    return (false) ; }

  // Compare constraint matrix, for both row-major and column-major orderings
  { const CoinPackedMatrix *rmm1=osi1->getMatrixByRow() ;
    const CoinPackedMatrix *rm  =osi2->getMatrixByRow() ;
    if (!rmm1->isEquivalent(*rm))
    { std::cerr
	<< "  Unequal constraint matrix, row-major ordering, "
	<< si1Name << " vs. " << si2Name << std::endl ;
      return (false) ; }
    const CoinPackedMatrix *cmm1=osi1->getMatrixByCol() ;
    const CoinPackedMatrix *cm  =osi2->getMatrixByCol() ;
    if (!cmm1->isEquivalent(*cm))
    { std::cerr
	<< "  Unequal constraint matrix, column-major ordering, "
	<< si1Name << " vs. " << si2Name << std::endl ;
      return (false) ; }
  }
  // Check column types
  { areEquiv = true ;
    for (int j = 0 ; j < colCnt && areEquiv == true ; j++)
    { if (osi1->isContinuous(j) != osi2->isContinuous(j))
        areEquiv = false ;
      if (osi1->isBinary(j) != osi2->isBinary(j))
	areEquiv = false ;
      if (osi1->isIntegerNonBinary(j) != osi2->isIntegerNonBinary(j))
	areEquiv = false ;
      if (osi1->isFreeBinary(j) != osi2->isFreeBinary(j))
	areEquiv = false ;
      if (osi1->isInteger(j) != osi2->isInteger(j))
	areEquiv = false ; }
    if (areEquiv == false)
    { std::cerr
	<< "  Unequal variable type, "
	<< si1Name << " vs. " << si2Name << std::endl ;
      return (false) ; }
  }
  return (true) ;
}


} // end OsiUnitTest namespace
