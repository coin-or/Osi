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

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "OsiSolverInterface.hpp"
/*
  #include "CoinFloatEqual.hpp"
  #include "CoinWarmStartBasis.hpp"
  #include "OsiRowCut.hpp"
  #include "OsiCuts.hpp"
  #include "OsiPresolve.hpp"
*/


using namespace OsiUnitTest ;

/*
  Helper methods for the OSI simplex API test.
*/
namespace {

/*
  Test that the given vector is a unit vector with a 1 in the specified index
  position
*/
bool isUnitVector (int ndx, int len, double *vec)
{ bool verbose = true ;
  bool retval = false ; 

  CoinAbsFltEq fltEq ;

  int nzCount = 0 ;
  int oneCount = 0 ;
  int onePosn = -1 ;

  for (int j = 0 ; j < len ; j++)
  { if (!fltEq(vec[j],0.0))
    { nzCount++ ;
      if (fltEq(vec[j],1.0))
      { oneCount++ ;
	onePosn = j ; } } }
  if (nzCount == 1 && oneCount == 1 && onePosn >= 0)
  { retval = true ; }
  
  if (verbose && !retval)
  { if (nzCount > oneCount)
    { std::cout
	<< "    Vector contains " << nzCount-oneCount
	<< " elements that are neither 1.0 or 0.0." << std::endl ; }
    if (oneCount > 1)
    { std::cout
	<< "    Vector contains " << oneCount
	<< " elements that are 1.0." << std::endl ; }
    if (oneCount < 1)
    { std::cout
	<< "    Vector contains no elements that are 1.0." << std::endl ; } }

  return (retval) ; }
  

/*
  Build a constraint system in basis order. Assumes that enableFactorization
  has been called. Indirect test of getBasics, because if the basis columns are
  loaded incorrectly, checks that B inv(B) = I will surely fail.
*/
CoinPackedMatrix *buildBasisMatrix (const OsiSolverInterface *si)
{ const bool verbose = false ;
  std::string solverName ;
  si->getStrParam(OsiSolverName,solverName) ;

  CoinPackedMatrix *basisMtx = new CoinPackedMatrix() ;
  const CoinPackedMatrix *mtx = si->getMatrixByCol() ;

  int m = si->getNumRows() ;
  int n = si->getNumCols() ;
  int *basicIndices = new int[m] ;
  si->getBasics(basicIndices) ;

  for (int i = 0 ; i < m ; i++)
  { int j = basicIndices[i] ;
    if (j < n)
    { if (verbose)
      { std::cout
	  << "  Retrieving column " << j << " for basis pos'n " << i << "."
	  << std::endl ; }
      CoinShallowPackedVector col = mtx->getVector(j) ;
      basisMtx->appendCol(col) ; }
    else
    { j -= n ;
      if (verbose)
      { std::cout
	  << "  Fabricating e<" << j << "> for basis pos'n " << i << "."
	  << std::endl ; }
      CoinPackedVector ei = CoinPackedVector(1,&j,1.0) ;
      basisMtx->appendCol(ei) ; } }

  return (basisMtx) ; }


/*
  Test columns beta<k> = inv(B)e<k> of the basis inverse by calculating
  B beta<k> and checking that the result is the appropriate unit vector.
  Also tests getBasics, because we use it to build the basis matrix.

  It's assumed that the si passed as a parameter is ready for tableau
  queries: a problem has been loaded and solved to optimality and
  enableFactorization has been called.
*/
int testBInvCol (const OsiSolverInterface *si)
{ const bool verbose = false ;
  std::string solverName ;
  si->getStrParam(OsiSolverName,solverName) ;

  int m = si->getNumRows() ;

  int errCnt = 0 ;

  std::cout << "  Testing getBInvCol ... " ;

  CoinPackedMatrix *basisMtx = buildBasisMatrix(si) ;
/*
  Fetch beta<k>, calculate B beta<k>, k = 0, ..., m-1, and check that we have
  the appropriate unit vector.
*/
  double *betak = new double[m] ;
  double *ek = new double[m] ;
  for (int k = 0 ; k < m ; k++)
  { 
    CoinFillN(betak,m,COIN_DBL_MAX) ;
    CoinFillN(ek,m,COIN_DBL_MAX) ;
    try
    { si->getBInvCol(k,betak) ; }
    catch (CoinError e)
    { std::cout
	<< "  getBInvCol threw on request for beta<" << k
	<< ">." << std::endl ;
      std::cout
	<< "  CoinError: " << e.fileName() << ":" << e.methodName() << ": "
	<< e.message() << std::endl ; }
    basisMtx->times(betak,ek) ;
    if (!isUnitVector(k,m,ek))
    { errCnt++ ;
      if (verbose)
      { std::cout
	  << "  " << "B beta<" << k
	  << "> != e<" << k << ">." << std::endl ; } } }
  delete[] betak ;
  delete[] ek ;
  delete basisMtx ;
/*
  Announce the result and we're done.
*/
  if (errCnt != 0)
  { std::cout << errCnt << " errors." << std::endl ; }
  else
  { std::cout << " ok." << std::endl ; }

  return (errCnt) ; }

/*
  Test rows beta<i> = e<i>inv(B) of the basis inverse by calculating beta<i>B
  and checking that the result is the appropriate unit vector. Also tests
  getBasics, because we need it to build the basis matrix.

  It's assumed that the si passed as a parameter is ready for tableau queries:
  a problem has been loaded and solved to optimality and enableFactorization
  has been called.
*/
int testBInvRow (const OsiSolverInterface *si)
{ const bool verbose = false ;
  std::string solverName ;
  si->getStrParam(OsiSolverName,solverName) ;

  int m = si->getNumRows() ;

  int errCnt = 0 ;

  std::cout << "  Testing getBInvRow ... " ;

/*
  Should construct in row-major form for transposeTimes, but efficiency is
  not all that big an issue here.
*/
  CoinPackedMatrix *basisMtx = buildBasisMatrix(si) ;
/*
  Fetch beta<i>, calculate beta<i> B, i = 0, ..., m-1, and check that we have
  the appropriate unit vector.
*/
  double *betai = new double[m] ;
  double *ei = new double[m] ;
  for (int i = 0 ; i < m ; i++)
  { CoinFillN(betai,m,COIN_DBL_MAX) ;
    CoinFillN(ei,m,COIN_DBL_MAX) ;
    try
    { si->getBInvRow(i,betai) ; }
    catch (CoinError e)
    { std::cout
	<< "  getBInvRow threw on request for beta<" << i
	<< ">." << std::endl ;
      std::cout
	<< "  CoinError: " << e.fileName() << ":" << e.methodName() << ": "
	<< e.message() << std::endl ; }
    basisMtx->transposeTimes(betai,ei) ;
    if (!isUnitVector(i,m,ei))
    { errCnt++ ;
      if (verbose)
      { std::cout
	  << "  " << "beta<" << i
	  << ">B != e<" << i << ">." << std::endl ; } } }
  delete[] betai ;
  delete[] ei ;
  delete basisMtx ;
/*
  Announce the result and we're done.
*/
  if (errCnt != 0)
  { std::cout << errCnt << " errors." << std::endl ; }
  else
  { std::cout << " ok." << std::endl ; }

  return (errCnt) ; }


/*
  Test columns abar<j> = inv(B) a<j> by checking that B abar<j> = a<j>.

  It's assumed that the si passed as a parameter is ready for tableau
  queries: a problem has been loaded and solved to optimality and
  enableFactorization has been called.
*/
int testBInvACol (const OsiSolverInterface *si)
{ const bool verbose = false ;
  std::string solverName ;
  si->getStrParam(OsiSolverName,solverName) ;

  int n = si->getNumCols() ;
  int m = si->getNumRows() ;

  int errCnt = 0 ;

  std::cout << "  Testing getBInvACol ... " ;

  CoinPackedMatrix *basisMtx = buildBasisMatrix(si) ;
  const CoinPackedMatrix *mtx = si->getMatrixByCol() ;
/*
  Fetch abar<j>, calculate B abar<j>, k = 0, ..., n-1, and check that the
  result is a<j>.
*/
  double *abarj = new double[m] ;
  double *aj = new double[m] ;
  for (int j = 0 ; j < n ; j++)
  { CoinFillN(abarj,m,COIN_DBL_MAX) ;
    CoinFillN(aj,m,COIN_DBL_MAX) ;
    try
    { si->getBInvACol(j,abarj) ; }
    catch (CoinError e)
    { std::cout
	<< "  getBInvACol threw on request for abar<" << j
	<< ">." << std::endl ;
      std::cout
	<< "  CoinError: " << e.fileName() << ":" << e.methodName() << ": "
	<< e.message() << std::endl ; }
    basisMtx->times(abarj,aj) ;
    const CoinShallowPackedVector pv = mtx->getVector(j) ;
    if (isEquivalent(pv,m,aj) == false)
    { errCnt++ ;
      if (verbose)
      { std::cout
	  << "  " << "B abar<" << j
	  << "> != a<" << j << ">." << std::endl ; } } }
  delete[] abarj ;
  delete[] aj ;
  delete basisMtx ;
/*
  Announce the result and we're done.
*/
  if (errCnt != 0)
  { std::cout << errCnt << " errors." << std::endl ; }
  else
  { std::cout << " ok." << std::endl ; }

  return (errCnt) ; }


/*
  Test rows abar<i> = e<i>(inv(B)(A I)). This is an awkward thing to check,
  because there's no analog to the column identity B abar<j> = a<j>. Go with
  brute force: Build inv(B)A row by row with getBInvARow and check it against
  inv(B)A built column by column with getBInvACol. (Clearly, testBInvACol
  should run first.)
  
  e<i>(inv(B)I) is of course beta<i> = e<i>inv(B), so we can just check that
  beta<i>B = e<i>.
  
  It's assumed that the si passed as a parameter is ready for tableau
  queries: a problem has been loaded and solved to optimality and
  enableFactorization has been called.
*/
int testBInvARow (const OsiSolverInterface *si)
{ const bool verbose = false ;
  std::string solverName ;
  si->getStrParam(OsiSolverName,solverName) ;

  int n = si->getNumCols() ;
  int m = si->getNumRows() ;

  int errCnt = 0 ;

  std::cout << "  Testing getBInvARow ... " ;

  CoinPackedMatrix *basisMtx = buildBasisMatrix(si) ;
/*
  Construct inv(B)A by column, then change over to row-major ordering so we
  can compare with the version build from tableau rows.
  
  Interesting quirk here: Turns out p0033's tableau has no nonzeros in rows
  6 or 15. Because of that, when abarj is converted to row-major, it has only
  14 rows, and that causes isEquivalent2 to fail. So force the full size.
*/
  CoinPackedMatrix abarjMtx ;
  double *abarj = new double[m] ;
  for (int j = 0 ; j < n ; j++)
  { si->getBInvACol(j,abarj) ;
    CoinPackedVector pkv ;
    pkv.setFullNonZero(m,abarj) ;
    abarjMtx.appendCol(pkv) ; }
  delete[] abarj ;
  abarjMtx.reverseOrdering() ;
  abarjMtx.setDimensions(m,n) ;
  if (verbose)
  { std::cout
      << "  Col-major tableau is " << abarjMtx.getNumRows() << " x "
      << abarjMtx.getNumCols() << " with " << abarjMtx.getNumElements()
      << " elements." << std::endl ; }
/*
  Construct inv(B)A by row. Check the vectors returned for inv(B)I = beta<i>
  as we go.
*/
  CoinPackedMatrix abariMtx ;
  abariMtx.reverseOrdering() ;
  double *abari = new double[n] ;
  double *betai = new double[m] ;
  double *ei = new double[m] ;
  for (int i = 0 ; i < m ; i++)
  { CoinFillN(abari,n,COIN_DBL_MAX) ;
    CoinFillN(betai,m,COIN_DBL_MAX) ;
    try
    { si->getBInvARow(i,abari,betai) ; }
    catch (CoinError e)
    { std::cout
	<< "  getBInvARow threw on request for abar<" << i
	<< ">." << std::endl ;
      std::cout
	<< "  CoinError: " << e.fileName() << ":" << e.methodName() << ": "
	<< e.message() << std::endl ; }
    CoinPackedVector pkv ;
    pkv.setFullNonZero(n,abari) ;
    if (verbose)
    { std::cout << "  Adding" ;
      const int *indices = pkv.getIndices() ;
      for (int v = 0 ; v < pkv.getNumElements() ; v++)
      { std::cout << " (" << i << "," << indices[v] << ")" ; }
      std::cout << std::endl ;
      if (!isEquivalent(pkv,n,abari))
        std::cout << "  !! packed abari != full abari !!" << std::endl ; }
    abariMtx.appendRow(pkv) ;
    basisMtx->transposeTimes(betai,ei) ;
    if (!isUnitVector(i,m,ei))
    { errCnt++ ;
      if (verbose)
      { std::cout
	  << "  " << "beta<" << i
	  << ">B != e<" << i << ">." << std::endl ; } } }
  abariMtx.setDimensions(m,n) ;
  if (verbose)
  { std::cout
      << "  Row-major tableau is " << abariMtx.getNumRows() << " x "
      << abariMtx.getNumCols() << " with " << abariMtx.getNumElements()
      << " elements." << std::endl ; }
  delete[] abari ;
  delete[] betai ;
  delete[] ei ;
  delete basisMtx ;
/*
  Check that the two matrices are equivalent. isEquivalent2 will report
  differences, but we won't get a good count. But then, one error is all we
  need to report.
*/
  if (!abariMtx.isEquivalent2(abarjMtx))
  { std::cout
      << "  Tableau built by rows does not match tableau built by columns."
      << std::endl ;
    errCnt++ ; }
/*
  Announce the result and we're done.
*/
  if (errCnt != 0)
  { std::cout << errCnt << " errors." << std::endl ; }
  else
  { std::cout << " ok." << std::endl ; }

  return (errCnt) ; }

/*
  Test the mode 2 portion of the simplex API.
  Solve an lp by hand
*/
int testSimplexMode2 (const OsiSolverInterface *emptySi, std::string mpsDir)
{ OsiSolverInterface * si = emptySi->clone();
  std::string solverName;
  si->getStrParam(OsiSolverName,solverName);
  std::string fn = mpsDir+"p0033";
  si->readMps(fn.c_str(),"mps");
  si->setObjSense(-1.0);
  si->initialSolve();
  si->setObjSense(1.0);
  // enable special mode
  si->enableSimplexInterface(true);
  // we happen to know that variables are 0-1 and rows are L
  int numberIterations=0;
  int numberColumns = si->getNumCols();
  int numberRows = si->getNumRows();
  double * fakeCost = new double[numberColumns];
  double * duals = new double [numberRows];
  double * djs = new double [numberColumns];
  const double * solution = si->getColSolution();
  memcpy(fakeCost,si->getObjCoefficients(),numberColumns*sizeof(double));
  while (1) {
    const double * dj;
    const double * dual;
    if ((numberIterations&1)==0) {
      // use given ones
      dj = si->getReducedCost();
      dual = si->getRowPrice();
    } else {
      // create
      dj = djs;
      dual = duals;
      si->getReducedGradient(djs,duals,fakeCost);
    }
    int i;
    int colIn=9999;
    int direction=1;
    double best=1.0e-6;
    // find most negative reduced cost
    // Should check basic - but should be okay on this problem
    for (i=0;i<numberRows;i++) {
      double value=dual[i];
      if (value>best) {
	direction=-1;
	best=value;
	colIn=-i-1;
      }
    }
    for (i=0;i<numberColumns;i++) {
      double value=dj[i];
      if (value<-best&&solution[i]<1.0e-6) {
	direction=1;
	best=-value;
	colIn=i;
      } else if (value>best&&solution[i]>1.0-1.0e-6) {
	direction=-1;
	best=value;
	colIn=i;
      }
    }
    if (colIn==9999)
      break; // should be optimal
    int colOut;
    int outStatus;
    double theta;
    assert(!si->primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL));
    printf("out %d, direction %d theta %g\n",
	   colOut,outStatus,theta);
    numberIterations++;
  }
  delete [] fakeCost;
  delete [] duals;
  delete [] djs;
  // exit special mode
  si->disableSimplexInterface();
  si->resolve();
  assert (!si->getIterationCount());
  si->setObjSense(-1.0);
  si->initialSolve();
  std::cout<<solverName<<" passed OsiSimplexInterface test"<<std::endl;
  delete si;

  return (0) ; }

/*
  Test Simplex API mode 1 (tableau access) methods.
*/
int testSimplexMode1 (const OsiSolverInterface *emptySi, std::string mpsDir)
{ bool verbose = false ;

  OsiSolverInterface * si = emptySi->clone();
  std::string solverName;
  si->getStrParam(OsiSolverName,solverName);
  si->setHintParam(OsiDoReducePrint,false,OsiHintDo) ;

  int errCnt = 0 ;
/*
  Read and solve p0033 (maximise) to provide a baseline.
*/
  std::string fn = mpsDir+"p0033";
  si->readMps(fn.c_str(),"mps");
  si->setObjSense(-1.0);
/*
  Check that there's no optimal basis available at this point.
*/
  bool testVal = si->basisIsAvailable() ;
  if (testVal)
  { failureMessage(*si,"Optimal basis available before initial solve.") ;
    errCnt++ ; }
  else
  if (verbose)
  { std::cout
      << "  " << solverName << " shows no optimal basis before initial solve."
      << std::endl ; }
  si->initialSolve();
  if (!si->isProvenOptimal())
  { failureMessage(*si,"Failed to solve p0033 to optimality.") ;
    return (1) ; }
  if (verbose)
  { std::cout
      << "  " << solverName << " solved p0033 z = " << si->getObjValue()
      << "." << std::endl ; }
  testVal = si->basisIsAvailable() ;
  if (!testVal)
  { failureMessage(*si,"No optimal basis available after initial solve.") ;
    errCnt++ ; }
  else
  if (verbose)
  { std::cout
      << "  " << solverName << " shows optimal basis after initial solve."
      << std::endl ; }
/*
  Enable simplex mode 1.
*/
  si->enableFactorization() ;
/*
  Test the various methods.
*/
  errCnt += testBInvCol(si) ;
  errCnt += testBInvRow(si) ;
  errCnt += testBInvACol(si) ;
  errCnt += testBInvARow(si) ;
/*
  Disable simplex mode 1.
*/
  si->disableFactorization() ;
/*
  Trash this solver and we're finished.
*/
  delete si ;

  return (errCnt) ;
}

} // end file-local namespace

/*
  Test a solver's implementation of the OSI simplex API.
*/

int testSimplexAPI (const OsiSolverInterface *emptySi,
		    const std::string &mpsDir)
{ OsiSolverInterface *si = emptySi->clone() ;
  std::string solverName;
  si->getStrParam(OsiSolverName,solverName);
/*
  Do the tests only if the solver implements the simplex API.
*/
  if (si->canDoSimplexInterface() == 0)
  { std::cout
      << solverName << " has no OsiSimplex API." << std::endl ;
    return (0) ; }
/*
  Test the mode 1 (tableau access) portion of the API.
*/
  int totalErrs = 0 ;
  if (si->canDoSimplexInterface() >= 1)
  { std::cout
      << "Testing Simplex API mode 1 for " << solverName << " ... "
      << std::endl ;
    int errCnt = testSimplexMode1(emptySi,mpsDir) ;
    totalErrs += errCnt ;
    if (errCnt > 0)
    { std::cout
        << "  Simplex API mode 1 tests incurred " << errCnt << " errors."
	<< std::endl ; }
    else
    { std::cout << "  Simplex API mode 1 ok." << std::endl ; }
  }
/*
  Test the mode 2 (pivot-by-pivot control) portion of the API.
*/
  if (si->canDoSimplexInterface() >= 2)
  { std::cout
      << "Testing Simplex API mode 2 for " << solverName << " ... "
      << std::endl ;
    int errCnt = testSimplexMode2(emptySi,mpsDir) ;
    totalErrs += errCnt ;
    if (errCnt > 0)
    { std::cout
        << "  Simplex API mode 1 tests incurred " << errCnt << " errors."
	<< std::endl ; }
    else
    { std::cout << "  Simplex API mode 2 ok." << std::endl ; }
  }

  return (totalErrs) ; }
