// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "CoinPackedMatrix.hpp"
#include "OsiOslSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

//#############################################################################
// A helper function that converts a packed matrix into a columnordered
// gap-free matrix.
//#############################################################################

static const CoinPackedMatrix *
toColumnOrderedGapFree(const CoinPackedMatrix& matrix)
{
   CoinPackedMatrix * m = 0;
   if (matrix.isColOrdered()) {
      const CoinBigIndex * start = matrix.getVectorStarts();
      const int * length = matrix.getVectorLengths();
      int i;
      for (i = matrix.getNumCols() - 1; i >= 0; --i)
	 if (start[i+1] - start[i] != length[i])
	    break;
      if (i >= 0) {
	 // got to get rid of the gaps
        m = new CoinPackedMatrix(true,matrix.getNumRows(),matrix.getNumCols(),
                                 matrix.getNumElements(),matrix.getElements(),
                                 matrix.getIndices(),matrix.getVectorStarts(),
                                 matrix.getVectorLengths(),0.0,0.0);
      }
   } else {
      // must create a column ordered copy without gaps
      m = new CoinPackedMatrix();
      m->setExtraGap(0.0);
      m->setExtraMajor(0.0);
      m->reverseOrderedCopyOf(matrix);
   }
   return m ? m : &matrix;
}

//#############################################################################
// Solve methods
//#############################################################################
void OsiOslSolverInterface::initialSolve()
{
  EKKModel* model = getMutableModelPtr();
  int rtcod;
  // Switch off printing if asked to
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = getHintParam(OsiDoReducePrint,takeHint,strength);
  assert(gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    if (!messageHandler()->logLevel())
      ekk_messagesPrintOff(model,1,5999);
    else if (messageHandler()->logLevel()==1)
      ekk_messagePrintOff(model,317);
  }
#if 0
  ekk_crash(model,1); 
  rtcod=ekk_primalSimplex(model,1);
#elif 0
  ekk_setIiternum(model, 0);
  rtcod=ekk_simplex(model, 2); // *FIXME* : why not 0 (eitherSimplex)
#else
  ekk_setIiternum(model, 0);
  rtcod=ekk_simplex(model, 0);
#endif
  // If abandoned use B&B status which says abandoned
  if (rtcod>200)
    ekk_setIprobstat(model,6);
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::resolve()
{
  EKKModel* model = getMutableModelPtr();
  int rtcod;

  ekk_mergeBlocks(model, 1);
  ekk_setIiternum(model, 0);
  // Switch off printing if asked to
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = getHintParam(OsiDoReducePrint,takeHint,strength);
  assert(gotHint);
  if (strength!=OsiHintIgnore&&takeHint) {
    if (!messageHandler()->logLevel())
      ekk_messagesPrintOff(model,1,5999);
    else if (messageHandler()->logLevel()==1)
      ekk_messagePrintOff(model,317);
  }
#if 0
  rtcod=ekk_dualSimplex(model); // *FIXME* : why not 0 (eitherSimplex)
#else
  rtcod=ekk_simplex(model, 256 + 32); // no presolve and no scaling
#endif
  // If abandoned use B&B status which says abandoned
  if (rtcod>200)
    ekk_setIprobstat(model,6);
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::branchAndBound()
{
  EKKModel* model = getMutableModelPtr();

  ekk_primalSimplex(model,1);
#if 0
  EKKCuts cuts;
  cuts.numberCuts = 0;
  cuts.maxCuts = 0;
  cuts.cut = NULL;
  
  EKKIntegerPresolve info;
  info.number01 = 0;
  info.sequence = NULL;
  info.numberClique = 0;
  info.cliqueType = NULL;
  info.cliqueMember = NULL;
  info.cliqueStart = NULL;
  info.implicationStart = NULL;
  info.implication = NULL;
  info.numberChains = 0;
  info.chainInformation = NULL;
  
  ekk_integerPresolve(model, &info, &cuts, 0);
  // OSL computes a threshold to stop adding cuts if change in objective
  // is less than this, but it may be too small.
  ekk_setRthreshold(model, CoinMax(1.0e-3, ekk_getRthreshold(model)));
  ekk_branchAndCut(model, NULL, NULL, &info, &cuts, 0, 0);
  ekk_deleteCuts(&cuts);
  ekk_deleteIntegerPresolve(&info);
#else
  //strong branching
  //  const int old_strategy = ekk_getIstrategy(model);
  //  ekk_setIstrategy(model, old_strategy|1024);
  ekk_branchAndBound(model, NULL, NULL);
  //  ekk_setIstrategy(model, old_strategy);
#endif
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiOslSolverInterface::setIntParam(OsiIntParam key, int value)
{
  switch (key) {
  case OsiMaxNumIteration:
    if (value < 0)
      return false;
    // should use getModelPtr(), but setting max iterations does not
    // invalidate anything that's cached.
    ekk_setImaxiter(getMutableModelPtr(), value);
    break;
  case OsiMaxNumIterationHotStart:
    if (value < 0)
      return false;
    OsiSolverInterface::setIntParam(key, value);
    break;
  case OsiLastIntParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiOslSolverInterface::setDblParam(OsiDblParam key, double value)
{
  int retval;
  // *THINK* : how to set opt limits for OSL?
  switch (key) {
  case OsiDualObjectiveLimit:
    retval = ekk_setRbbcutoff(getMutableModelPtr(), value);
    // OsiSolverInterface::setDblParam(key, value);
    return retval == 0;

  case OsiPrimalObjectiveLimit:
    return OsiSolverInterface::setDblParam(key, value);

  case OsiDualTolerance: // return false if not succeeded (out of range)
    retval = ekk_setRtoldinf(getMutableModelPtr(), value);
    return retval == 0;

  case OsiPrimalTolerance: // return false if not succeeded (out of range)
    retval = ekk_setRtolpinf(getMutableModelPtr(), value);
    return retval == 0;

  case OsiObjOffset: 
    retval = ekk_setRobjectiveOffset(getMutableModelPtr(), value);
    return retval == 0;

  case OsiLastDblParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiOslSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  int retval;
  switch (key) {
  case OsiProbName:
    retval = ekk_setCname(getMutableModelPtr(), value.c_str());
    return retval == 0;

  case OsiSolverName:
    return false;
  case OsiLastStrParam:
    return false;
  }
  return false;
}


//-----------------------------------------------------------------------------

bool
OsiOslSolverInterface::getIntParam(OsiIntParam key, int& value) const 
{
  switch (key) {
  case OsiMaxNumIteration:
    value = ekk_getImaxiter(getMutableModelPtr());
    break;
  case OsiMaxNumIterationHotStart:
    OsiSolverInterface::getIntParam(key, value);
    break;
  case OsiLastIntParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiOslSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  // *THINK* : how to set opt limits for OSL?
  switch (key) {
  case OsiDualObjectiveLimit:
    value = ekk_getRbbcutoff(getMutableModelPtr());
    //OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiPrimalObjectiveLimit:
    OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiDualTolerance:
    value = ekk_getRtoldinf(getMutableModelPtr());
    break;
  case OsiPrimalTolerance:
    value = ekk_getRtolpinf(getMutableModelPtr());
    break;
  case OsiObjOffset:
    value = ekk_getRobjectiveOffset(getMutableModelPtr());
    break;
  case OsiLastDblParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiOslSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  switch (key) {
  case OsiProbName:
    value = ekk_getCname(getMutableModelPtr());
    break;
  case OsiSolverName:
    value = "osl";
    break;
  case OsiLastStrParam:
    return false;
  }
  return true;
}


//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiOslSolverInterface::isAbandoned() const
{
  EKKModel* model = getMutableModelPtr();
  return (ekk_getIprobstat(model)==6);
}

bool OsiOslSolverInterface::isProvenOptimal() const
{
  EKKModel* model = getMutableModelPtr();
  const int stat = ekk_getIprobstat(model);
  return (stat == 0);
}

bool OsiOslSolverInterface::isProvenPrimalInfeasible() const
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  const int stat = ekk_getIprobstat(model);
  if (stat != 1)
     return false;
  if (ekk_lastAlgorithm(model) == 2) { // dual
     const int stat2 = ekk_getIprobstat2(model);
     if (stat2 == 11) // over dual limit
	return false;
  }
  return true;
}

bool OsiOslSolverInterface::isProvenDualInfeasible() const
{
  EKKModel* model = getMutableModelPtr();
  const int stat = ekk_getIprobstat(model);
  return stat == 2;
}

bool OsiOslSolverInterface::isPrimalObjectiveLimitReached() const
{
  // *TEST*
  double limit = 0.0;
  getDblParam(OsiPrimalObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  EKKModel* model = getMutableModelPtr();
  if (ekk_getIprobstat(model)==6)
    return false;
  const int lastalgo = ekk_lastAlgorithm(model);
  const double obj = ekk_getRobjvalue(model);
  const double maxmin = ekk_getRmaxmin(model);

  switch (lastalgo) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
   case 2: // dual simplex
     if (ekk_getIprobstat(model) == 0) // optimal
	return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
     return false;
   case 1: // primal simplex
     return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiOslSolverInterface::isDualObjectiveLimitReached() const
{
  // *TEST*
  double limit = 0.0;
  getDblParam(OsiDualObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  EKKModel* model = getMutableModelPtr();
  if (ekk_getIprobstat(model)==6)
    return false;
  const int lastalgo = ekk_lastAlgorithm(model);
  const double obj = ekk_getRobjvalue(model);
  const double maxmin = ekk_getRmaxmin(model);

  switch (lastalgo) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
   case 1: // primal simplex
     if (ekk_getIprobstat(model) == 0) // optimal
	return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
     return false;
   case 2: // dual simplex
     if (ekk_getIprobstat(model) != 0 && ekk_getIprobstat2(model) == 11)
	// over dual limit
	return true;
     return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiOslSolverInterface::isIterationLimitReached() const
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  if (ekk_getIprobstat(model)==6)
    return false;
  const int stat = ekk_getIprobstat(model);
  return (stat == 3);
}

//#############################################################################
// WarmStart related methods
//#############################################################################

CoinWarmStart* OsiOslSolverInterface::getWarmStart() const
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  const int numcols = getNumCols();
  const int numrows = getNumRows();

  CoinWarmStartBasis* ws = new CoinWarmStartBasis;
  ws->setSize(numcols, numrows);

  int i;

  for (i = 0; i < numrows; ++i) {
    switch (ekk_rowStatus(model, i)) {
    case -2:
    case -1:
      ws->setArtifStatus(i, CoinWarmStartBasis::atLowerBound);
      break;
    case  0:
      ws->setArtifStatus(i, CoinWarmStartBasis::basic);
      break;
    case  1:
    case  2:
      ws->setArtifStatus(i, CoinWarmStartBasis::atUpperBound);
      break;
    }
  }

  for (i = 0; i < numcols; ++i) {
    switch (ekk_columnStatus(model, i)) {
    case -2:
    case -1:
      ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
      break;
    case  0:
      ws->setStructStatus(i, CoinWarmStartBasis::basic);
      break;
    case  1: 
    case  2:
      ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
      break;
    }
  }
  
  return ws;
}

//-----------------------------------------------------------------------------

bool OsiOslSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  // *TEST*
  const CoinWarmStartBasis* ws =
    dynamic_cast<const CoinWarmStartBasis*>(warmstart);

  if (! ws)
    return false;

  const int numcols = ws->getNumStructural();
  const int numrows = ws->getNumArtificial();

  EKKModel* model = getMutableModelPtr();

  const int maxdim = CoinMax(numcols, numrows);
  int * atLower = new int[maxdim];
  int * atUpper = new int[maxdim];
  int * basic = new int[maxdim];
  int numLower = 0;
  int numUpper = 0;
  int numBasic = 0;
  int i;

  for (i = 0; i < numrows; ++i) {
    switch (ws->getArtifStatus(i)) {
    case CoinWarmStartBasis::atUpperBound:
      atUpper[numUpper++] = i;
      break;
    case CoinWarmStartBasis::atLowerBound:
      atLower[numLower++] = i;
      break;
    case CoinWarmStartBasis::basic:
      basic[numBasic++] = i;
      break;
    case CoinWarmStartBasis::isFree: // superbasic
      basic[numBasic++] = i;
      break;
    }
  }
  if (ekk_setRowsNonBasicAtLower(model, numLower, atLower) != 0)
    goto setWarmStart_Error;
  if (ekk_setRowsNonBasicAtUpper(model, numUpper, atUpper) != 0)
    goto setWarmStart_Error;
  if (ekk_markRowsAsBasic(model, numBasic, basic) != 0)
    goto setWarmStart_Error;

  numLower = 0;
  numUpper = 0;
  numBasic = 0;
  for (i = 0; i < numcols; ++i) {
    switch (ws->getStructStatus(i)) {
    case CoinWarmStartBasis::atUpperBound:
      atUpper[numUpper++] = i;
      break;
    case CoinWarmStartBasis::atLowerBound:
      atLower[numLower++] = i;
      break;
    case CoinWarmStartBasis::basic:
      basic[numBasic++] = i;
      break;
    case CoinWarmStartBasis::isFree: // superbasic
      basic[numBasic++] = i;
      break;
    }
  }
  if (ekk_setColumnsNonBasicAtLower(model, numLower, atLower) != 0)
    goto setWarmStart_Error;
  if (ekk_setColumnsNonBasicAtUpper(model, numUpper, atUpper) != 0)
    goto setWarmStart_Error;
  if (ekk_markColumnsAsBasic(model, numBasic, basic) != 0)
    goto setWarmStart_Error;

  delete[] basic;
  delete[] atUpper;
  delete[] atLower;
  return true;

  setWarmStart_Error:
  delete[] basic;
  delete[] atUpper;
  delete[] atLower;
  return false;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiOslSolverInterface::markHotStart()
{
  // *TEST*
  //EKKModel* model = getMutableModelPtr();
  delete ws_;
  ws_ = dynamic_cast<CoinWarmStartBasis*>(getWarmStart());
  //  ekk_startFastDualSimplex(model, itlim);
}

void OsiOslSolverInterface::solveFromHotStart()
{
  EKKModel* model = getMutableModelPtr();
  // *TEST*
  ekk_setIiternum(model, 0);
  setWarmStart(ws_);
  //  ekk_fastDualSimplex(model);
  itlimOrig_ = ekk_getImaxiter(model);
  int itlim;
  OsiSolverInterface::getIntParam(OsiMaxNumIterationHotStart, itlim);
  ekk_setImaxiter(model, itlim);
  resolve();
  ekk_setImaxiter(model, itlimOrig_);
}

void OsiOslSolverInterface::unmarkHotStart()
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  //  ekk_stopFastDualSimplex(model);
  ekk_setImaxiter(model, itlimOrig_);
  delete ws_;
  ws_ = NULL;
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows, columns and elements
//------------------------------------------------------------------
int OsiOslSolverInterface::getNumCols() const
{
  return ekk_getInumcols(getMutableModelPtr());
}
int OsiOslSolverInterface::getNumRows() const
{
  return ekk_getInumrows(getMutableModelPtr());
}
int OsiOslSolverInterface::getNumElements() const
{
  return ekk_getInumels(getMutableModelPtr());
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  
const double * OsiOslSolverInterface::getColLower() const
{
  return ekk_collower(getMutableModelPtr());
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getColUpper() const
{
  return ekk_colupper(getMutableModelPtr());
}
//------------------------------------------------------------------
const char * OsiOslSolverInterface::getRowSense() const
{
  extractSenseRhsRange();
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getRightHandSide() const
{
  extractSenseRhsRange();
  return rhs_;
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getRowRange() const
{
  extractSenseRhsRange();
  return rowrange_;
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getRowLower() const
{
  return ekk_rowlower(getMutableModelPtr());
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getRowUpper() const
{
  return ekk_rowupper(getMutableModelPtr());
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getObjCoefficients() const
{
  return ekk_objective(getMutableModelPtr());
}
//------------------------------------------------------------------
double OsiOslSolverInterface::getObjSense() const
{
  return ekk_getRmaxmin(getMutableModelPtr());
}

//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------
bool OsiOslSolverInterface::isContinuous(int colNumber) const
{
  const char * intType = ekk_integerType(getMutableModelPtr());
  if ( intType==NULL ) return true;
  if ( intType[colNumber]==0 ) return true;
  return false;
}
//------------------------------------------------------------------
#if 0
bool OsiOslSolverInterface::isInteger(
      int columnNumber ) const
{
  return !(isContinuous(columnNumber));
}
//------------------------------------------------------------------
bool OsiOslSolverInterface::isBinary(
  int columnNumber ) const
{
  if ( isContinuous(columnNumber) ) return false;  
  EKKModel * m=getMutableModelPtr();
  const double * upper = ekk_colupper(m);
  const double * lower = ekk_collower(m);
  if (
    (upper[columnNumber]== 1 || upper[columnNumber]== 0) && 
    (lower[columnNumber]== 0 || lower[columnNumber]==1)
    ) return true;
  else return false;
}
//------------------------------------------------------------------
bool OsiOslSolverInterface::isIntegerNonBinary(
  int columnNumber ) const
{
  if ( isInteger(columnNumber) && !isBinary(columnNumber) )
    return true; 
  else return false;
}
//------------------------------------------------------------------
bool OsiOslSolverInterface::isFreeBinary(
  int columnNumber ) const
{
  if ( isContinuous(columnNumber) ) return false;
  EKKModel * m=getMutableModelPtr();
  const double * upper = ekk_colupper(m);
  const double * lower = ekk_collower(m);
  if (
    (upper[columnNumber]== 1) &&
    (lower[columnNumber]== 0)
    ) return true;
  else return false;
}
#endif

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------
const CoinPackedMatrix * OsiOslSolverInterface::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL ) {
    EKKMatrixCopy rowCopy;
    ekk_createRowCopy(getMutableModelPtr(), &rowCopy);
    matrixByRow_ = new CoinPackedMatrix();
    matrixByRow_->copyOf(false /* not column ordered */,
			 getNumCols(), getNumRows(), getNumElements(),
			 rowCopy.element, rowCopy.index, (CoinBigIndex *) rowCopy.start,
			 NULL /* compute lengths */);
    ekk_free(rowCopy.element);
    ekk_free(rowCopy.index);
    ekk_free(rowCopy.start);
  }
  return matrixByRow_;
}
//------------------------------------------------------------------
const CoinPackedMatrix * OsiOslSolverInterface::getMatrixByCol() const
{
  if ( matrixByColumn_ == NULL ) {
    EKKMatrixCopy colCopy;
    ekk_createColumnCopy(getMutableModelPtr(), &colCopy);
    matrixByColumn_ = new CoinPackedMatrix();
    matrixByColumn_->copyOf(true /* column ordered */,
			    getNumRows(), getNumCols(), getNumElements(),
			    colCopy.element, colCopy.index, (CoinBigIndex *) colCopy.start,
			    NULL /* compute lengths */);
    ekk_free(colCopy.element);
    ekk_free(colCopy.index);
    ekk_free(colCopy.start);
  }
  return matrixByColumn_;
}
//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double OsiOslSolverInterface::getInfinity() const
{
  return OSL_INFINITY;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

const double * OsiOslSolverInterface::getColSolution() const
{
  return ekk_colsol(getMutableModelPtr());
}
//------------------------------------------------------------------

const double * OsiOslSolverInterface::getRowPrice() const
{
  return ekk_rowduals(getMutableModelPtr());
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getReducedCost() const
{
  // *TEST*
  return ekk_colrcosts(getMutableModelPtr());
}
//------------------------------------------------------------------
const double * OsiOslSolverInterface::getRowActivity() const
{
  // *TEST*
  return ekk_rowacts(getMutableModelPtr());
}
//------------------------------------------------------------------
double OsiOslSolverInterface::getObjValue() const
{
#if 0
  // This does not pass unitTest if getObjValue is called
  // before solving.
  return ekk_getRobjvalue(getMutableModelPtr());
#else
  return OsiSolverInterface::getObjValue();
#endif
}
//------------------------------------------------------------------
int OsiOslSolverInterface::getIterationCount() const
{
  return ekk_getIiternum(getMutableModelPtr());
}
//------------------------------------------------------------------
std::vector<double*> OsiOslSolverInterface::getDualRays(int maxNumRays) const
{
  const int m = getNumRows();
  double* ray = new double[m];
  const double* negray = ekk_rowaux(getMutableModelPtr());
  for (int i = 0; i < m; ++i) {
     ray[i] = -negray[i];
  }
  return std::vector<double*>(1, ray);
}
//------------------------------------------------------------------
std::vector<double*> OsiOslSolverInterface::getPrimalRays(int maxNumRays) const
{
  const int n = getNumCols();
  double* ray = new double[n];
  CoinDisjointCopyN(ekk_colaux(getMutableModelPtr()), n, ray);
  return std::vector<double*>(1, ray);
}
//------------------------------------------------------------------
#if 0
OsiVectorInt
OsiOslSolverInterface::getFractionalIndices(const double etol) const
{
  OsiVectorInt retVal;
  EKKModel * m = getMutableModelPtr();
  int numInts = ekk_getInumints(m);
  int * intInx = ekk_listOfIntegers(m);
  const double * colSolVec = ekk_colsol(m); 
  OsiAbsFltEq eq(etol);
  for ( int i=0; i<numInts; i++ ) {
    const double colSolElem = colSolVec[intInx[i]];
    const double distanceFromInteger = colSolElem - floor(colSolElem + 0.5);
    if ( !eq( distanceFromInteger, 0.0 ) ) {
      retVal.push_back(intInx[i]);
    }
  }
  ekk_free( intInx );
  return retVal;
}
#endif

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################
void
OsiOslSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  ekk_copyObjective(model, &elementValue, elementIndex, elementIndex+1 );
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setColLower( int elementIndex, double elementValue )
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  const double newVal =
    forceIntoRange(elementValue, -OSL_INFINITY, OSL_INFINITY);
  ekk_copyCollower(model, &newVal, elementIndex, elementIndex+1 );
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setColUpper( int elementIndex, double elementValue )
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  const double newVal =
    forceIntoRange(elementValue, -OSL_INFINITY, OSL_INFINITY);
  ekk_copyColupper(model, &newVal, elementIndex, elementIndex+1 );
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  // *TEST*
  const int numcols = getNumCols();
  if (indexLast - indexFirst < numcols/8) {
    // It's probably more effective to invoke setColBounds() iteratively
    while (indexFirst != indexLast) {
      setColBounds(*indexFirst, boundList[0], boundList[1]);
      ++indexFirst;
      boundList += 2;
    }
  } else {
    // probably faster to set everything at once
    EKKModel* model = getMutableModelPtr();
    double * lower = ekk_getCollower(model);
    double * upper = ekk_getColupper(model);
    while (indexFirst != indexLast) {
      const int iCol=*indexFirst++;
      lower[iCol]= forceIntoRange(*boundList++, -OSL_INFINITY, OSL_INFINITY);
      upper[iCol]= forceIntoRange(*boundList++, -OSL_INFINITY, OSL_INFINITY);
    }
    ekk_setCollower(model, lower);
    ekk_setColupper(model, upper);
    ekk_free(lower);
    ekk_free(upper);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setRowLower( int i, double elementValue )
{
  // *TEST*
  EKKModel* model = getMutableModelPtr();
  const double newVal =
    forceIntoRange(elementValue, -OSL_INFINITY, OSL_INFINITY);
  ekk_copyRowlower(model, &newVal, i, i+1 );
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(newVal, getRowUpper()[i],
			rowsense_[i], rhs_[i], rowrange_[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setRowUpper( int i, double elementValue )
{
  // *TEST*
  const double newVal =
    forceIntoRange(elementValue, -OSL_INFINITY, OSL_INFINITY);
  ekk_copyRowupper(getMutableModelPtr(), &newVal, i, i+1 );
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(getRowLower()[i], newVal,
			rowsense_[i], rhs_[i], rowrange_[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setRowBounds( int i,
				     double lower, double upper )
{
  // *TEST*
  const double newLower = forceIntoRange(lower, -OSL_INFINITY, OSL_INFINITY);
  const double newUpper = forceIntoRange(upper, -OSL_INFINITY, OSL_INFINITY);
  EKKModel* model = getMutableModelPtr();
  ekk_copyRowlower(model, &newLower, i, i+1 );
  ekk_copyRowupper(model, &newUpper, i, i+1 );
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    convertBoundToSense(newLower, newUpper,
			rowsense_[i], rhs_[i], rowrange_[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  // *TEST*
  double lower, upper;
  convertSenseToBound(sense, rightHandSide, range, lower, upper);
  setRowBounds(i, lower, upper);
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  // *TEST*
  const int numrows = getNumRows();
  if (indexLast - indexFirst < numrows / 16) {
    // It's probably more effective to invoke setRowBounds() iteratively
    while (indexFirst != indexLast) {
      setRowBounds(*indexFirst, boundList[0], boundList[1]);
      ++indexFirst;
      boundList += 2;
    }
  } else {
    // probably faster to set everything at once
    EKKModel* model = getMutableModelPtr();
    double * lower = ekk_getRowlower(model);
    double * upper = ekk_getRowupper(model);
    const int len = indexLast - indexFirst;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      lower[iRow]= forceIntoRange(*boundList++, -OSL_INFINITY, OSL_INFINITY);
      upper[iRow]= forceIntoRange(*boundList++, -OSL_INFINITY, OSL_INFINITY);
    }
    if (rowsense_ != NULL) {
      assert ((rhs_ != NULL) && (rowrange_ != NULL));
      indexFirst -= len;
      while (indexFirst != indexLast) {
	const int iRow=*indexFirst++;
	convertBoundToSense(lower[iRow], upper[iRow],
			    rowsense_[iRow], rhs_[iRow], rowrange_[iRow]);
      }
    }
    ekk_setRowlower(model, lower);
    ekk_setRowupper(model, upper);
    ekk_free(lower);
    ekk_free(upper);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  // *TEST*
  const int numrows = getNumRows();
  if (indexLast - indexFirst < numrows / 16) {
    // It's probably more effective to invoke setRowType() iteratively
    while (indexFirst != indexLast) {
      if (rangeList){
	setRowType(*indexFirst++, *senseList++, *rhsList++, *rangeList++);
      } else {
	setRowType(*indexFirst++, *senseList++, *rhsList++, 0);
      }	
    }
  } else {
    // probably faster to set everything at once
    EKKModel* model = getMutableModelPtr();
    double * lower = ekk_getRowlower(model);
    double * upper = ekk_getRowupper(model);
    const int len = indexLast - indexFirst;
    while (indexFirst != indexLast) {
      const int iRow= *indexFirst++;
      if (rangeList){
	convertSenseToBound(*senseList++, *rhsList++, *rangeList++,
			    lower[iRow], upper[iRow]);
      } else {
	convertSenseToBound(*senseList++, *rhsList++, 0,
			    lower[iRow], upper[iRow]);
      }
    }
    if (rowsense_ != NULL) {
      assert ((rhs_ != NULL) && (rowrange_ != NULL));
      indexFirst -= len;
      senseList -= len;
      rhsList -= len;
      if (rangeList)
	rangeList -= len;
      while (indexFirst != indexLast) {
	const int iRow=*indexFirst++;
	rowsense_[iRow] = *senseList++;
	rhs_[iRow] = *rhsList++;
	if (rangeList)
	  rowrange_[iRow] = *rangeList++;
      }
    }
    ekk_setRowlower(model, lower);
    ekk_setRowupper(model, upper);
    ekk_free(lower);
    ekk_free(upper);
  }
}
//#############################################################################
void
OsiOslSolverInterface::setContinuous(int index)
{
  // *FIXME*: This is really painful, but this is the only way...
  EKKModel* model = getMutableModelPtr();
  int intnum = ekk_getInumints(model);
  if (intnum > 0) {
    int * intlist = ekk_listOfIntegers(model);
    int * inspoint = std::lower_bound(intlist, intlist + intnum, index);
    if (*inspoint == index) {
      ekk_deleteIntegerInformation(model);
      if (intnum > 1) {
	// *inspoint = intlist[intnum - 1];
	CoinCopy(inspoint + 1, intlist + intnum, inspoint);
	ekk_addIntegerSet(model, 100, intnum - 1, intlist, NULL, NULL);
      }
    }
    ekk_free(intlist);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setInteger(int index)
{
  ekk_markAsInteger(getMutableModelPtr(), index);
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setContinuous(const int* indices, int len)
{
  // *FIXME*: This is really painful, but this is the only way...
  EKKModel* model = getMutableModelPtr();
  int intnum = ekk_getInumints(model);
  if (intnum > 0 && len > 0) {
    int * intlist = ekk_listOfIntegers(model);
    int * sorted = new int[len];
    CoinDisjointCopyN(indices, len, sorted);
    std::sort(sorted, sorted+len);
    int * diff = new int[intnum];
    int numdiff = std::set_difference(intlist, intlist+intnum,
				      indices, indices+len, diff) - diff;
    if (numdiff < intnum) {
      ekk_deleteIntegerInformation(model);
      if (numdiff > 0) {
	ekk_addIntegerSet(model, 100, numdiff, diff, NULL, NULL);
      }
    }
    delete[] diff;
    delete[] sorted;
    ekk_free(intlist);
  }
}
//-----------------------------------------------------------------------------
void
OsiOslSolverInterface::setInteger(const int* indices, int len)
{
  ekk_addIntegerSet(getMutableModelPtr(), 100, len, indices, NULL, NULL);
}
//#############################################################################
void OsiOslSolverInterface::setObjSense(double s) 
{
  ekk_setRmaxmin(getMutableModelPtr(),s);
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::setColSolution(const double * cs) 
{
  EKKModel * m = getMutableModelPtr();
  ekk_copyColsol(m,cs,0,ekk_getInumcols(m));
}
//-----------------------------------------------------------------------------
void OsiOslSolverInterface::setRowPrice(const double * rs) 
{
  EKKModel * m = getMutableModelPtr();
  ekk_copyRowduals(m,rs,0,ekk_getInumrows(m));
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiOslSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  // *TEST*
  ekk_addOneColumn(getModelPtr(), obj, collb, colub,
		   vec.getNumElements(), vec.getIndices(), vec.getElements());
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  // *TEST*
  int i;
  int nz = 0;
  for (i = 0; i < numcols; ++i)
    nz += cols[i]->getNumElements();

  int* index = new int[nz];
  double* elem = new double[nz];
  int* start = new int[numcols+1];

  nz = 0;
  start[0] = 0;
  for (i = 0; i < numcols; ++i) {
    const CoinPackedVectorBase* col = cols[i];
    const int len = col->getNumElements();
    CoinDisjointCopyN(col->getIndices(), len, index+nz);
    CoinDisjointCopyN(col->getElements(), len, elem+nz);
    nz += len;
    start[i+1] = nz;
  }
  ekk_addColumns(getModelPtr(), numcols, obj, collb, colub,
		 start, index, elem);
  delete[] start;
  delete[] elem;
  delete[] index;
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  // *TEST*
  ekk_deleteColumns(getModelPtr(), num, columnIndices);
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  // *TEST*
  ekk_addOneRow(getModelPtr(), rowlb, rowub,
		vec.getNumElements(), vec.getIndices(), vec.getElements());
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  // *TEST*
  double rowlb, rowub;
  convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
  ekk_addOneRow(getModelPtr(), rowlb, rowub,
		vec.getNumElements(), vec.getIndices(), vec.getElements());
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  // *TEST*
  int i;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

  int* index = new int[nz];
  double* elem = new double[nz];
  int* start = new int[numrows+1];

  nz = 0;
  start[0] = 0;
  for (i = 0; i < numrows; ++i) {
    const CoinPackedVectorBase* row = rows[i];
    const int len = row->getNumElements();
    CoinDisjointCopyN(row->getIndices(), len, index+nz);
    CoinDisjointCopyN(row->getElements(), len, elem+nz);
    nz += len;
    start[i+1] = nz;
  }
  ekk_addRows(getModelPtr(), numrows, rowlb, rowub, start, index, elem);
  delete[] start;
  delete[] elem;
  delete[] index;
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  // *TEST*
  int i;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

  int* index = new int[nz];
  double* elem = new double[nz];
  int* start = new int[numrows+1];
  double* rowlb = new double[numrows];
  double* rowub = new double[numrows];

  nz = 0;
  start[0] = 0;
  for (i = 0; i < numrows; ++i) {
    convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rowlb[i], rowub[i]);
    const CoinPackedVectorBase* row = rows[i];
    const int len = row->getNumElements();
    CoinDisjointCopyN(row->getIndices(), len, index+nz);
    CoinDisjointCopyN(row->getElements(), len, elem+nz);
    nz += len;
    start[i+1] = nz;
  }
  ekk_addRows(getModelPtr(), numrows, rowlb, rowub, start, index, elem);
  delete[] rowub;
  delete[] rowlb;
  delete[] start;
  delete[] elem;
  delete[] index;
}
//-----------------------------------------------------------------------------
void 
OsiOslSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  // *TEST*
  ekk_deleteRows(getModelPtr(), num, rowIndices);
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiOslSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
   const CoinPackedMatrix * m = toColumnOrderedGapFree(matrix);
   loadProblem(m->getNumCols(), m->getNumRows(),
	       m->getVectorStarts(), m->getIndices(), m->getElements(),
	       collb, colub, obj, rowlb, rowub);
   if (m != &matrix)
      delete m;
}

//-----------------------------------------------------------------------------

void
OsiOslSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     double*& rowlb, double*& rowub)
{
   loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowlb;  rowlb = 0;
   delete[] rowub;  rowub = 0;
}

//-----------------------------------------------------------------------------

void
OsiOslSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
   const CoinPackedMatrix * m = toColumnOrderedGapFree(matrix);
   loadProblem(m->getNumCols(), m->getNumRows(),
	       m->getVectorStarts(), m->getIndices(), m->getElements(),
	       collb, colub, obj, rowsen, rowrhs, rowrng);
   if (m != &matrix)
      delete m;
}

//-----------------------------------------------------------------------------

void
OsiOslSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     char*& rowsen, double*& rowrhs,
				     double*& rowrng)
{
   loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowsen; rowsen = 0;
   delete[] rowrhs; rowrhs = 0;
   delete[] rowrng; rowrng = 0;
}

//-----------------------------------------------------------------------------

void
OsiOslSolverInterface::loadProblem(const int numcols, const int numrows,
				   const CoinBigIndex * start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
   // Use getMutableModelPtr(), so that the cached stuff is not immediately
   // deleted.
   ekk_loadRimModel(getMutableModelPtr(),
		    numrows, rowlb, rowub, numcols, obj, collb, colub);
   ekk_addColumnElementBlock(getMutableModelPtr(),
			     numcols, index, (const int *)start, value);
   // Now we can free the cached results
   freeCachedResults();
}
//-----------------------------------------------------------------------------

void
OsiOslSolverInterface::loadProblem(const int numcols, const int numrows,
				   const CoinBigIndex * start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
   assert( rowsen != NULL );
   assert( rowrhs != NULL );
   double * rowlb = new double[numrows];
   double * rowub = new double[numrows];
   for (int i = numrows-1; i >= 0; --i) {   
      convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
   }
   // Use getMutableModelPtr(), so that the cached stuff is not immediately
   // deleted.
   ekk_loadRimModel(getMutableModelPtr(),
		    numrows, rowlb, rowub, numcols, obj, collb, colub);
   ekk_addColumnElementBlock(getMutableModelPtr(),
			     numcols, index, (const int *) start, value);
   // Now we can free the cached results
   freeCachedResults();
   delete[] rowlb;
   delete[] rowub;
}

//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------

int OsiOslSolverInterface::readMps(const char *filename,
				   const char *extension)
{
   return OsiSolverInterface::readMps(filename, extension);
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------

void OsiOslSolverInterface::writeMps(const char * filename,
				     const char * extension,
				     double objSense) const
{
  std::string f(filename);
  std::string e(extension);
  std::string fullname;
  if (e!="") {
    fullname = f + "." + e;
  } else {
    // no extension so no trailing period
    fullname = f;
  }
  ekk_exportModel(getMutableModelPtr(), fullname.c_str(), 1, 2);
}

//#############################################################################
// OSL specific public interfaces
//#############################################################################

EKKModel * OsiOslSolverInterface::getModelPtr()
{
  freeCachedResults();
  return getMutableModelPtr();
}

//------------------------------------------------------------------- 

void OsiOslSolverInterface::incrementInstanceCounter()
{
  if ( numInstances_ == 0 ) {
    contextPtr_ = ekk_initializeContext();
    assert( contextPtr_ != NULL );    
    ekk_messagesPrintOff(ekk_baseModel(contextPtr_), 1, 5999);
    //FILE *oslMesgFILE = fopen("osl.log","a+");    
    //ekk_setLogfileFilePointer(ekk_baseModel(contextPtr_),stderr );    
    //ekk_messagesPrintOn(NULL ,1,2999);
  }
  numInstances_++;
}

//------------------------------------------------------------------- 

void OsiOslSolverInterface::decrementInstanceCounter()
{
  assert( getNumInstances() != 0 );
  numInstances_--;
  if ( numInstances_ == 0 ) {
    //FILE * oslMesgFile =
    //   ekk_getLogfileFilePointer(ekk_baseModel(contextPtr_));
    //fclose(oslMesgFile);
    ekk_endContext(contextPtr_);
    contextPtr_ = NULL;
  }
}

//------------------------------------------------------------------- 

unsigned int OsiOslSolverInterface::getNumInstances()
{
  return numInstances_;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiOslSolverInterface::OsiOslSolverInterface ()
  :
OsiSolverInterface(),
modelPtr_(NULL),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
matrixByRow_(NULL),
matrixByColumn_(NULL)
{
  incrementInstanceCounter();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * OsiOslSolverInterface::clone(bool CopyData) const
{
   if (CopyData) {
      return new OsiOslSolverInterface(*this);
   } else {
      return new OsiOslSolverInterface();
   }
}

#if 0
//-------------------------------------------------------------------
// Alternate Constructor 
//-------------------------------------------------------------------
OsiOslSolverInterface::OsiOslSolverInterface (EKKModel * m)
:
OsiSolverInterface(),
modelPtr_(m)
{
  // nothing to do here
}
#endif

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiOslSolverInterface::OsiOslSolverInterface (
                  const OsiOslSolverInterface & source)
:
OsiSolverInterface(source),
modelPtr_(NULL),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
matrixByRow_(NULL),
matrixByColumn_(NULL)
{
  incrementInstanceCounter();  
  if ( source.modelPtr_ !=NULL ) {
    ekk_copyModel(getModelPtr(),source.modelPtr_);
    // OK and proper to leave rowsense_, rhs_, and
    // rowrange_ to NULL.  They will be constructed
    // if they are required.
  }
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiOslSolverInterface::~OsiOslSolverInterface ()
{
  gutsOfDestructor();
  decrementInstanceCounter();
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiOslSolverInterface &
OsiOslSolverInterface::operator=(const OsiOslSolverInterface& rhs)
{
  if (this != &rhs) {    
    OsiSolverInterface::operator=(rhs);
    gutsOfDestructor();
    if ( rhs.modelPtr_ !=NULL ) {
      ekk_copyModel(getModelPtr(),rhs.modelPtr_);
    }
    delete ws_;
    ws_ = NULL;
  }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiOslSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  EKKModel * m = getModelPtr();
  const CoinPackedVector & row=rowCut.row();
  ekk_addOneRow(m, rowCut.lb(),rowCut.ub(),
		row.getNumElements(),row.getIndices(),row.getElements() ); 
}

//-----------------------------------------------------------------------------

void OsiOslSolverInterface::applyColCut( const OsiColCut & cc )
{
  EKKModel * m = getModelPtr();
  const double * oslColLB = ekk_collower(m);
  const double * oslColUB = ekk_colupper(m);
  const CoinPackedVector & lbs = cc.lbs();
  const CoinPackedVector & ubs = cc.ubs();
  int i;

  for ( i=0; i<lbs.getNumElements(); i++ ) {
    if ( lbs.getElements()[i] > oslColLB[lbs.getIndices()[i]] )
       //oslColLB[lbs.getIndices()[i]] = lbs.getElements()[i];
       ekk_copyCollower(m, &(lbs.getElements()[i]),
			lbs.getIndices()[i],lbs.getIndices()[i]+1);
  }
  for ( i=0; i<ubs.getNumElements(); i++ ) {
    if ( ubs.getElements()[i] < oslColUB[ubs.getIndices()[i]] )
      // oslColUB[ubs.getIndices()[i]] = ubs.getElements()[i];
      ekk_copyColupper(m, &(ubs.getElements()[i]),
		       ubs.getIndices()[i],ubs.getIndices()[i]+1);
  }
}

//#############################################################################
// Static methods and data
//#############################################################################

EKKContext * OsiOslSolverInterface::getContextPtr()
{
  assert( contextPtr_ != NULL );
  return contextPtr_;
}

//------------------------------------------------------------------- 

EKKContext * OsiOslSolverInterface::contextPtr_  = NULL;

//-------------------------------------------------------------------

unsigned int OsiOslSolverInterface::numInstances_ = 0;
 
//#############################################################################
// Private methods
//#############################################################################

//-------------------------------------------------------------------
// Get pointer to EKKModel.
// const methods should use getMutableModelPtr().
// non-const methods should use getModelPtr().
//------------------------------------------------------------------- 
EKKModel * OsiOslSolverInterface::getMutableModelPtr() const
{
  if ( modelPtr_ == NULL ) {
    modelPtr_ = ekk_newModel(getContextPtr(),NULL);
  }
  return modelPtr_;
}

//-------------------------------------------------------------------

#if 0
void OsiOslSolverInterface::gutsOfCopy( const OsiOslSolverInterface & source )
{
  modelPtr_ = source.modelPtr_;
}
#endif

//-------------------------------------------------------------------

void OsiOslSolverInterface::gutsOfDestructor()
{  
  if ( modelPtr_ != NULL ) {
    ekk_deleteModel(modelPtr_);
    modelPtr_=NULL;
    freeCachedResults();
  }
  assert( modelPtr_==NULL );
  assert( rowsense_==NULL );
  assert( rhs_==NULL );
  assert( rowrange_==NULL );
  assert( ws_==NULL );
  assert( matrixByRow_==NULL );
  assert( matrixByColumn_==NULL );
}

//------------------------------------------------------------------- 

void OsiOslSolverInterface::freeCachedResults()
{  
  delete [] rowsense_;
  delete [] rhs_;
  delete [] rowrange_;
  delete matrixByRow_;
  delete matrixByColumn_;
  delete ws_;
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
  matrixByRow_=NULL;
  matrixByColumn_=NULL;
  ws_ = NULL;
}

//------------------------------------------------------------------
void OsiOslSolverInterface::extractSenseRhsRange() const
{
  if (rowsense_ == NULL) {
    // all three must be NULL
    assert ((rhs_ == NULL) && (rowrange_ == NULL));
    
    EKKModel * m = getMutableModelPtr();
    
    int nr=ekk_getInumrows(m);
    if ( nr != 0 ) {
      rowsense_ = new char[nr];
      rhs_ = new double[nr];
      rowrange_ = new double[nr];
      std::fill(rowrange_,rowrange_+nr,0.0);
      
      const double * lb = ekk_rowlower(m);
      const double * ub = ekk_rowupper(m);
      
      int i;
      for ( i=0; i<nr; i++ ) {
        convertBoundToSense(lb[i], ub[i], rowsense_[i], rhs_[i], rowrange_[i]);
      }
    }
  }
}

//#############################################################################
// Resets as if default constructor
void 
OsiOslSolverInterface::reset()
{
  setInitialData(); // clear base class
  gutsOfDestructor();
  itlimOrig_=9999999;
}
