// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>

#include "CoinHelperFunctions.hpp"
#include "CoinMpsIO.hpp"
#include "CoinMessage.hpp"
#include "CoinWarmStart.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include <cassert>

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
// It is assumed that only bounds (on vars/constraints) can change between
// markHotStart() and unmarkHotStart()
//#############################################################################

void OsiSolverInterface::markHotStart()
{
  delete ws_;
  ws_ = getWarmStart();
}

void OsiSolverInterface::solveFromHotStart()
{
  setWarmStart(ws_);
  resolve();
}

void OsiSolverInterface::unmarkHotStart()
{
  delete ws_;
  ws_ = NULL;
}

//#############################################################################
// Get indices of solution vector which are integer variables presently at
// fractional values
//#############################################################################

OsiVectorInt
OsiSolverInterface::getFractionalIndices(const double etol) const
{
   const int colnum = getNumCols();
   OsiVectorInt frac;
   CoinAbsFltEq eq(etol);
   for (int i = 0; i < colnum; ++i) {
      if (isInteger(i)) {
	 const double ci = getColSolution()[i];
	 const double distanceFromInteger = ci - floor(ci + 0.5);
	 if (! eq(distanceFromInteger, 0.0))
	    frac.push_back(i);
      }
   }
   return frac;
}


int OsiSolverInterface::getNumElements() const
{
  return getMatrixByRow()->getNumElements();
}

//#############################################################################
// Methods for determining the type of column variable.
// The method isContinuous() is presently implemented in the derived classes.
// The methods below depend on isContinuous and the values
// stored in the column bounds.
//#############################################################################

bool 
OsiSolverInterface::isBinary(int colIndex) const
{
  if ( isContinuous(colIndex) ) return false; 
  const double * cu = getColUpper();
  const double * cl = getColLower();
  if (
    (cu[colIndex]== 1 || cu[colIndex]== 0) && 
    (cl[colIndex]== 0 || cl[colIndex]==1)
    ) return true;
  else return false;
}
//-----------------------------------------------------------------------------
bool 
OsiSolverInterface::isInteger(int colIndex) const
{
   return !isContinuous(colIndex);
}
//-----------------------------------------------------------------------------
bool 
OsiSolverInterface::isIntegerNonBinary(int colIndex) const
{
  if ( isInteger(colIndex) && !isBinary(colIndex) )
    return true; 
  else return false;
}
//-----------------------------------------------------------------------------
bool 
OsiSolverInterface::isFreeBinary(int colIndex) const
{
  if ( isContinuous(colIndex) ) return false;
  const double * cu = getColUpper();
  const double * cl = getColLower();
  if (
    (cu[colIndex]== 1) &&
    (cl[colIndex]== 0)
    ) return true;
  else return false;
}

//#############################################################################
// Built-in (i.e., slow) methods for problem modification
//#############################################################################

void
OsiSolverInterface::setObjCoeffSet(const int* indexFirst,
				  const int* indexLast,
				  const double* coeffList)
{
   const int cnt = indexLast - indexFirst;
   for (int i = 0; i < cnt; ++i) {
      setObjCoeff(indexFirst[i], coeffList[i]);
   }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::setColSetBounds(const int* indexFirst,
				    const int* indexLast,
				    const double* boundList)
{
  while (indexFirst != indexLast) {
    setColBounds(*indexFirst, boundList[0], boundList[1]);
    ++indexFirst;
    boundList += 2;
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::setRowSetBounds(const int* indexFirst,
				    const int* indexLast,
				    const double* boundList)
{
  while (indexFirst != indexLast) {
    setRowBounds(*indexFirst, boundList[0], boundList[1]);
    ++indexFirst;
    boundList += 2;
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::setRowSetTypes(const int* indexFirst,
				   const int* indexLast,
				   const char* senseList,
				   const double* rhsList,
				   const double* rangeList)
{
  while (indexFirst != indexLast) {
    setRowType(*indexFirst++, *senseList++, *rhsList++, *rangeList++);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::setContinuous(const int* indices, int len)
{
  for (int i = 0; i < len; ++i) {
    setContinuous(indices[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::setInteger(const int* indices, int len)
{
  for (int i = 0; i < len; ++i) {
    setInteger(indices[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::addCols(const int numcols,
			    const CoinPackedVectorBase * const * cols,
			    const double* collb, const double* colub,   
			    const double* obj)
{
  for (int i = 0; i < numcols; ++i) {
    addCol(*cols[i], collb[i], colub[i], obj[i]);
  }
}
/* Add a column (primal variable) to the problem. */
void 
OsiSolverInterface::addCol(int numberElements, const int * rows, const double * elements,
			   const double collb, const double colub,   
			   const double obj) 
{
  CoinPackedVector column(numberElements, rows, elements);
  addCol(column,collb,colub,obj);
}
/* Add a set of columns (primal variables) to the problem.
   
This default implementation simply makes repeated calls to addCol().
*/
void OsiSolverInterface::addCols(const int numcols,
				 const int * columnStarts, const int * rows, const double * elements,
				 const double* collb, const double* colub,   
				 const double* obj)
{
  double infinity = getInfinity();
  for (int i = 0; i < numcols; ++i) {
    int start = columnStarts[i];
    int number = columnStarts[i+1]-start;
    assert (number>=0);
    addCol(number, rows+start, elements+start, collb ? collb[i] : 0.0, 
	   colub ? colub[i] : infinity, 
	   obj ? obj[i] : 0.0);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::addRows(const int numrows,
			    const CoinPackedVectorBase * const * rows,
			    const double* rowlb, const double* rowub)
{
  for (int i = 0; i < numrows; ++i) {
    addRow(*rows[i], rowlb[i], rowub[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::addRows(const int numrows,
			    const CoinPackedVectorBase * const * rows,
			    const char* rowsen, const double* rowrhs,   
			    const double* rowrng)
{
  for (int i = 0; i < numrows; ++i) {
    addRow(*rows[i], rowsen[i], rowrhs[i], rowrng[i]);
  }
}
/* Add a row (constraint) to the problem. */
void 
OsiSolverInterface::addRow(int numberElements, const int * columns, const double * elements,
			   const double rowlb, const double rowub) 
{
  CoinPackedVector row(numberElements, columns, elements);
  addRow(row,rowlb,rowub);
}
/* Add a set of rows (constraints) to the problem.
   
The default implementation simply makes repeated calls to addRow().
*/
void 
OsiSolverInterface::addRows(const int numrows,
			    const int * rowStarts, const int * columns, const double * elements,
			    const double* rowlb, const double* rowub)
{
  double infinity = getInfinity();
  for (int i = 0; i < numrows; ++i) {
    int start = rowStarts[i];
    int number = rowStarts[i+1]-start;
    assert (number>=0);
    addRow(number, columns+start, elements+start, rowlb ? rowlb[i] : -infinity, 
	   rowub ? rowub[i] : infinity);
  }
}


//#############################################################################
// Implement getObjValue in a simple way that the derived solver interfaces
// can use if the choose.
//#############################################################################
double OsiSolverInterface::getObjValue() const
{
  int nc = getNumCols();
  const double * objCoef = getObjCoefficients();
  const double * colSol  = getColSolution();
  double objOffset=0.0;
  getDblParam(OsiObjOffset,objOffset);
  
  // Compute dot product of objCoef and colSol and then adjust by offset
  double retVal = CoinPackedVector(nc,objCoef).dotProduct(colSol)-objOffset;
  return retVal;
}

//#############################################################################
// Apply Cuts
//#############################################################################

OsiSolverInterface::ApplyCutsReturnCode
OsiSolverInterface::applyCuts( const OsiCuts & cs, double effectivenessLb ) 
{
  OsiSolverInterface::ApplyCutsReturnCode retVal;
  int i;

  // Loop once for each column cut
  for ( i=0; i<cs.sizeColCuts(); i ++ ) {
    if ( cs.colCut(i).effectiveness() < effectivenessLb ) {
      retVal.incrementIneffective();
      continue;
    }
    if ( !cs.colCut(i).consistent() ) {
      retVal.incrementInternallyInconsistent();
      continue;
    }
    if ( !cs.colCut(i).consistent(*this) ) {
      retVal.incrementExternallyInconsistent();
      continue;
    }
    if ( cs.colCut(i).infeasible(*this) ) {
      retVal.incrementInfeasible();
      continue;
    }
    applyColCut( cs.colCut(i) );
    retVal.incrementApplied();
  }

  // Loop once for each row cut
  for ( i=0; i<cs.sizeRowCuts(); i ++ ) {
    if ( cs.rowCut(i).effectiveness() < effectivenessLb ) {
      retVal.incrementIneffective();
      continue;
    }
    if ( !cs.rowCut(i).consistent() ) {
      retVal.incrementInternallyInconsistent();
      continue;
    }
    if ( !cs.rowCut(i).consistent(*this) ) {
      retVal.incrementExternallyInconsistent();
      continue;
    }
    if ( cs.rowCut(i).infeasible(*this) ) {
      retVal.incrementInfeasible();
      continue;
    }
    applyRowCut( cs.rowCut(i) );
    retVal.incrementApplied();
  }
  
  return retVal;
}
/* Apply a collection of row cuts which are all effective.
   applyCuts seems to do one at a time which seems inefficient.
   The default does slowly, but solvers can override.
*/
void 
OsiSolverInterface::applyRowCuts(int numberCuts, const OsiRowCut * cuts)
{
  int i;
  for (i=0;i<numberCuts;i++) {
    applyRowCut(cuts[i]);
  }
}

//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// from the solverInterface.
// This field is the application to optionally define and use.
//#############################################################################

void OsiSolverInterface::setApplicationData(void * appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void * OsiSolverInterface::getApplicationData() const
{
  return appData_;
}

//#############################################################################
// Methods related to Row Cut Debuggers
//#############################################################################

//-------------------------------------------------------------------
// Activate Row Cut Debugger<br>
// If the model name passed is on list of known models
// then all cuts are checked to see that they do NOT cut
// off the known optimal solution.  

void OsiSolverInterface::activateRowCutDebugger (const char * modelName)
{
  delete rowCutDebugger_;
  rowCutDebugger_ = new OsiRowCutDebugger(*this,modelName);
}

//-------------------------------------------------------------------
// Get Row Cut Debugger<br>
// If there is a row cut debugger object associated with
// model AND if the known optimal solution is within the
// current feasible region then a pointer to the object is
// returned which may be used to test validity of cuts.
// Otherwise NULL is returned

const OsiRowCutDebugger * OsiSolverInterface::getRowCutDebugger() const
{
  if (rowCutDebugger_&&rowCutDebugger_->onOptimalPath(*this)) {
    return rowCutDebugger_;
  } else {
    return NULL;
  }
}

//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiSolverInterface::OsiSolverInterface () :
  rowCutDebugger_(NULL),
  appData_(NULL), 
  ws_(NULL),
  handler_(NULL),
  defaultHandler_(true)
{
  setInitialData();
}
// Set data for default constructor
void 
OsiSolverInterface::setInitialData()
{
  delete rowCutDebugger_;
  rowCutDebugger_ = NULL;
  delete ws_;
  ws_ = NULL;
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  defaultHandler_=true;
  appData_=NULL;
  intParam_[OsiMaxNumIteration] = 9999999;
  intParam_[OsiMaxNumIterationHotStart] = 9999999;

  dblParam_[OsiDualObjectiveLimit] = DBL_MAX;
  dblParam_[OsiPrimalObjectiveLimit] = DBL_MAX;
  dblParam_[OsiDualTolerance] = 1e-6;
  dblParam_[OsiPrimalTolerance] = 1e-6;
  dblParam_[OsiObjOffset] = 0.0;

  strParam_[OsiProbName] = "OsiDefaultName";
  strParam_[OsiSolverName] = "Unknown Solver";
  handler_ = new CoinMessageHandler();
  messages_ = CoinMessage();

  // initialize all hints
  int hint;
  for (hint=OsiDoPresolveInInitial;hint<OsiLastHintParam;hint++) {
    hintParam_[hint] = false;
    hintStrength_[hint] = OsiHintIgnore;
  }
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiSolverInterface::OsiSolverInterface (const OsiSolverInterface & rhs) :
  rowCutDebugger_(NULL),
  appData_(rhs.appData_),
  ws_(NULL),
  defaultHandler_(true)
{  
  if ( rhs.rowCutDebugger_!=NULL )
    rowCutDebugger_ = new OsiRowCutDebugger(*rhs.rowCutDebugger_);
  defaultHandler_ = rhs.defaultHandler_;
  if (defaultHandler_)
    handler_ = new CoinMessageHandler(*rhs.handler_);
  else
    handler_ = rhs.handler_;
  messages_ = CoinMessages(rhs.messages_);
  CoinDisjointCopyN(rhs.intParam_, OsiLastIntParam, intParam_);
  CoinDisjointCopyN(rhs.dblParam_, OsiLastDblParam, dblParam_);
  CoinDisjointCopyN(rhs.strParam_, OsiLastStrParam, strParam_);
  CoinDisjointCopyN(rhs.hintParam_, OsiLastHintParam, hintParam_);
  CoinDisjointCopyN(rhs.hintStrength_, OsiLastHintParam, hintStrength_);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiSolverInterface::~OsiSolverInterface ()
{
  // delete debugger - should be safe as only ever returned const
  delete rowCutDebugger_;
  rowCutDebugger_ = NULL;
  delete ws_;
  ws_ = NULL;
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiSolverInterface &
OsiSolverInterface::operator=(const OsiSolverInterface& rhs)
{
  if (this != &rhs) {
    appData_ = rhs.appData_;
    delete rowCutDebugger_;
    if ( rhs.rowCutDebugger_!=NULL )
      rowCutDebugger_ = new OsiRowCutDebugger(*rhs.rowCutDebugger_);
    else
      rowCutDebugger_ = NULL;
    CoinDisjointCopyN(rhs.intParam_, OsiLastIntParam, intParam_);
    CoinDisjointCopyN(rhs.dblParam_, OsiLastDblParam, dblParam_);
    CoinDisjointCopyN(rhs.strParam_, OsiLastStrParam, strParam_);
    CoinDisjointCopyN(rhs.hintParam_, OsiLastHintParam, hintParam_);
    CoinDisjointCopyN(rhs.hintStrength_, OsiLastHintParam, hintStrength_);
    delete ws_;
    ws_ = NULL;
     if (defaultHandler_) {
      delete handler_;
      handler_ = NULL;
    }
    defaultHandler_ = rhs.defaultHandler_;
    if (defaultHandler_)
      handler_ = new CoinMessageHandler(*rhs.handler_);
    else
      handler_ = rhs.handler_;
 }
  return *this;
}

//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------

int OsiSolverInterface::readMps(const char * filename,
    const char * extension)
{
  CoinMpsIO m;
  m.setInfinity(getInfinity());
  
  int numberErrors = m.readMps(filename,extension);
  handler_->message(COIN_SOLVER_MPS,messages_)
    <<m.getProblemName()<< numberErrors <<CoinMessageEol;
  if (!numberErrors) {

    // set objective function offest
    setDblParam(OsiObjOffset,m.objectiveOffset());

    // set problem name
    setStrParam(OsiProbName,m.getProblemName());

    // no errors
    loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		m.getObjCoefficients(),m.getRowSense(),m.getRightHandSide(),
		m.getRowRange());
    const char * integer = m.integerColumns();
    if (integer) {
      int i,n=0;
      int nCols=m.getNumCols();
      int * index = new int [nCols];
      for (i=0;i<nCols;i++) {
	if (integer[i]) {
	  index[n++]=i;
	}
      }
      setInteger(index,n);
      delete [] index;
    }
  }
  return numberErrors;
}
/* Read a problem in MPS format from the given full filename.
   
This uses CoinMpsIO::readMps() to read
the MPS file and returns the number of errors encountered.
It also may return an array of set information
*/
int 
OsiSolverInterface::readMps(const char *filename, const char*extension,
			    int & numberSets, CoinSet ** & sets)
{
  CoinMpsIO m;
  m.setInfinity(getInfinity());
  
  int numberErrors = m.readMps(filename,extension,numberSets,sets);
  handler_->message(COIN_SOLVER_MPS,messages_)
    <<m.getProblemName()<< numberErrors <<CoinMessageEol;
  if (!numberErrors) {

    // set objective function offest
    setDblParam(OsiObjOffset,m.objectiveOffset());

    // set problem name
    setStrParam(OsiProbName,m.getProblemName());

    // no errors
    loadProblem(*m.getMatrixByCol(),m.getColLower(),m.getColUpper(),
		m.getObjCoefficients(),m.getRowSense(),m.getRightHandSide(),
		m.getRowRange());
    const char * integer = m.integerColumns();
    if (integer) {
      int i,n=0;
      int nCols=m.getNumCols();
      int * index = new int [nCols];
      for (i=0;i<nCols;i++) {
	if (integer[i]) {
	  index[n++]=i;
	}
      }
      setInteger(index,n);
      delete [] index;
    }
  }
  return numberErrors;
}

int 
OsiSolverInterface::writeMpsNative(const char *filename, 
			     const char ** rowNames, 
			     const char ** columnNames,
			     int formatType,
			     int numberAcross,
			     double objSense) const
{
   const int numcols = getNumCols();
   char* integrality = new char[numcols];
   bool hasInteger = false;
   for (int i = 0; i < numcols; ++i) {
      if (isInteger(i)) {
	 integrality[i] = 1;
	 hasInteger = true;
      } else {
	 integrality[i] = 0;
      }
   }

   // Get multiplier for objective function - default 1.0
   double * objective = new double[numcols];
   memcpy(objective,getObjCoefficients(),numcols*sizeof(double));
   if (objSense*getObjSense()<0.0) {
     for (int i = 0; i < numcols; ++i) 
       objective [i] = - objective[i];
   }

   CoinMpsIO writer;
   writer.setInfinity(getInfinity());
   writer.passInMessageHandler(handler_);
   writer.setMpsData(*getMatrixByCol(), getInfinity(),
		     getColLower(), getColUpper(),
		     objective, hasInteger ? integrality : 0,
		     getRowLower(), getRowUpper(),
		     (const char **)0, (const char **)0);
   delete [] objective;
   delete[] integrality;
   return writer.writeMps(filename, 1 /*gzip it*/, formatType, numberAcross);
}

// Pass in Message handler (not deleted at end)
void 
OsiSolverInterface::passInMessageHandler(CoinMessageHandler * handler)
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  defaultHandler_=false;
  handler_=handler;
}
// Set language
void 
OsiSolverInterface::newLanguage(CoinMessages::Language language)
{
  messages_ = CoinMessage(language);
}
// copy all parameters in this section from one solver to another
void 
OsiSolverInterface::copyParameters(OsiSolverInterface & rhs)
{
  appData_ = rhs.appData_;
  if ( rhs.rowCutDebugger_!=NULL )
    rowCutDebugger_ = new OsiRowCutDebugger(*rhs.rowCutDebugger_);
  else
    rowCutDebugger_ = NULL;
  defaultHandler_ = rhs.defaultHandler_;
  if (defaultHandler_)
    handler_ = new CoinMessageHandler(*rhs.handler_);
  else
    handler_ = rhs.handler_;
   CoinDisjointCopyN(rhs.intParam_, OsiLastIntParam, intParam_);
  CoinDisjointCopyN(rhs.dblParam_, OsiLastDblParam, dblParam_);
  CoinDisjointCopyN(rhs.strParam_, OsiLastStrParam, strParam_);
  CoinDisjointCopyN(rhs.hintParam_, OsiLastHintParam, hintParam_);
  CoinDisjointCopyN(rhs.hintStrength_, OsiLastHintParam, hintStrength_);
}
// Resets as if default constructor
void 
OsiSolverInterface::reset()
{
  // Throw an exception
  throw CoinError("Needs coding for this interface", "reset",
		  "OsiSolverInterface");
}
