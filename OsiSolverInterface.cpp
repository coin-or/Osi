// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>

#include "CoinHelperFunctions.hpp"
#include "OsiWarmStart.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiMpsReader.hpp"
#include "OsiOsiMessage.hpp"

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
   OsiAbsFltEq eq(etol);
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
			    const OsiPackedVectorBase * const * cols,
			    const double* collb, const double* colub,   
			    const double* obj)
{
  for (int i = 0; i < numcols; ++i) {
    addCol(*cols[i], collb[i], colub[i], obj[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::addRows(const int numrows,
			    const OsiPackedVectorBase * const * rows,
			    const double* rowlb, const double* rowub)
{
  for (int i = 0; i < numrows; ++i) {
    addRow(*rows[i], rowlb[i], rowub[i]);
  }
}
//-----------------------------------------------------------------------------
void
OsiSolverInterface::addRows(const int numrows,
			    const OsiPackedVectorBase * const * rows,
			    const char* rowsen, const double* rowrhs,   
			    const double* rowrng)
{
  for (int i = 0; i < numrows; ++i) {
    addRow(*rows[i], rowsen[i], rowrhs[i], rowrng[i]);
  }
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
  appData_(NULL), 
  rowCutDebugger_(NULL),
  ws_(NULL),
  defaultHandler_(true)
{
  intParam_[OsiMaxNumIteration] = 9999999;
  intParam_[OsiMaxNumIterationHotStart] = 9999999;

  dblParam_[OsiDualObjectiveLimit] = DBL_MAX;
  dblParam_[OsiPrimalObjectiveLimit] = DBL_MAX;
  dblParam_[OsiDualTolerance] = 1e-6;
  dblParam_[OsiPrimalTolerance] = 1e-6;
  dblParam_[OsiObjOffset] = 0.0;

  strParam_[OsiProbName] = "OsiDefaultName";
  handler_ = new OsiMessageHandler();
  messages_ = OsiOsiMessage();
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiSolverInterface::OsiSolverInterface (const OsiSolverInterface & rhs) :
  appData_(rhs.appData_),
  rowCutDebugger_(NULL),
  ws_(NULL),
  defaultHandler_(true)
{  
  if ( rhs.rowCutDebugger_!=NULL )
    rowCutDebugger_ = new OsiRowCutDebugger(*rhs.rowCutDebugger_);
  defaultHandler_ = rhs.defaultHandler_;
  if (defaultHandler_)
    handler_ = new OsiMessageHandler(*rhs.handler_);
  else
    handler_ = rhs.handler_;
  messages_ = OsiOsiMessage();
  CoinDisjointCopyN(rhs.intParam_, OsiLastIntParam, intParam_);
  CoinDisjointCopyN(rhs.dblParam_, OsiLastDblParam, dblParam_);
  CoinDisjointCopyN(rhs.strParam_, OsiLastStrParam, strParam_);
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
    delete ws_;
    ws_ = NULL;
     if (defaultHandler_) {
      delete handler_;
      handler_ = NULL;
    }
    defaultHandler_ = rhs.defaultHandler_;
    if (defaultHandler_)
      handler_ = new OsiMessageHandler(*rhs.handler_);
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
  std::string f(filename);
  std::string e(extension);
  std::string fullname;
  if (e!="") {
    fullname = f + "." + e;
  } else {
    // no extension so no trailing period
    fullname = f;
  }
  OsiMpsReader m;
  int numberErrors;
  m.setInfinity(getInfinity());
  
  numberErrors = m.readMps(filename,extension);
  handler_->message(OSI_SOLVER_MPS,messages_)
    <<m.getProblemName()<< numberErrors <<OsiMessageEol;
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

// Pass in Message handler (not deleted at end)
void 
OsiSolverInterface::passInMessageHandler(OsiMessageHandler * handler)
{
  defaultHandler_=false;
  handler_=handler;
}
// Set language
void 
OsiSolverInterface::newLanguage(OsiMessages::Language language)
{
  messages_ = OsiOsiMessage(language);
}

// Function to return number in most efficient way
// Also creates row name field
/* formatType is
   0 - normal and 8 character names
   1 - extra accuracy
   2 - IEEE hex
   4 - normal but free format
*/
static void convertDouble(int formatType,double value,char outputValue[20],
			  const char * name, char outputRow[100])
{
  assert (formatType!=2);
  if ((formatType&3)==0) {
    if (value<1.0e20) {
      int power10;
      if (value>=0.0)
	power10 =(int) log10(value);
      else
	power10 =(int) log10(-value)+1;
      if (power10<5&&power10>-2) {
	char format[7];
	int decimal = 10-power10;
	if (decimal>10)
	  decimal=10;
	sprintf(format,"%%12.%df",decimal);
	sprintf(outputValue,format,value);
	// take off trailing 0
	int j;
	for (j=11;j>=0;j--) {
	  if (outputValue[j]=='0')
	    outputValue[j]=' ';
	  else
	    break;
	}
      } else {
	sprintf(outputValue,"%12g",value);
      }
    } else {
      outputValue[0]= '\0'; // needs no value
    }
  } else {
    if (value<1.0e20) {
      sprintf(outputValue,"%19g",value);
      // take out blanks
      int i=0;
      int j;
      for (j=0;j<19;j++) {
	if (outputValue[j]!=' ')
	  outputValue[i++]=outputValue[j];
      }
      outputValue[i]='\0';
    } else {
      outputValue[0]= '\0'; // needs no value
    }
  }
  strcpy(outputRow,name);
  if (!formatType) {
    int i;
    // pad out to 12 and 8
    for (i=0;i<12;i++) {
      if (outputValue[i]=='\0')
	break;
    }
    for (;i<12;i++) 
      outputValue[i]=' ';
    outputValue[12]='\0';
    for (i=0;i<8;i++) {
      if (outputRow[i]=='\0')
	break;
    }
    for (;i<8;i++) 
      outputRow[i]=' ';
    outputRow[8]='\0';
  }
}
// Put out card image
static void outputCard(int formatType,int numberFields,
		       FILE *fp, std::string head,
		       const char * name,
		       const char outputValue[2][20],
		       const char outputRow[2][100])
{
  fprintf(fp,"%s",head.c_str());
  int i;
  if (!formatType) {
    char outputColumn[9];
    strcpy(outputColumn,name);
    for (i=0;i<8;i++) {
      if (outputColumn[i]=='\0')
	break;
    }
    for (;i<8;i++) 
      outputColumn[i]=' ';
    outputColumn[8]='\0';
    fprintf(fp,"%s  ",outputColumn);
    for (i=0;i<numberFields;i++) {
      fprintf(fp,"%s  %s",outputRow[i],outputValue[i]);
      if (i<numberFields-1)
	fprintf(fp,"   ");
    }
  } else {
    fprintf(fp,"%s",name);
    for (i=0;i<numberFields;i++) {
      fprintf(fp," %s %s",outputRow[i],outputValue[i]);
    }
  }
  fprintf(fp,"\n");
}
/* Write the problem into an mps file of the given filename,
   names may be null.  formatType is
   0 - normal
   1 - extra accuracy
   2 - IEEE hex
*/
int 
OsiSolverInterface::writeMps(const char *filename, 
			     const char ** rowNames, 
			     const char ** columnNames,
			     int formatType,
			     int numberAcross) const
{
  FILE * fp = fopen(filename,"w");
  if (!fp)
    return -1;
  int numberRows = getNumRows();
  int numberColumns = getNumCols();
  char ** rowNameTemp=NULL;
  char ** columnNameTemp=NULL;
  // If long names free format
  unsigned int length = 8;
  bool freeFormat=(formatType!=0);
  int i;
  if (rowNames) {
    for (i=0;i<numberRows;i++) {
      if (strlen(rowNames[i])>length)
	length = strlen(rowNames[i]);
    }
  } else {
    rowNameTemp = new char * [numberRows];
    for (i=0;i<numberRows;i++) {
      char * name = new char[9];
      sprintf(name,"R%7.7d",i);
      rowNameTemp[i]=name;
    }
    rowNames=(const char **) rowNameTemp;
  }
  if (columnNames) {
    for (i=0;i<numberColumns;i++) {
      if (strlen(columnNames[i])>length)
	length = strlen(columnNames[i]);
    }
  } else {
    columnNameTemp = new char * [numberColumns];
    for (i=0;i<numberColumns;i++) {
      char * name = new char[9];
      sprintf(name,"C%7.7d",i);
      columnNameTemp[i]=name;
    }
    columnNames=(const char **) columnNameTemp;
  }
  if (length>8&&!freeFormat) {
    freeFormat = true;
    formatType=4;
  }

  // NAME card

  fprintf(fp,"NAME          ");
  std::string problemName;
  getStrParam(OsiProbName, problemName);
  if (!problemName.size()) 
    problemName="BLANK   ";
  if (problemName.size()>=8) {
    for (i=0;i<8;i++) 
      fprintf(fp,"%c",problemName[i]);
  } else {
    unsigned int i;
    for (i=0;i<problemName.size();i++) 
      fprintf(fp,"%c",problemName[i]);
    for (;i<8;i++) 
      fprintf(fp,"%c",' ');
  }

  if (freeFormat)
    fprintf(fp,"  FREE");

  // finish off name and do ROWS card and objective 

  fprintf(fp,"\nROWS\n N  OBJROW\n");

  // Rows section
  // Sense array
  const char * sense = getRowSense();
  
  for (i=0;i<numberRows;i++) {
    if (sense[i]!='R')
      fprintf(fp," %c  %s\n",sense[i],rowNames[i]);
    else
      fprintf(fp," L  %s\n",rowNames[i]);
  }
  
  // COLUMNS card
  fprintf(fp,"COLUMNS\n");

  bool ifBounds=false;
  double largeValue = getInfinity();

  const double * columnLower = getColLower();
  const double * columnUpper = getColUpper();
  const double * objective = getObjCoefficients();
  const OsiPackedMatrix * matrix = getMatrixByCol();
  const double * elements = matrix->getElements();
  const int * rows = matrix->getIndices();
  const int * starts = matrix->getVectorStarts();
  const int * lengths = matrix->getVectorLengths();

  char outputValue[2][20];
  char outputRow[2][100];

  // Through columns (only put out if elements or objective value)
  for (i=0;i<numberColumns;i++) {
    if (objective[i]||lengths[i]) {
      // see if bound will be needed
      if (columnLower[i]||columnUpper[i]<largeValue)
	ifBounds=true;
      int numberFields=0;
      if (objective[i]) {
	convertDouble(formatType,objective[i],outputValue[0],
		      "OBJROW",outputRow[0]);
	numberFields=1;
      }
      if (numberFields==numberAcross) {
	// put out card
	outputCard(formatType, numberFields,
		   fp, "    ",
		   columnNames[i],
		   outputValue,
		   outputRow);
	numberFields=0;
      }
      int j;
      for (j=0;j<lengths[i];j++) {
	convertDouble(formatType,elements[starts[i]+j],
		      outputValue[numberFields],
		      rowNames[rows[starts[i]+j]],
		      outputRow[numberFields]);
	numberFields++;
	if (numberFields==numberAcross) {
	  // put out card
	  outputCard(formatType, numberFields,
		     fp, "    ",
		     columnNames[i],
		     outputValue,
		     outputRow);
	  numberFields=0;
	}
      }
      if (numberFields) {
	// put out card
	outputCard(formatType, numberFields,
		   fp, "    ",
		   columnNames[i],
		   outputValue,
		   outputRow);
      }
    }
  }
  

  bool ifRange=false;
  // RHS
  fprintf(fp,"RHS\n");

  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();
  
  int numberFields = 0;
  for (i=0;i<numberRows;i++) {
    double value;
    switch (sense[i]) {
    case 'E':
      value=rowLower[i];
      break;
    case 'R':
      value=rowUpper[i];
      ifRange=true;
      break;
    case 'L':
      value=rowUpper[i];
      break;
    case 'G':
      value=rowLower[i];
      break;
    default:
      value=0.0;
      break;
    }
    if (value) {
      convertDouble(formatType,value,
		    outputValue[numberFields],
		    rowNames[i],
		    outputRow[numberFields]);
      numberFields++;
      if (numberFields==numberAcross) {
	// put out card
	outputCard(formatType, numberFields,
		   fp, "    ",
		   "RHS",
		   outputValue,
		   outputRow);
	numberFields=0;
      }
    }
  }
  if (numberFields) {
    // put out card
    outputCard(formatType, numberFields,
	       fp, "    ",
	       "RHS",
	       outputValue,
	       outputRow);
  }
  

  if (ifRange) {
    // RANGE
    fprintf(fp,"RANGE\n");

    numberFields = 0;
    for (i=0;i<numberRows;i++) {
      if (sense[i]=='R') {
	double value =rowUpper[i]-rowLower[i];
	convertDouble(formatType,value,
		      outputValue[numberFields],
		      rowNames[i],
		      outputRow[numberFields]);
	numberFields++;
	if (numberFields==numberAcross) {
	  // put out card
	  outputCard(formatType, numberFields,
		     fp, "    ",
		     "RANGE",
		     outputValue,
		     outputRow);
	  numberFields=0;
	}
      }
    }
    if (numberFields) {
      // put out card
      outputCard(formatType, numberFields,
		 fp, "    ",
		 "RANGE",
		 outputValue,
		 outputRow);
    }
  }
  

  if (ifBounds) {

    // BOUNDS
    fprintf(fp,"BOUNDS\n");

    for (i=0;i<numberColumns;i++) {
      if (objective[i]||lengths[i]) {
	// see if bound will be needed
	if (columnLower[i]||columnUpper[i]<largeValue) {
	  int numberFields=1;
	  std::string header[2];
	  double value[2];
	  if (columnLower[i]<=-largeValue) {
	    // FR or MI
	    if (columnUpper[i]>=largeValue) {
	      header[0]=" FR ";
	      value[0] = largeValue;
	    } else {
	      header[0]=" MI ";
	      value[0] = largeValue;
	      header[0]=" UP ";
	      value[0] = columnUpper[i];
	      numberFields=2;
	    }
	  } else if (fabs(columnUpper[i]-columnLower[i])<1.0e-8) {
	    header[0]=" FX ";
	    value[0] = columnLower[i];
	  } else {
	    // do LO if needed
	    if (columnLower[i]) {
	      // LO
	      header[0]=" LO ";
	      value[0] = columnLower[i];
	      if (isInteger(i)) {
		// Integer variable so UI
		header[1]=" UI ";
		value[1] = columnUpper[i];
		numberFields=2;
	      } else if (columnUpper[i]<largeValue) {
		// UP
		header[1]=" UP ";
		value[1] = columnUpper[i];
		numberFields=2;
	      }
	    } else {
	      if (isInteger(i)) {
		// Integer variable so BV or UI
		if (fabs(columnUpper[i]-1.0)<1.0e-8) {
		  // BV
		  header[0]=" BV ";
		  value[0] = largeValue;
		} else {
		  // UI
		  header[0]=" UI ";
		  value[0] = columnUpper[i];
		}
	      } else {
		// UP
		header[0]=" UP ";
		value[0] = columnUpper[i];
	      }
	    }
	  }
	  // put out fields
	  int j;
	  for (j=0;j<numberFields;j++) {
	    convertDouble(formatType,value[j],
			  outputValue[0],
			  columnNames[i],
			  outputRow[0]);
	    // put out card
	    outputCard(formatType, 1,
		       fp, header[j],
		       "BOUND",
		       outputValue,
		       outputRow);
	  }
	}
      }
    }
  }

  // and finish

  fprintf(fp,"ENDATA\n");

  fclose(fp);

  if (rowNameTemp) {
    for (i=0;i<numberRows;i++) {
      delete [] rowNameTemp[i];
    }
    delete [] rowNameTemp;
  }
  if (columnNameTemp) {
    for (i=0;i<numberColumns;i++) {
      delete [] columnNameTemp[i];
    }
    delete [] columnNameTemp;
  }
  return 0;
}
