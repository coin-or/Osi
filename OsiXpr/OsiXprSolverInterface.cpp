// copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <numeric>
#include <strstream>

#define __ANSIC_
#include <xpresso.h>
#undef  __ANSIC_

#include "OsiXprSolverInterface.hpp"

// xpresso.h is redefining range in away that
// is in conflict with other usages
#define rangeTemp range
#undef range
#include "CoinHelperFunctions.hpp"
#include "OsiCuts.hpp"
#include "OsiColCut.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiRowCut.hpp"
#include "OsiWarmStartBasis.hpp"

// Let Xpress have its own definition of range again
#define range rangeTemp

//#############################################################################
// Solve methods
//#############################################################################

void
OsiXprSolverInterface::initialSolve(){
   activateMe();

   freeSolution();

   if ( objsense_ == 1.0 ) {
     if (getLogFilePtr()!=NULL) {
       fprintf(getLogFilePtr(),"minim(\"l\");\n");
     }
      minim("l");
   }
   else if ( objsense_ == -1.0 ) {
     if (getLogFilePtr()!=NULL) {
       fprintf(getLogFilePtr(),"maxim(\"l\");\n");
     }
      maxim("l");
   }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::resolve(){
   activateMe();

   freeSolution();

   if ( objsense_ == 1.0 ) {
     if (getLogFilePtr()!=NULL) {
       fprintf(getLogFilePtr(),"minim(\"dl\");\n");
     }
      minim("dl");
   }
   else if ( objsense_ == -1.0 ) {
     if (getLogFilePtr()!=NULL) {
       fprintf(getLogFilePtr(),"maxim(\"dl\");\n");
     }
      maxim("dl");
   }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::branchAndBound(){
  activateMe();
  
  freeSolution();

  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"global;\n");
  }
  global();
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiXprSolverInterface::setIntParam(OsiIntParam key, int value)
{
  bool retval = false;
  
  switch (key) {
  case OsiMaxNumIteration:
    retval = seticv(N_ITRLIM, value) == 0;
    break;
  case OsiMaxNumIterationHotStart:
    retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::setDblParam(OsiDblParam key, double value)
{
  bool retval = false;
  switch (key) {
  case OsiDualObjectiveLimit:
    retval = setdcv(N_CUTOFF, value) == 0;
    break;
  case OsiPrimalObjectiveLimit:
    retval = false;	
    break;
  case OsiDualTolerance:
    retval = false;	
    break;
  case OsiPrimalTolerance:
    retval = false;	
    break;
  case OsiObjOffset: 
    return OsiSolverInterface::setDblParam(key, value);
  case OsiLastDblParam:
    retval = false;
    break;
  }
  return retval;
}
//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  bool retval=false;
  switch (key) {

  case OsiProbName:
    OsiSolverInterface::setStrParam(key,value);
    return retval = true;

  case OsiLastStrParam:
    return false;
  }
  return false;
}
//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  bool retval = false;

  switch (key) {
  case OsiMaxNumIteration:
    retval = geticv(N_ITRLIM, &value)==1 ? true : false;
    break;
  case OsiMaxNumIterationHotStart:
    retval = false;
    break;
  case OsiLastIntParam:
    retval = false;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  bool retval = false;

  switch (key) {
  case OsiDualObjectiveLimit:
    retval = getdcv(N_CUTOFF, &value)==1 ? true : false;
    break;
  case OsiPrimalObjectiveLimit:
    retval = false;
    break;
  case OsiDualTolerance:
    retval = false;
    break;
  case OsiPrimalTolerance:
    retval = false;
    break;
  case OsiObjOffset:
    retval = OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiLastDblParam:
    retval = false;
    break;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  switch (key) {
  case OsiProbName:
    OsiSolverInterface::getStrParam(key, value);
    break;
  case OsiLastStrParam:
    return false;
  }
  return true;
}
//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiXprSolverInterface::isAbandoned() const
{
  activateMe();

  int  status, glstat;

  getipv(N_STATUS, &status);
  getipv(N_GLSTAT, &glstat);

  return status == 4 || glstat == 3 || glstat == 4;
  // LP unfinished || global search incomplete -- no int sol ||
  // global search incomplete -- int sol found
}

bool OsiXprSolverInterface::isProvenOptimal() const
{
  activateMe();

  int  status, glstat;

  getipv(N_STATUS, &status);
  getipv(N_GLSTAT, &glstat);

  return status == 1 || status == 3 || glstat == 6;
  // LP optimal || LP obj worse than cutoff || 
  // global search complete -- int found
}

bool OsiXprSolverInterface::isProvenPrimalInfeasible() const
{
  activateMe();

  int status;

  getipv(N_STATUS, &status);

  return status == 2; 		// LP infeasible
}

bool OsiXprSolverInterface::isProvenDualInfeasible() const
{
  activateMe();

  int status;

  getipv(N_STATUS, &status);

  return status == 5;		// LP Unbounded
}

bool OsiXprSolverInterface::isPrimalObjectiveLimitReached() const
{
  activateMe();

  return false;			// N/A in XOSL
}

bool OsiXprSolverInterface::isDualObjectiveLimitReached() const
{
  activateMe();

  int status;

  getipv(N_STATUS, &status);

  return status == 6;		// LP cut off in dual
}

bool OsiXprSolverInterface::isIterationLimitReached() const
{
  activateMe();

  int itrlim, itcnt;

  geticv(N_ITRLIM, &itrlim);
  getipv(N_ITCNT, &itcnt);

  return itcnt >= itrlim;
}

//#############################################################################
// WarmStart related methods
//#############################################################################

OsiWarmStart* OsiXprSolverInterface::getWarmStart() const
{
  activateMe();

  int pstat, retstat;

  getipv(N_PSTAT, &pstat);
  if ( pstat != 7 ) return NULL;

  OsiWarmStartBasis *ws = NULL;
  int numcols = getNumCols();
  int numrows = getNumRows();
  const double *lb = getColLower();
  double infty = getInfinity();
  int *rstatus = new int[numrows];
  int *cstatus = new int[numcols];

  retstat = getbasis(rstatus, cstatus);

  if ( retstat == 0 ) {
    int i;

    ws = new OsiWarmStartBasis;
    ws->setSize( numcols, numrows );
      
    for( i = 0;  i < numrows;  i++ ) {
      switch( rstatus[i] ) {
      case 0:
	ws->setArtifStatus(i, OsiWarmStartBasis::atLowerBound);
	break;
      case 1:
	ws->setArtifStatus(i, OsiWarmStartBasis::basic);
	break;
      case 2:
	ws->setArtifStatus(i, OsiWarmStartBasis::atUpperBound);
	break;
      default:  // unknown row status
	delete ws;
	ws = NULL;
	goto TERMINATE;
      }
    }

    for( i = 0;  i < numcols;  i++ ) {
      switch( cstatus[i] ) {
      case 0:
	if ( lb[i] <= -infty )
	  ws->setStructStatus(i, OsiWarmStartBasis::isFree);
	else
	  ws->setStructStatus(i, OsiWarmStartBasis::atLowerBound);
	break;
      case 1:
	ws->setStructStatus( i, OsiWarmStartBasis::basic );
	break;
      case 2:
	ws->setStructStatus( i, OsiWarmStartBasis::atUpperBound );
	break;
      default:  // unknown column status
	delete ws;
	ws = NULL;
	goto TERMINATE;
      }
    }
  }

 TERMINATE:
  delete[] cstatus;
  delete[] rstatus;

  return ws;
}

//-----------------------------------------------------------------------------

bool OsiXprSolverInterface::setWarmStart(const OsiWarmStart* warmstart)
{
  const OsiWarmStartBasis* ws = dynamic_cast<const OsiWarmStartBasis*>(warmstart);

  if ( !ws ) return false;

  activateMe();

  int numcols = ws->getNumStructural();
  int numrows = ws->getNumArtificial();

  if ( numcols != getNumCols() || numrows != getNumRows() )
    return false;

  bool retval;
  int retstat;
  int *cstatus = new int[numcols];
  int *rstatus = new int[numrows];
  int i;

  for ( i = 0;  i < numrows;  i++ ) {
    switch( ws->getArtifStatus(i) ) {
    case OsiWarmStartBasis::atLowerBound:
      rstatus[i] = 0;
      break;
    case OsiWarmStartBasis::basic:
      rstatus[i] = 1;
      break;
    case OsiWarmStartBasis::atUpperBound:
      rstatus[i] = 2;
      break;
    default:  // unknown row status
      retval = false;
      goto TERMINATE;
    }
  }

  for( i = 0;  i < numcols;  i++ ) {
    switch( ws->getStructStatus(i) ) {
    case OsiWarmStartBasis::atLowerBound: 
    case OsiWarmStartBasis::isFree:
      cstatus[i] = 0;
      break;
    case OsiWarmStartBasis::basic:
      cstatus[i] = 1;
      break;
    case OsiWarmStartBasis::atUpperBound:
      cstatus[i] = 2;
      break;
    default:  // unknown row status
      retval = false;
      goto TERMINATE;
    }
  }

  retstat = loadbasis(rstatus, cstatus);
  retval = (retstat == 0);

 TERMINATE:
  delete[] cstatus;
  delete[] rstatus;
  return retval;
  
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiXprSolverInterface::markHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::markHotStart();
}

void OsiXprSolverInterface::solveFromHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::solveFromHotStart();
}

void OsiXprSolverInterface::unmarkHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::unmarkHotStart();
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows and columns
//------------------------------------------------------------------
int
OsiXprSolverInterface::getNumCols() const
{
  activateMe();

  if ( !isDataLoaded() ) return 0;

  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"  int ncols;\n");
    fprintf(getLogFilePtr(),"  getipv(N_NCOL,&ncols);\n");
    fprintf(getLogFilePtr(),"}\n");
  }

  int     ncols;

  getipv(N_NCOL, &ncols);

  return ncols;
}
//-----------------------------------------------------------------------------
int
OsiXprSolverInterface::getNumRows() const
{
   activateMe();

   if ( !isDataLoaded() ) return 0;

   if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"  int nrows;\n");
    fprintf(getLogFilePtr(),"  getipv(N_NROW,&nrows);\n");
    fprintf(getLogFilePtr(),"}\n");
  }

   int     nrows;

   getipv(N_NROW, &nrows);

   return nrows;
}
//-----------------------------------------------------------------------------
int
OsiXprSolverInterface::getNumElements() const
{
   activateMe();

   if ( !isDataLoaded() ) return 0;

   if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"  int nels;\n");
    fprintf(getLogFilePtr(),"  getipv(N_NELEM,&nels);\n");
    fprintf(getLogFilePtr(),"}\n");
  }

   int     retVal;

   getipv(N_NELEM, &retVal);

   return retVal;
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  
const double *
OsiXprSolverInterface::getColLower() const
{
   if ( collower_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     ncols = getNumCols();
         
         if (getLogFilePtr()!=NULL) {
           fprintf(getLogFilePtr(),"{\n");
           fprintf(getLogFilePtr(),"  double * cl[%d];\n",ncols);
           fprintf(getLogFilePtr(),"  getbdl(cl,0,ncols-1);\n");
           fprintf(getLogFilePtr(),"}\n");
         }

	 collower_ = new double[ncols];
	 getbdl(collower_, 0, ncols - 1);
      }
   }

   return collower_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getColUpper() const
{
   if ( colupper_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     ncols = getNumCols();   
         
         if (getLogFilePtr()!=NULL) {
           fprintf(getLogFilePtr(),"{\n");
           fprintf(getLogFilePtr(),"  double * cu[%d];\n",ncols);
           fprintf(getLogFilePtr(),"  getbdu(cu,0,ncols-1);\n");
           fprintf(getLogFilePtr(),"}\n");
         }

	 colupper_ = new double[ncols];
	 getbdu(colupper_, 0, ncols - 1);
      }
   }

   return colupper_;
}
//-----------------------------------------------------------------------------
const char *
OsiXprSolverInterface::getRowSense() const
{
   if ( rowsense_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     nrows = getNumRows();
       
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  char rowsense[%d];\n",nrows);
          fprintf(getLogFilePtr(),"  getrowtype(rowsense, 0, %d);\n",nrows-1);
          fprintf(getLogFilePtr(),"}\n");
        }
        rowsense_ = new char[nrows];
        getrowtype(rowsense_, 0, nrows - 1);
      }
   }

   return rowsense_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getRightHandSide() const
{
   if ( rhs_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     nrows = getNumRows();  
         
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  double rhs[%d];\n",nrows);
          fprintf(getLogFilePtr(),"  getrhs(rhs, 0, %d);\n",nrows-1);
          fprintf(getLogFilePtr(),"}\n");
        }

	 rhs_ = new double[nrows];
	 getrhs(rhs_, 0, nrows - 1);

         // Make sure free rows have rhs of zero
         const char * rs = getRowSense();
         int nr = getNumRows();
         int i;
         for ( i = 0;  i < nr;  i++ ) {
           if ( rs[i] == 'N' ) rhs_[i]=0.0;
         }
      }
   }

   return rhs_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getRowRange() const
{
   if ( rowrange_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     nrows = getNumRows();

        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  double rowrange[%d];\n",nrows);
          fprintf(getLogFilePtr(),"  getrng(rowrange, 0, %d);\n",nrows-1);
          fprintf(getLogFilePtr(),"}\n");
        }
	 rowrange_ = new double[nrows];
	 getrng(rowrange_, 0, nrows - 1);

         // Make sure non-R rows have range of 0.0
         // XPRESS seems to set N and L rows to a range of Infinity
         const char * rs = getRowSense();
         int nr = getNumRows();
         int i;
         for ( i = 0;  i < nr;  i++ ) {
           if ( rs[i] != 'R' ) rowrange_[i] = 0.0;
         }
      }
   }

   return rowrange_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getRowLower() const
{
   if ( rowlower_ == NULL ) {
      int     nrows = getNumRows();
      const   char    *rowsense = getRowSense();
      const   double  *rhs      = getRightHandSide();
      const   double  *rowrange = getRowRange();

      if ( nrows > 0 ) {
	 rowlower_ = new double[nrows];

         double dum1;
	 for ( int i = 0;  i < nrows;  i++ ) {
           convertSenseToBound(rowsense[i], rhs[i], rowrange[i],
			       rowlower_[i], dum1);
	 }
      }
   }

   return rowlower_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getRowUpper() const
{
   if ( rowupper_ == NULL ) {
      int     nrows = getNumRows();
      const   char    *rowsense = getRowSense();
      const   double  *rhs      = getRightHandSide();
      const   double  *rowrange = getRowRange();

      if ( nrows > 0 ) {
	 rowupper_ = new double[nrows];

         double dum1;
	 for ( int i = 0;  i < nrows;  i++ ) {
           convertSenseToBound(rowsense[i], rhs[i], rowrange[i],
			       dum1, rowupper_[i]);
	 }
      }
   }

   return rowupper_;
}
//-----------------------------------------------------------------------------
const double *
OsiXprSolverInterface::getObjCoefficients() const
{
   if ( objcoeffs_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {
	 int     ncols = getNumCols();

        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  double objc[%d];\n",ncols);
          fprintf(getLogFilePtr(),"  getobj(objc,0,%d);\n",ncols);
          fprintf(getLogFilePtr(),"}\n");
        }
	 objcoeffs_ = new double[ncols];
	 getobj(objcoeffs_, 0, ncols - 1);
      }
   }

   return objcoeffs_;
}
//-----------------------------------------------------------------------------
double
OsiXprSolverInterface::getObjSense() const
{
   return objsense_;
}

//-----------------------------------------------------------------------------
// Return information on integrality
//-----------------------------------------------------------------------------

bool
OsiXprSolverInterface::isContinuous(int colNumber) const
{
   getVarTypes();

   //std::cerr <<"OsiXprSolverInterface::isContinuous " <<vartype_[colNumber] <<std::endl;
   if ( vartype_ == NULL ) return true;
   if ( vartype_[colNumber] == 'C' ) return true;
   return false;
}
//-----------------------------------------------------------------------------
#if 0
bool
OsiXprSolverInterface::isInteger( int colNumber ) const
{
   return !(isContinuous(colNumber));
}
//-----------------------------------------------------------------------------
bool
OsiXprSolverInterface::isBinary( int colNumber ) const
{
   const double *cu = colupper();
   const double *cl = collower();
  
   getVarTypes();

   if ( vartype_ == NULL ) return false;
   return (vartype_[colNumber] == 'I' || vartype_[colNumber] == 'B') && 
      (cu[colNumber] == 0.0 || cu[colNumber] == 1.0) && 
      (cl[colNumber] == 0.0 || cl[colNumber] == 1.0);
}
//-----------------------------------------------------------------------------
bool
OsiXprSolverInterface::isIntegerNonBinary( int colNumber ) const
{
   getVarTypes();

   if ( vartype_ == NULL ) return false;
   return (vartype_[colNumber] == 'I' || vartype_[colNumber] == 'B') &&
      !isBinary(colNumber);  
}
//-----------------------------------------------------------------------------
bool
OsiXprSolverInterface::isFreeBinary( int colNumber ) const
{
   const   double  *colupper = this->colupper();
   const   double  *collower = this->collower();

   getVarTypes();

   return isBinary(colNumber) && colupper[colNumber] != collower[colNumber];
}
#endif

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------

const OsiPackedMatrix *
OsiXprSolverInterface::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL ) {
    activateMe();

    if ( isDataLoaded() ) {

      int     nrows = getNumRows();
      int     ncols = getNumCols();
      int     nelems;

      getrows(NULL, NULL, NULL, 0, &nelems, 0, nrows - 1);
         
      if (getLogFilePtr()!=NULL) {
	fprintf(getLogFilePtr(),"{\n");
	fprintf(getLogFilePtr(),"   int start[%d];\n",nrows+1);
	fprintf(getLogFilePtr(),"   int length[%d];\n",nrows);
	fprintf(getLogFilePtr(),"   int index[%d];\n",nelems);
	fprintf(getLogFilePtr(),"   double element[%d];\n",nelems);
	fprintf(getLogFilePtr(),"   int nelems;\n");
	fprintf(getLogFilePtr(),"   getrows(NULL, NULL, NULL, 0, &nelems, 0, %d);\n",nrows-1);       
	fprintf(getLogFilePtr(),"   getrows(start, index, element, nelems, &nelems, 0, %d);\n",nrows-1);
	fprintf(getLogFilePtr(),"}\n");
      }
         
      int     *start   = new int   [nrows + 1];
      int     *length  = new int   [nrows];
      int     *index   = new int   [nelems];
      double  *element = new double[nelems];

      getrows(start, index, element, nelems, &nelems, 0, nrows - 1);

      std::adjacent_difference(start + 1, start + (nrows+1), length);
      
      matrixByRow_ = new OsiPackedMatrix();
      matrixByRow_->assignMatrix(false /* not column ordered */,
				 ncols, nrows, nelems,
				 element, index, start, length);
    } else {
      matrixByRow_ = new OsiPackedMatrix();
      matrixByRow_->reverseOrdering();
    }
  }

  return matrixByRow_;
} 

//-----------------------------------------------------------------------------
const OsiPackedMatrix *
OsiXprSolverInterface::getMatrixByCol() const
{
   if ( matrixByCol_ == NULL ) {
      matrixByCol_ = new OsiPackedMatrix(*getMatrixByRow());
      matrixByCol_->reverseOrdering();
   }

   return matrixByCol_;
}

//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double
OsiXprSolverInterface::getInfinity() const
{
   return DPLINF;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

const double *
OsiXprSolverInterface::getColSolution() const
{
   if ( colsol_ == NULL ) {
      activateMe();

      if ( isDataLoaded() ) {

        int nc = getNumCols();

        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  double colsol[%d];\n",nc);
          fprintf(getLogFilePtr(),"  solution(colsol, NULL, NULL, NULL);\n");
          fprintf(getLogFilePtr(),"}\n");
        }

	 colsol_ = new double[nc];
	 solution(colsol_, NULL, NULL, NULL);
      }
   }

   return colsol_;
}

//-----------------------------------------------------------------------------

const double *
OsiXprSolverInterface::getRowPrice() const
{
   if ( rowprice_ == NULL ) {
      activateMe();
      int nr = getNumRows();

      if ( isDataLoaded() ) {
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  double rowprice[%d];\n",nr);
          fprintf(getLogFilePtr(),"  solution(NULL, NULL, rowprice, NULL);\n");
          fprintf(getLogFilePtr(),"}\n");
        }
	 rowprice_ = new double[nr];
	 solution(NULL, NULL, rowprice_, NULL);
      }
   }
   return rowprice_;
}

//-----------------------------------------------------------------------------

const double * OsiXprSolverInterface::getReducedCost() const
{
  if ( colprice_ == NULL ) {
    activateMe();
    int nc = getNumCols();

    if ( isDataLoaded() ) {
      if (getLogFilePtr()!=NULL) {
	fprintf(getLogFilePtr(),"{\n");
	fprintf(getLogFilePtr(),"  double rowprice[%d];\n",nc);
	fprintf(getLogFilePtr(),"  solution(NULL, NULL, NULL, colprice_);\n");
	fprintf(getLogFilePtr(),"}\n");
      }
      colprice_ = new double[nc];
      solution(NULL, NULL, NULL, colprice_);
    }
  }
  return colprice_;
}

//-----------------------------------------------------------------------------

const double * OsiXprSolverInterface::getRowActivity() const
{
  if( rowact_ == NULL ) {
    activateMe();

    if ( isDataLoaded() ) {
      int nrows = getNumRows();
      const double *rhs = getRightHandSide();
      if( nrows > 0 ) {
	int status;

	getipv(N_PSTAT, &status);

	if ( status == 7 ) {
	  int i;

	  rowact_ = new double[nrows];
	  solution(NULL, rowact_, NULL, NULL);

	  for ( i = 0;  i < nrows;  i++ )
	    rowact_[i] = rhs[i] - rowact_[i];
	} else {
	  CoinFillN(rowact_, nrows, 0.0);
	}
      }
    }
  }
  return rowact_;
}

//-----------------------------------------------------------------------------

double
OsiXprSolverInterface::getObjValue() const
{
   activateMe();

   double  objvalue = 0;

   if ( isDataLoaded() ) {
     if (getLogFilePtr()!=NULL) {
       fprintf(getLogFilePtr(),"{\n");
       fprintf(getLogFilePtr(),"  double objvalue;\n");
       fprintf(getLogFilePtr(),"  getdpv(N_DOBJVL, &objvalue);\n");
       fprintf(getLogFilePtr(),"}\n");
     }
     getdpv(N_DOBJVL, &objvalue);
   }

   return objvalue;
}

//-----------------------------------------------------------------------------

int OsiXprSolverInterface::getIterationCount() const
{
  int itcnt;

  getipv(N_ITCNT, &itcnt);

  return itcnt;
}

//-----------------------------------------------------------------------------

std::vector<double*> OsiXprSolverInterface::getDualRays(int maxNumRays) const
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getDualRays",
		 "OsiXprSolverInterface");
  return std::vector<double*>();
}

//-----------------------------------------------------------------------------

std::vector<double*> OsiXprSolverInterface::getPrimalRays(int maxNumRays) const
{
#if 0
  // *FIXME* : Still need to expand column into full ncols-length vector

  const int nrows = getNumRows();
  int nrspar;
  getipv(N_NRSPAR, &nrspar);
  int junb;
  int retcode;

  retcode = getunb(&junb);

  if ( retcode != 0 ) 
    return std::vector<double *>(0, (double *) NULL);;

  double *ray = new double[nrows];


  if ( junb < nrows ) {		// it's a slack
    int i;

    for ( i = 0;  i < nrows;  i++ ) ray[i] = 0.0; 
    ray[junb] = 1.0; 
    retcode = ftran(ray);
  } else if ( junb >= nrows + nrspar && 
	      junb < nrows + nrspar + getNumCols() ){			
    				// it's a structural variable
    int *mstart = new int[nrows];
    int *mrowind = new int[nrows];
    double *dmatval = new double[nrows];
    int nelt;
    int jcol = junb - nrows - nrspar;

    retcode = getcols(mstart, mrowind, dmatval, nrows, &nelt, 
		      jcol, jcol); 
    /* Unpack into the zeroed array y */ 
    int i, ielt;

    for ( i = 0;  i < nrows;  i++ ) ray[i] = 0.0; 
    for ( ielt = 0;  ielt < nelt;  ielt++ ) 
      ray[mrowind[ielt]] = dmatval[ielt]; 
    retcode = ftran(ray);

    delete [] mstart;
    delete [] mrowind;
    delete [] dmatval;
  } else { 			// it's an error
    retcode = 1;
  }

  if ( retcode == 0 ) return std::vector<double *>(1, ray);
  else {
    delete ray;
    return std::vector<double *>(0, (double *) NULL); 
  }
#endif

  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getPrimalRays",
		 "OsiXprSolverInterface");
  return std::vector<double*>();
}

//-----------------------------------------------------------------------------

#if 0
OsiVectorInt
OsiXprSolverInterface::getFractionalIndices(const double etol) const
{
   OsiVectorInt retVal;
   int     numInts = numintvars();
   const   double  *sol = colsol();

   getVarTypes();
  
   OsiRelFltEq eq(etol);

   for ( int i = 0;  i < numInts;  i++ ) {
      double colSolElem = sol[ivarind_[i]];
      double distanceFromInteger = colSolElem - floor(colSolElem + 0.5);

      if ( !eq( distanceFromInteger, 0.0 ) )
	 retVal.push_back(ivarind_[i]);
   }

   return retVal;
}
#endif

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

void
OsiXprSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
   activateMe();

   if ( isDataLoaded() ) {
      chgobj(1, &elementIndex, &elementValue);
      freeCachedResults();
   }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setColLower( int elementIndex, double elementValue )
{
   activateMe();

   if ( isDataLoaded() ) {
      char boundType = 'L';

      getVarTypes();
      
      if (getLogFilePtr()!=NULL) {
        fprintf(getLogFilePtr(),"chgbds(1, %d, %c, %f );\n",
          elementIndex,boundType,elementValue);
      }

      chgbds(1, &elementIndex, &boundType, &elementValue);
      if ( vartype_[elementIndex] == 'B' && 
	   (elementValue != 0.0 && elementValue != 1.0) ) {
        char elementType = 'I';
        
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"chgcoltype(1, %d, %c );\n",
            elementIndex,elementType);
        }
	chgcoltype(1, &elementIndex, &elementType);
      }
      freeCachedResults();
      //    delete [] collower_;
      //    collower_ = NULL;
   }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setColUpper( int elementIndex, double elementValue )
{
   activateMe();

   if ( isDataLoaded() ) {
      char boundType = 'U';

      getVarTypes();      
      if (getLogFilePtr()!=NULL) {
        fprintf(getLogFilePtr(),"chgbds(1, %d, %c, %f );\n",
		elementIndex,boundType,elementValue);
      }
      chgbds(1, &elementIndex, &boundType, &elementValue);
      if ( vartype_[elementIndex] == 'B' && 
	   (elementValue != 0.0 && elementValue != 1.0) ) {
	 char elementType = 'I';  
         
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"chgcoltype(1, %d, %c );\n",
		  elementIndex,elementType);
        }

	 chgcoltype(1, &elementIndex, &elementType);
      }
      freeCachedResults();
      //    delete [] colupper_;
      //    colupper_ = NULL;
   } 
}

//-----------------------------------------------------------------------------

void OsiXprSolverInterface::setColBounds(const int elementIndex, double lower, double upper )
{
   if ( isDataLoaded() ) {
     char qbtype[2] = { 'L', 'U' };
     int mindex[2];
     double bnd[2];

     mindex[0] = elementIndex;
     mindex[1] = elementIndex;
     bnd[0] = lower;
     bnd[1] = upper;

     chgbds(2, mindex, qbtype, bnd);
     if ( vartype_[mindex[0]] == 'B' && 
	  !((lower == 0.0 && upper == 0.0) ||
	    (lower == 1.0 && upper == 1.0) ||
	    (lower == 0.0 && upper == 1.0)) ) {
       char elementType = 'I';  
         
       chgcoltype(1, &mindex[0], &elementType);
     }
     freeCachedResults();
     //    delete [] colupper_;
     //    colupper_ = NULL;
   }
}

//-----------------------------------------------------------------------------

void OsiXprSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  OsiSolverInterface::setColSetBounds(indexFirst, indexLast, boundList);
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowLower( int elementIndex, double elementValue )
{
  // activateMe();

  double rhs   = getRightHandSide()[elementIndex];
  double range = getRowRange()[elementIndex];
  char   sense = getRowSense()[elementIndex];
  double lower, upper;

  convertSenseToBound(sense, rhs, range, lower, upper);
  if( lower != elementValue ) {
    convertBoundToSense(elementValue, upper, sense, rhs, range);
    setRowType(elementIndex, sense, rhs, range);
    // freeCachedResults(); --- invoked in setRowType()
  }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowUpper( int elementIndex, double elementValue )
{
  // activateMe();

  double rhs   = getRightHandSide()[elementIndex];
  double range = getRowRange()[elementIndex];
  char   sense = getRowSense()[elementIndex];
  double lower, upper;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( upper != elementValue ) {
    convertBoundToSense(lower, elementValue, sense, rhs, range);
    setRowType(elementIndex, sense, rhs, range);
    // freeCachedResults(); --- invoked in setRowType()
  }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
  double rhs, range;
  char sense;
  
  convertBoundToSense( lower, upper, sense, rhs, range );
  setRowType( elementIndex, sense, rhs, range );
  // freeCachedRowRim(); --- invoked in setRowType()
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowType(int index, char sense, double rightHandSide,
				  double range)
{
  activateMe();

  if ( isDataLoaded() ) {
    int mindex[1] = {index};
    char qrtype[1] = {sense}; 
    double rhs[1] = {rightHandSide};
    double rng[1] = {range};

    chgrowtype(1, mindex, qrtype);
    chgrhs(1, mindex, rhs);
    chgrng(1, mindex, rng);

    freeCachedResults();
  }
}

//-----------------------------------------------------------------------------

void OsiXprSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  OsiSolverInterface::setRowSetBounds(indexFirst, indexLast, boundList);
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  OsiSolverInterface::setRowSetTypes(indexFirst, indexLast, senseList, rhsList, rangeList);
}

//#############################################################################
void 
OsiXprSolverInterface::setContinuous(int index) 
{
  activateMe();

  if ( isDataLoaded() ) {
    int pstat;

    getipv(N_PSTAT, &pstat);

    if ( pstat & 6 == 0 ) { 		// not presolved
      char qctype = 'C';

      chgcoltype(1, &index, &qctype);
      freeCachedResults();
    }
  }
}

void 
OsiXprSolverInterface::setInteger(int index) 
{
  activateMe();

  if ( isDataLoaded() ) {
    int pstat;

    getipv(N_PSTAT, &pstat);

    if ( pstat & 6 == 0 ) { 		// not presolved
      char qctype;

      if ( getColLower()[index] == 0.0 && 
	   getColUpper()[index] == 1.0 ) 
	qctype = 'B';
      else
	qctype = 'I';

      chgcoltype(1, &index, &qctype);
      freeCachedResults();
    }
  }
}

void 
OsiXprSolverInterface::setContinuous(const int* indices, int len) 
{
  activateMe();

  if ( isDataLoaded() ) {
    int pstat;

    getipv(N_PSTAT, &pstat);

    if ( pstat & 6 == 0 ) { 		// not presolved
      char *qctype = new char[len];

      CoinFillN(qctype, len, 'C');
      chgcoltype(1, const_cast<int *>(indices), qctype);
      freeCachedResults();
    }
  }
}

void 
OsiXprSolverInterface::setInteger(const int* indices, int len) 
{
  activateMe();

  if ( isDataLoaded() ) {
    int pstat;

    getipv(N_PSTAT, &pstat);

    if ( pstat & 6 == 0 ) { 		// not presolved
      char *qctype = new char[len];
      const double* clb = getColLower();
      const double* cub = getColUpper();

      for ( int i = 0;  i < len;  i++ ) {
	if ( clb[indices[i]] == 0.0 && cub[indices[i]] == 1.0 )
	  qctype[i] = 'B';
	else 
	  qctype[i] = 'I';
      }

      chgcoltype(1, const_cast<int *>(indices), qctype);
      freeCachedResults();
    }
  }
}

//#############################################################################

void
OsiXprSolverInterface::setObjSense(double s) 
{
   objsense_ = s;
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setColSolution(const double *colsol)
{
   activateMe();

   freeSolution();

   colsol_ = new double[getNumCols()];

   for ( int i = 0;  i < getNumCols();  i++ )
      colsol_[i] = colsol[i];
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::setRowPrice(const double *rowprice)
{
   activateMe();

   freeSolution();

   rowprice_ = new double[getNumRows()];

   for ( int i = 0;  i < getNumRows();  i++ )
      rowprice_[i] = rowprice[i];
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################

void 
OsiXprSolverInterface::addCol(const OsiPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  activateMe();

  if ( isDataLoaded() ) {
    freeCachedResults();

    int mstart = 0;

    addcols(1, vec.getNumElements(), const_cast<double*>(&obj),
	    &mstart,
	    const_cast<int*>(vec.getIndices()),
	    const_cast<double*>(vec.getElements()),
	    const_cast<double*>(&collb),
	    const_cast<double*>(&colub));
  }
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::addCols(const int numcols,
			       const OsiPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  // activateMe();
  // freeCachedResults();

  for( int i = 0;  i < numcols;  i++ )
    addCol( *(cols[i]), collb[i], colub[i], obj[i] );
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::deleteCols(const int num, const int *columnIndices)
{
  activateMe();
  freeCachedResults();
  delcols(num, const_cast<int *>(columnIndices));
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  // activateMe(); -- will be invoked
  // freeCachedResults(); -- will be invoked

  char sense;
  double rhs, range;

  convertBoundToSense(rowlb, rowub, sense, rhs, range);
  addRow(vec, sense, rhs, range);
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  activateMe();
  freeCachedResults();

  int mstart[1] = {0};

  addrows(1, vec.getNumElements(), 
	  const_cast<char *>(&rowsen), const_cast<double *>(&rowrhs), 
	  const_cast<double *>(&rowrng), mstart,
	  const_cast<int *>(vec.getIndices()),
	  const_cast<double *>(vec.getElements()));
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "addRows",
		 "OsiXprSolverInterface");
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "addRows",
		 "OsiXprSolverInterface");
}
//-----------------------------------------------------------------------------
void 
OsiXprSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  activateMe();
  freeCachedResults();

  delrows(num, const_cast<int *>(rowIndices));
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiXprSolverInterface::loadProblem(const OsiPackedMatrix& matrix,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  const double inf = getInfinity();
  
  char   * rowSense = new char  [matrix.getNumRows()];
  double * rowRhs   = new double[matrix.getNumRows()];
  double * rowRange = new double[matrix.getNumRows()];
  
  int i;
  for ( i = matrix.getNumRows() - 1; i >= 0; --i) {
    
    double rlb;
    if ( rowlb!=NULL )
      rlb = rowlb[i];
    else
      rlb = -inf;
    
     double rub;
     if ( rowub!=NULL )
       rub = rowub[i];
     else
       rub = inf;
     
     convertBoundToSense(rlb,rub,rowSense[i],rowRhs[i],rowRange[i]);
#if 0
     if ( rlb==rub ) {
       rowSense[i]='E';
       rowRhs[i]  =rlb;
       rowRange[i]=0.0;
       continue;
     }
     if ( rlb<=-inf && rub>=inf ) {
       rowSense[i]='N';
       rowRhs[i]  =inf;
       rowRange[i]=0.0;
       continue;
     }
     if ( rlb<=-inf && !(rub>=inf) ) {
       rowSense[i]='L';
       rowRhs[i]  =rub;
       rowRange[i]=0.0;
       continue;
     }
     if ( !(rlb<=-inf) && rub>=inf ) {
       rowSense[i]='G';
       rowRhs[i]  =rlb;
       rowRange[i]=0.0;
       continue;
     }
     if ( !(rlb<=-inf) && !(rub>=inf) ) {
       rowSense[i]='R';
       rowRhs[i]  =rub;
       rowRange[i]=rub-rlb;
       continue;
     }
#endif
   }
  
   loadProblem(matrix, collb, colub, obj, rowSense, rowRhs, rowRange ); 
   
   delete [] rowSense;
   delete [] rowRhs;
   delete [] rowRange;
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::assignProblem(OsiPackedMatrix*& matrix,
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

// #define OSIXPR_ADD_OBJ_ROW

void
OsiXprSolverInterface::loadProblem(const OsiPackedMatrix& matrix,
                                   const double* collb, const double* colub,
                                   const double* obj,
                                   const char* rowsen, const double* rowrhs,   
                                   const double* rowrng)
{
  assert( rowsen != NULL );
  assert( rowrhs != NULL );

  activateMe();
  freeCachedResults();
  int i;
 
  // Set column values to defaults if NULL pointer passed
  int nc=matrix.getNumCols();
  double * clb;  
  double * cub;
  double * ob;
  if ( collb!=NULL ) {
    clb=const_cast<double*>(collb);
  }
  else {
    clb = new double[nc];
    for( i=0; i<nc; i++ ) clb[i]=0.0;
  }
  if ( colub!=NULL ) 
    cub=const_cast<double*>(colub);
  else {
    cub = new double[nc];
    for( i=0; i<nc; i++ ) cub[i]=DPLINF;
  }
  if ( obj!=NULL ) 
    ob=const_cast<double*>(obj);
  else {
    ob = new double[nc];
    for( i=0; i<nc; i++ ) ob[i]=0.0;
  }
  
  bool freeMatrixRequired = false;
  OsiPackedMatrix * m = NULL;
  if ( !matrix.isColOrdered() ) {
    m = new OsiPackedMatrix();
    m->reverseOrderedCopyOf(matrix);
    freeMatrixRequired = true;
  } else {
    m = const_cast<OsiPackedMatrix *>(&matrix);
  }
  
  // Generate a problem name
  char probName[256];
  sprintf(probName, "Prob%i", osiSerial_);

  
  nc = m->getNumCols();
  int nr =  m->getNumRows();
  
  if ( getLogFilePtr()!=NULL ) {   
    fprintf(getLogFilePtr(),"{\n"); 

    fprintf(getLogFilePtr(),"  char rowsen[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowsen[%d]='%c';\n",i,rowsen[i]);

    fprintf(getLogFilePtr(),"  double rowrhs[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowrhs[%d]=%f;\n",i,rowrhs[i]);
    
    fprintf(getLogFilePtr(),"  double rowrng[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowrng[%d]=%f;\n",i,rowrng[i]);

    fprintf(getLogFilePtr(),"  double ob[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  ob[%d]=%f;\n",i,ob[i]);

    fprintf(getLogFilePtr(),"  double clb[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  clb[%d]=%f;\n",i,clb[i]);

    fprintf(getLogFilePtr(),"  double cub[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  cub[%d]=%f;\n",i,cub[i]);

    fprintf(getLogFilePtr(),"  int vectorStarts[%d];\n",nc+1);
    for ( i=0; i<=nc; i++ )
      fprintf(getLogFilePtr(),"  vectorStarts[%d]=%d;\n",i,m->getVectorStarts()[i]);

    fprintf(getLogFilePtr(),"  int vectorLengths[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  vectorLengths[%d]=%d;\n",i,m->getVectorLengths()[i]);
    
    fprintf(getLogFilePtr(),"  int indices[%d];\n",m->getVectorStarts()[nc]);
    for ( i=0; i<m->getVectorStarts()[nc]; i++ )
      fprintf(getLogFilePtr(),"  indices[%d]=%d;\n",i,m->getIndices()[i]);

    fprintf(getLogFilePtr(),"  double elements[%d];\n",m->getVectorStarts()[nc]);
    for ( i=0; i<m->getVectorStarts()[nc]; i++ )
      fprintf(getLogFilePtr(),"  elements[%d]=%f;\n",i,m->getElements()[i]);

    fprintf(getLogFilePtr(),
            "  int iret = loadprob(\"%s\",\n"
            "                      %d,\n"
            "                      %d,\n"
            "                      rowsen,\n"
            "                      rowrhs,\n"
            "                      rowrng,\n"
            "                      ob,\n"
            "                      vectorStarts,\n"
            "                      vectorLengths,\n"
            "                      indices,\n"
            "                      elements,\n"
            "                      clb,\n"
            "                      cub );\n",probName,nc,nr );    
    fprintf(getLogFilePtr(),"}\n");
  }
  // Need to cast away const'ness
  int iret = loadprob(probName,
    nc,
    nr,
    rowsen,
    const_cast<double*>(rowrhs),
    const_cast<double*>(rowrng),
    ob,
    const_cast<int*>(m->getVectorStarts()),
    const_cast<int*>(m->getVectorLengths()),
    const_cast<int*>(m->getIndices()),
    const_cast<double*>(m->getElements()),
    clb,
    cub );
  setStrParam(OsiProbName,probName);
  
  if ( iret != 0 )
    getipv(N_ERRNO, &iret);
  assert( iret == 0 );
  
   
  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"   char pname[256];\n");
    fprintf(getLogFilePtr(),"   getprob(pname);\n");
    fprintf(getLogFilePtr(),"}\n");
  }

  char pname[256];      // Problem names can be 200 chars in XPRESS 12
  getprob(pname);
  xprProbname_ = pname;
  
  if ( collb==NULL ) delete[] clb;
  if ( colub==NULL ) delete[] cub;
  if ( obj  ==NULL ) delete[] ob;
  
  if (freeMatrixRequired) {
    delete m;
  }
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::assignProblem(OsiPackedMatrix*& matrix,
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
OsiXprSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub )
{
  const double inf = getInfinity();
  
  char   * rowSense = new char  [numrows];
  double * rowRhs   = new double[numrows];
  double * rowRange = new double[numrows];
  
  for ( int i = numrows - 1; i >= 0; --i ) {
    const double lower = rowlb ? rowlb[i] : -inf;
    const double upper = rowub ? rowub[i] : inf;
    convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], rowRange[i] );
  }

  loadProblem(numcols, numrows, start, index, value, collb, colub, obj,
	      rowSense, rowRhs, rowRange);
  delete [] rowSense;
  delete [] rowRhs;
  delete [] rowRange;

}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const char* rowsen, const double* rowrhs,
				   const double* rowrng )
{
  assert( rowsen != NULL );
  assert( rowrhs != NULL );

  activateMe();
  freeCachedResults();
  int i;
 
  // Set column values to defaults if NULL pointer passed
  int nc = numcols;
  int nr = numrows;
  int * len = new int[nc];
  double * clb;  
  double * cub;
  double * ob;

  std::adjacent_difference(start, start + (nc+1), len);
  
  if ( collb!=NULL ) {
    clb=const_cast<double*>(collb);
  }
  else {
    clb = new double[nc];
    for( i=0; i<nc; i++ ) clb[i]=0.0;
  }
  if ( colub!=NULL ) 
    cub=const_cast<double*>(colub);
  else {
    cub = new double[nc];
    for( i=0; i<nc; i++ ) cub[i]=DPLINF;
  }
  if ( obj!=NULL ) 
    ob=const_cast<double*>(obj);
  else {
    ob = new double[nc];
    for( i=0; i<nc; i++ ) ob[i]=0.0;
  }

  // Generate a problem name
  char probName[256];
  sprintf(probName, "Prob%i", osiSerial_);

  if ( getLogFilePtr()!=NULL ) {   
    fprintf(getLogFilePtr(),"{\n"); 

    fprintf(getLogFilePtr(),"  char rowsen[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowsen[%d]='%c';\n",i,rowsen[i]);

    fprintf(getLogFilePtr(),"  double rowrhs[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowrhs[%d]=%f;\n",i,rowrhs[i]);
    
    fprintf(getLogFilePtr(),"  double rowrng[%d];\n",nr);
    for ( i=0; i<nr; i++ )
      fprintf(getLogFilePtr(),"  rowrng[%d]=%f;\n",i,rowrng[i]);

    fprintf(getLogFilePtr(),"  double ob[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  ob[%d]=%f;\n",i,ob[i]);

    fprintf(getLogFilePtr(),"  double clb[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  clb[%d]=%f;\n",i,clb[i]);

    fprintf(getLogFilePtr(),"  double cub[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  cub[%d]=%f;\n",i,cub[i]);

    fprintf(getLogFilePtr(),"  int vectorStarts[%d];\n",nc+1);
    for ( i=0; i<=nc; i++ )
      fprintf(getLogFilePtr(),"  vectorStarts[%d]=%d;\n",i,start[i]);

    fprintf(getLogFilePtr(),"  int vectorLengths[%d];\n",nc);
    for ( i=0; i<nc; i++ )
      fprintf(getLogFilePtr(),"  vectorLengths[%d]=%d;\n",i,len[i]);
    
    fprintf(getLogFilePtr(),"  int indices[%d];\n",start[nc]);
    for ( i=0; i<start[nc]; i++ )
      fprintf(getLogFilePtr(),"  indices[%d]=%d;\n",i,index[i]);

    fprintf(getLogFilePtr(),"  double elements[%d];\n",start[nc]);
    for ( i=0; i<start[nc]; i++ )
      fprintf(getLogFilePtr(),"  elements[%d]=%f;\n",i,value[i]);

    fprintf(getLogFilePtr(),
            "  int iret = loadprob(\"%s\",\n"
            "                      %d,\n"
            "                      %d,\n"
            "                      rowsen,\n"
            "                      rowrhs,\n"
            "                      rowrng,\n"
            "                      ob,\n"
            "                      vectorStarts,\n"
            "                      vectorLengths,\n"
            "                      indices,\n"
            "                      elements,\n"
            "                      clb,\n"
            "                      cub );\n",probName,nc,nr );    
    fprintf(getLogFilePtr(),"}\n");
  }
  // Need to cast away const'ness
  int iret = loadprob(probName,
    nc,
    nr,
    rowsen,
    const_cast<double*>(rowrhs),
    const_cast<double*>(rowrng),
    ob,
    const_cast<int*>(start),
    const_cast<int*>(len),
    const_cast<int*>(index),
    const_cast<double*>(value),
    clb,
    cub );
  setStrParam(OsiProbName,probName);
  
  if ( iret != 0 )
    getipv(N_ERRNO, &iret);
  assert( iret == 0 );
  
   
  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"   char pname[256];\n");
    fprintf(getLogFilePtr(),"   getprob(pname);\n");
    fprintf(getLogFilePtr(),"}\n");
  }

  char pname[256];      // Problem names can be 200 chars in XPRESS 12
  getprob(pname);
  xprProbname_ = pname;
  
  if ( collb == NULL ) delete[] clb;
  if ( colub == NULL ) delete[] cub;
  if ( obj   == NULL ) delete[] ob;
  delete[] len;
}



//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------

int
OsiXprSolverInterface::readMps(const char *filename, const char *extension) 
{
  activateMe();
   
#if 0
  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"{\n");
    fprintf(getLogFilePtr(),"  input(\"%s\");\n",filename);
    fprintf(getLogFilePtr(),"  int namlen;\n");
    fprintf(getLogFilePtr(),"  geticv(N_NAMSIZ,&namlen);\n");
    fprintf(getLogFilePtr(),"  namlen *= 8;\n");
  }
 
  // Read Mps file.
  // XPRESS generates its own file extensions, so we ignore any supplied.
  int iret = input(filename);
  if ( iret != 0 )  
    getipv(N_ERRNO, &iret);
  assert( iret == 0 );

  // Get length of Mps names
  int namlen;
  getipv(N_NAMSIZ, &namlen);
  namlen *= 8;

  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"  char objRowName[%d],rowName[%d];\n",namlen,namlen);
    fprintf(getLogFilePtr(),"  getccv(N_OBJNAM, objRowName);\n");
    fprintf(getLogFilePtr(),"  int nr;\n");
    fprintf(getLogFilePtr(),"  seticv(N_NROW, &nr);");
  }

  // Allocate space to hold row names.
  char * objRowName = new char[namlen+1];
  char * rowName    = new char[namlen+1];

  // Get name of objective row.
  // If "" returned, then first N row is objective row
  getccv(N_OBJNAM,objRowName);

  // Get number of rows
  int nr;
  getipv(N_NROW, &nr);
    
  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"  char rs[%d];\n",nr);
    fprintf(getLogFilePtr(),"  getrowtype(rs, 0, %d);\n",nr-1);
  }

  // Get row sense.
  char * rs = new char[nr];
  getrowtype(rs, 0, nr - 1);

  // Loop once for each row looking for objective row
  int i;
  for ( i=0; i<nr; i++ ) {
    
    // Objective row must be an N row
    if ( rs[i]=='N' ) {      
      
      if (getLogFilePtr()!=NULL) {
        fprintf(getLogFilePtr(),"  getnames(1,rowName,%d,%d);\n",i,i);
      }
      
      // Get name of this row
      getnames(1,rowName,i,i);
      
      // Is this the objective row?
      if( strcmp(rowName,objRowName)==0 ||
        strcmp("",objRowName) == 0 ) {
        
        if (getLogFilePtr()!=NULL) {
          fprintf(getLogFilePtr(),"  int rowToDelete[1];\n");
          fprintf(getLogFilePtr(),"  rowToDelete[0]=%d;\n",i);
          fprintf(getLogFilePtr(),"  delrows(1,rowToDelete);\n");
        }
        
        // found objective row. now delete it
        int rowToDelete[1];
        rowToDelete[0] = i;
        delrows(1,rowToDelete);
        break;
      }
    }
  }

  delete [] rs;
  delete [] rowName;
  delete [] objRowName;
  
  if (getLogFilePtr()!=NULL) {
    fprintf(getLogFilePtr(),"  char pname[256];\n");
    fprintf(getLogFilePtr(),"  getprob(pname);\n");
    fprintf(getLogFilePtr(),"}\n");
  }
  
  char pname[256];      // Problem names can be 200 chars in XPRESS 12
  getprob(pname);
  xprProbname_ = pname;
  return 0;
#else
  int retVal = OsiSolverInterface::readMps(filename,extension);
  getStrParam(OsiProbName,xprProbname_);
  return retVal;
#endif
}


//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------
void OsiXprSolverInterface::writeMps(const char *filename,
				     const char *extension) const
{
  activateMe();

  // Note: XPRESS insists on ignoring the extension and 
  // adding ".mat" instead.  Forewarned is forearmed.
  // Note: Passing an empty filename produces a file named after 
  // the problem.

  output(filename, "");
}

//#############################################################################
// Static data and methods
//#############################################################################

void
OsiXprSolverInterface::incrementInstanceCounter()
{
  if ( numInstances_ == 0 ) {          
    if ( getLogFilePtr()!=NULL ) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"  setoptlog(\"xpress.log\");\n");
      fprintf(getLogFilePtr(),"  initlz(NULL,0);\n");
      fprintf(getLogFilePtr(),"}\n");
    }
    const char *logfile = "xpress.log";
    setoptlog(logfile);
    int iret = initlz(NULL, 0);
    
    // Student Version returns 32
    if ( iret != 32 ) {
      if ( iret != 0 ) getipv(N_ERRNO, &iret);
      assert(iret == 0);
    }
  }

   numInstances_++;
   osiSerial_++;

   //  cout << "Instances = " << numInstances_ << "; Serial = " << osiSerial_ << endl;
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::decrementInstanceCounter()
{
  assert( numInstances_ != 0 );
  numInstances_--;
  
  //  cout << "Instances = " << numInstances_ << endl;
  
  if ( numInstances_ == 0 ) {
    
    if ( getLogFilePtr()!=NULL ) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"  freexo();\n");
      fprintf(getLogFilePtr(),"}\n");
    }

    freexo();
  }
}

//-----------------------------------------------------------------------------

unsigned int
OsiXprSolverInterface::getNumInstances()
{
   return numInstances_;
}

//-----------------------------------------------------------------------------

// Return XPRESS-MP Version number
int OsiXprSolverInterface::version()
{
  int retVal;
  geticv( N_VERNO, &retVal );
  return retVal;
}

//-----------------------------------------------------------------------------

FILE * OsiXprSolverInterface::getLogFilePtr()
{
  if (logFilePtr_!=NULL){fclose(logFilePtr_);logFilePtr_=NULL;}
  if (logFileName_!=NULL) logFilePtr_ = fopen(logFileName_,"a+");
  return logFilePtr_;
}

//-----------------------------------------------------------------------------

void OsiXprSolverInterface::setLogFileName( const char * filename )
{
  logFileName_ = filename;
}

//-----------------------------------------------------------------------------

unsigned int OsiXprSolverInterface::numInstances_ = 0;

unsigned int OsiXprSolverInterface::osiSerial_ = 0;  

FILE * OsiXprSolverInterface::logFilePtr_ = NULL;  

const char * OsiXprSolverInterface::logFileName_ = NULL;

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiXprSolverInterface::OsiXprSolverInterface (int newrows, int newnz) :
OsiSolverInterface(),
xprSaved_(false),
xprMatrixId_(-1),
matrixByRow_(NULL),
matrixByCol_(NULL),
colupper_(NULL),
collower_(NULL),
rowupper_(NULL),
rowlower_(NULL),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
objcoeffs_(NULL),
objsense_(1),
colsol_(NULL),
rowsol_(NULL),
rowact_(NULL),
rowprice_(NULL),
colprice_(NULL),
ivarind_(NULL),
ivartype_(NULL),
vartype_(NULL)
{
  incrementInstanceCounter();
  
  xprProbname_ = "";
  
  // newrows and newnz specify room to leave in the solver's matrix
  // structure for cuts.  Note that this is a *global* parameter.
  // The value in effect when the problem is loaded pertains.
  
  if ( newrows > 0 && newnz > 0 ) {         
    if ( getLogFilePtr()!=NULL ) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"  seticv(N_NRXTRA, %d);\n",newrows);
      fprintf(getLogFilePtr(),"  seticv(N_NMXTRA, %d);\n",newnz);
      fprintf(getLogFilePtr(),"}\n");
    }
    seticv(N_NRXTRA, newrows);
    seticv(N_NMXTRA, newnz);
  }
}

//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiSolverInterface *
OsiXprSolverInterface::clone(bool copyData) const
{  
   return (new OsiXprSolverInterface(*this));
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiXprSolverInterface::
OsiXprSolverInterface (const OsiXprSolverInterface & source) :
   OsiSolverInterface(source),
   xprSaved_(false),
   xprMatrixId_(-1),
   matrixByRow_(NULL),
   matrixByCol_(NULL),
   colupper_(NULL),
   collower_(NULL),
   rowupper_(NULL),
   rowlower_(NULL),
   rowsense_(NULL),
   rhs_(NULL),
   rowrange_(NULL),
   objcoeffs_(NULL),
   objsense_(1),
   colsol_(NULL),
   rowsol_(NULL),
   rowact_(NULL),
   rowprice_(NULL),
   colprice_(NULL),
   ivarind_(NULL),
   ivartype_(NULL),
   vartype_(NULL)
{
   incrementInstanceCounter();
  
   xprProbname_ = "";

   gutsOfCopy(source);

   // Other values remain NULL until requested
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiXprSolverInterface::~OsiXprSolverInterface ()
{
   if ( xprSaved_ ) {
      //    cout << "Problem " << xprProbname_ << " deleted from matrix " << xprMatrixId_ << "." << endl; 
      //    *** Temporarily no matrix deletes until XPRESS bug resolved.
      //    int iret = delmat(xprMatrixId_);
      //    if ( iret != 0 ) { 
      //      getipv(N_ERRNO, &iret);
      //      cout << "Deletion reported error " << iret << endl;
      //    }
      //    assert( iret == 0 );
   } else if ( xprCurrentProblem_ == this ) xprCurrentProblem_ = NULL;

   gutsOfDestructor();

   decrementInstanceCounter();
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiXprSolverInterface &
OsiXprSolverInterface::operator=(const OsiXprSolverInterface& rhs)
{
   if ( this != &rhs ) {    
      OsiSolverInterface::operator=(rhs);
      osiSerial_++;       // even though numInstances_ doesn't change
      gutsOfDestructor();
      gutsOfCopy(rhs);
   }
   return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void
OsiXprSolverInterface::applyColCut( const OsiColCut & cc )
{
   activateMe();

   const double  *collower = getColLower();
   const double  *colupper = getColUpper();
   const OsiPackedVector & lbs = cc.lbs();
   const OsiPackedVector & ubs = cc.ubs();

   int     *index = new int   [lbs.getNumElements() + ubs.getNumElements()];
   char    *btype = new char  [lbs.getNumElements() + ubs.getNumElements()];
   double  *value = new double[lbs.getNumElements() + ubs.getNumElements()];

   int     i, nbds;

   for ( i = nbds = 0;  i < lbs.getNumElements();  i++, nbds++ ) {
      index[nbds] = lbs.getIndices()[i];
      btype[nbds] = 'L';
      value[nbds] = (collower[index[nbds]] > lbs.getElements()[i]) ?
	 collower[index[nbds]] : lbs.getElements()[i];
   }

   for ( i = 0;  i < ubs.getNumElements();  i++, nbds++ ) {
      index[nbds] = ubs.getIndices()[i];
      btype[nbds] = 'U';
      value[nbds] = (colupper[index[nbds]] < ubs.getElements()[i]) ?
	 colupper[index[nbds]] : ubs.getElements()[i];
   }
   
    if (getLogFilePtr()!=NULL) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"   int index[%d];\n",nbds);
      fprintf(getLogFilePtr(),"   char btype[%d];\n",nbds);
      fprintf(getLogFilePtr(),"   double value[%d];\n",nbds);
      for ( i=0; i<nbds; i++ ) {
        fprintf(getLogFilePtr(),"   index[%d]=%d;\n",i,index[i]);
        fprintf(getLogFilePtr(),"   btype[%d]='%c';\n",i,btype[i]);
        fprintf(getLogFilePtr(),"   value[%d]=%f;\n",i,value[i]);
      }
      fprintf(getLogFilePtr(),"    chgbds(%d, index, btype, value);\n",nbds);
      fprintf(getLogFilePtr(),"}\n");
    }

   chgbds(nbds, index, btype, value);

   delete [] index;
   delete [] btype;
   delete [] value;

   freeCachedResults();
   //  delete [] colupper_;      colupper_ = NULL;
   //  delete [] collower_;      collower_ = NULL;
}

//-----------------------------------------------------------------------------

void
OsiXprSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
   activateMe();

   const   OsiPackedVector & row=rowCut.row();
   int     start[2] = {0, row.getNumElements()};
   char    sense = rowCut.sense();
   double  rhs   = rowCut.rhs();

   // XPRESS-MP header file xpresso.h is defining range
   // This is causing the rowCut range method name to
   // be incorretly changed.
   // Play some preprocessor games to get around.
#define rangeTemp range
#undef range
   double  r = rowCut.range();
#define range rangeTemp

    if (getLogFilePtr()!=NULL) {
      int i;
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"  char sense = '%c';\n",sense);
      fprintf(getLogFilePtr(),"  double rhs = '%f';\n",rhs);
      fprintf(getLogFilePtr(),"  double r = '%f';\n",r);
      fprintf(getLogFilePtr(),"  int start[2] = {0,%d};\n",start[1]);
      fprintf(getLogFilePtr(),"  int indices[%d];\n",row.getNumElements());
      fprintf(getLogFilePtr(),"  double elements[%d];\n",row.getNumElements());
      for ( i=0; i<row.getNumElements(); i++ ) {
        fprintf(getLogFilePtr(),"  indices[%d]=%d;\n",i,row.getIndices()[i]);
        fprintf(getLogFilePtr(),"  elements[%d]=%f;\n",i,row.getElements()[i]);
      }
      fprintf(getLogFilePtr(),
        "  int rc = addrows(1, %d, &sense, &rhs, &r, start,indices, elements);\n",row.getNumElements());
      fprintf(getLogFilePtr(),"}\n");
    }

   // In XPRESS addrows() prototype, indices and elements should be const, but
   // they're not. 
   int rc = addrows(1, row.getNumElements(), &sense, &rhs, &r,
		    start, const_cast<int *>(row.getIndices()),
		    const_cast<double *>(row.getElements())); 
   assert( rc == 0 );

   freeCachedResults();
}

//#############################################################################
// Private methods
//#############################################################################

void
OsiXprSolverInterface::gutsOfCopy( const OsiXprSolverInterface & source )
{
   source.activateMe();

   if ( source.xprProbname_ != "" ) {    // source has data
     std::ostrstream pname;
     pname << xprProbname_ << "#" << osiSerial_ <<'\0';
     xprProbname_ = pname.str();
     //    sprintf(xprProbname_, "%s#%d", source.xprProbname_, osiSerial_);
     
     //    cout << "Problem " << xprProbname_ << " copied to matrix ";
     
     
     if ( getLogFilePtr()!=NULL ) {
       fprintf(getLogFilePtr(),"{\n");
       fprintf(getLogFilePtr(),"  int matrixId = %d;\n",xprMatrixId_);
       fprintf(getLogFilePtr(),"  int iret = cpymat(\"%s\", &matrixId_);\n",xprProbname_.c_str());
       fprintf(getLogFilePtr(),"}\n");
     }
     
     int iret = cpymat(xprProbname_.c_str(), &xprMatrixId_);
     if ( iret != 0 ) getipv(N_ERRNO, &iret);
     assert( iret == 0 );
     //    cout << xprMatrixId_ << "." << endl;
     xprSaved_ = true;
   }
}


//-------------------------------------------------------------------
void
OsiXprSolverInterface::gutsOfDestructor()
{
   freeCachedResults();

   assert(matrixByRow_ == NULL);
   assert(matrixByCol_ == NULL);
   assert(colupper_    == NULL);
   assert(collower_    == NULL);
   assert(rowupper_    == NULL);
   assert(rowlower_    == NULL);
                    
   assert(rowsense_    == NULL);
   assert(rhs_         == NULL);
   assert(rowrange_    == NULL);
                    
   assert(objcoeffs_   == NULL);
                    
   assert(colsol_      == NULL);
   assert(rowsol_      == NULL);
   assert(rowprice_    == NULL);
   assert(colprice_    == NULL);
   assert(rowact_      == NULL);
   assert(vartype_     == NULL);
}

//-------------------------------------------------------------------

void
OsiXprSolverInterface::freeSolution()
{
   delete [] colsol_;	    colsol_      = NULL;
   delete [] rowsol_;       rowsol_      = NULL;
   delete [] rowact_;       rowact_      = NULL;
   delete [] rowprice_;	    rowprice_    = NULL;
   delete [] colprice_;	    colprice_    = NULL;
}

//-------------------------------------------------------------------

void
OsiXprSolverInterface::freeCachedResults()
{
   delete matrixByRow_;     matrixByRow_ = NULL;
   delete matrixByCol_;     matrixByCol_ = NULL;
   delete [] colupper_;     colupper_    = NULL;
   delete [] collower_;	    collower_    = NULL;
   delete [] rowupper_;	    rowupper_    = NULL;
   delete [] rowlower_;	    rowlower_    = NULL;

   delete [] rowsense_;	    rowsense_    = NULL;
   delete [] rhs_;	    rhs_         = NULL;
   delete [] rowrange_;	    rowrange_    = NULL;

   delete [] objcoeffs_;    objcoeffs_   = NULL;

   freeSolution();

   delete [] ivarind_;      ivarind_     = NULL;
   delete [] ivartype_;     ivartype_    = NULL;
   delete [] vartype_;      vartype_     = NULL;
}

//-------------------------------------------------------------------
// Set up lists of integer variables
//-------------------------------------------------------------------
int
OsiXprSolverInterface::getNumIntVars() const
{
  activateMe();
  
  int     nintvars = 0, nsets = 0;
  
  if ( isDataLoaded() ) {  
    
    if (getLogFilePtr()!=NULL) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"   int nintvars,nsets;\n");
      fprintf(getLogFilePtr(),"   getglobal(&nintvars, &nsets, NULL, NULL, NULL, NULL, NULL, NULL, NULL);\n");
      fprintf(getLogFilePtr(),"}\n");
    }
    getglobal(&nintvars, &nsets, 
      NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  }
  
  return nintvars;
}

void
OsiXprSolverInterface::getVarTypes() const
{
   int     nintvars = getNumIntVars();
   int     nsets;
   int     ncols = getNumCols();

   if ( vartype_ == NULL && nintvars > 0 ) {
      activateMe();
   
    if (getLogFilePtr()!=NULL) {
      fprintf(getLogFilePtr(),"{\n");
      fprintf(getLogFilePtr(),"   int nintvars = %d;\n",nintvars);
      fprintf(getLogFilePtr(),"   int nsets;\n");
      fprintf(getLogFilePtr(),"   char ivartype[%d];\n",nintvars);
      fprintf(getLogFilePtr(),"   char ivarind[%d];\n",nintvars);
      fprintf(getLogFilePtr(),"   getglobal(&nintvars, &nsets, ivartype, ivarind, NULL, NULL, NULL, NULL, NULL);\n");
      fprintf(getLogFilePtr(),"}\n");
    }

      ivartype_ = new char[nintvars];
      ivarind_  = new int[nintvars];

      getglobal(&nintvars, &nsets, 
		ivartype_, ivarind_, NULL, NULL, NULL, NULL, NULL);
      // Currently, only binary and integer vars are supported.

      vartype_  = new char[ncols];

      int     i, j;
      for ( i = j = 0;  j < ncols;  j++ ) {
	 if ( i < nintvars && j == ivarind_[i] ) {
	    vartype_[j] = ivartype_[i];
	    i++;
	 } else 
	    vartype_[j] = 'C';
      }
   }
}

//-------------------------------------------------------------------
// Make the active matrix in XPRESS-MP refer to the current solver 
// object.  If another object is current, then it needs to be saved.
//------------------------------------------------------------------- 
const OsiXprSolverInterface *
OsiXprSolverInterface::xprCurrentProblem_ = NULL;

void
OsiXprSolverInterface::activateMe() const
{
   if ( xprCurrentProblem_ == this ) return; // we're already active

   if ( xprCurrentProblem_ ) {   // someone else is active...
      if ( xprCurrentProblem_->xprProbname_ != "" ) {   // ...and has data
      
        if ( getLogFilePtr()!=NULL ) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  int matrixId;\n");
          fprintf(getLogFilePtr(),"  int iret = savmat(& matrixId);\n");
        }

	 //      cout << "Problem " << xprCurrentProblem_->xprProbname_
	 //	   << " saved as matrix ";
	 int iret = savmat(& xprCurrentProblem_->xprMatrixId_);
	 if ( iret != 0 ) getipv(N_ERRNO, &iret);
	 assert( iret == 0 );
	 //      cout << xprCurrentProblem_->xprMatrixId_ << "." << endl;
	 xprCurrentProblem_->xprSaved_ = true;
     
        if ( getLogFilePtr()!=NULL ) {
          fprintf(getLogFilePtr(),"  // matrixId returned was %d;\n",xprCurrentProblem_->xprMatrixId_);
          fprintf(getLogFilePtr(),"}\n");
        }

      }
   }

   if ( xprSaved_ ) {           // we were saved before
      //    cout << "Problem " << xprProbname_ << " restored from  matrix "
      //         << xprMatrixId_ << "." << endl;
           
        if ( getLogFilePtr()!=NULL ) {
          fprintf(getLogFilePtr(),"{\n");
          fprintf(getLogFilePtr(),"  int iret = resmat(%d);\n",xprMatrixId_);
          fprintf(getLogFilePtr(),"}\n");
        }

      int iret = resmat(xprMatrixId_);
      if ( iret != 0 ) getipv(N_ERRNO, &iret);
      assert( iret == 0 );
      xprSaved_ = false;
   }

   xprCurrentProblem_ = this;
}

//-------------------------------------------------------------------

bool
OsiXprSolverInterface::isDataLoaded() const
{
   return xprProbname_ != "";
}

//#############################################################################
