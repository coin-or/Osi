//  LAST EDIT: Thu Mar  8 14:17:15 2001 by Tobias Pfender (opt14!bzfpfend) 
//-----------------------------------------------------------------------------
// name:     OSI Interface for CPLEX
// author:   Tobias Pfender
//           Konrad-Zuse-Zentrum Berlin (Germany)
//           email: pfender@zib.de
// date:     09/25/2000
// comments: please scan this file for '???' and read the comments
//-----------------------------------------------------------------------------
// Copyright (C) 2000, Tobias Pfender, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifdef COIN_USE_CPX
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <string>

#include "CoinError.hpp"

#include "OsiCpxSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiWarmStartBasis.hpp"

//#############################################################################
// A couple of helper functions
//#############################################################################

inline void freeCacheDouble( double*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheChar( char*& ptr )
{
  if( ptr != NULL )
    {
      delete [] ptr;
      ptr = NULL;
    }
}

inline void freeCacheMatrix( OsiPackedMatrix*& ptr )
{
  if( ptr != NULL )
    {
      delete ptr;
      ptr = NULL;
    }
}

inline void checkCPXerror( int err, std::string cpxfuncname, std::string osimethod )
{
  if( err != 0 )
    {
      char s[100];
      sprintf( s, "%s returned error %d", cpxfuncname.c_str(), err );
      cout << "ERROR: " << s << " (" << osimethod << " in OsiCpxSolverInterface)" << endl;
      throw CoinError( s, osimethod.c_str(), "OsiCpxSolverInterface" );
    }
}

//#############################################################################
// Solve methods
//#############################################################################

void OsiCpxSolverInterface::initialSolve()
{
  CPXLPptr lp = getMutableLpPtr();
  // If mip problem (ie integer data exits),
  // then must change problem type.
  // CPLEX will return an error condition if this is not done
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_MIP ) 
    CPXchgprobtype(env_,lp,CPXPROB_RELAXED);

  CPXprimopt( env_, lp );

  freeCachedResults();
}
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::resolve()
{
  CPXLPptr lp = getMutableLpPtr();
  // If mip problem (ie integer data exits),
  // then must change problem type.
  // CPLEX will return an error condition if this is not done
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_MIP ) 
    CPXchgprobtype(env_,lp,CPXPROB_RELAXED);

  CPXdualopt( env_, lp );   

  freeCachedResults();
}
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::branchAndBound()
{
  CPXLPptr lp = getMutableLpPtr();
  // Problem type might have been changed from Mip to relaxed.
  // It this was done, then change it back.
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_RELAXED ) 
    CPXchgprobtype(env_,lp,CPXPROB_MIP);
  CPXmipopt( env_, lp );

  freeCachedResults();
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiCpxSolverInterface::setIntParam(OsiIntParam key, int value)
{
  bool retval = false;
  switch (key)
    {
    case OsiMaxNumIteration:
      retval = ( CPXsetintparam( env_, CPX_PARAM_ITLIM, value ) == 0 );  // ??? OsiMaxNumIteration == #Simplex-iterations ???
      break;
    case OsiMaxNumIterationHotStart:
      /* FIXME */
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
OsiCpxSolverInterface::setDblParam(OsiDblParam key, double value)
{
  bool retval = false;
  switch (key)
    {
    case OsiDualObjectiveLimit:
      if( getObjSense() == +1 )
	retval = ( CPXsetdblparam( env_, CPX_PARAM_OBJULIM, value ) == 0 ); // min
      else
	retval = ( CPXsetdblparam( env_, CPX_PARAM_OBJLLIM, value ) == 0 ); // max
      break;
    case OsiPrimalObjectiveLimit:
      if( getObjSense() == +1 )
	retval = ( CPXsetdblparam( env_, CPX_PARAM_OBJLLIM, value ) == 0 ); // min
      else
	retval = ( CPXsetdblparam( env_, CPX_PARAM_OBJULIM, value ) == 0 ); // max
      break;
    case OsiDualTolerance:
      retval = ( CPXsetdblparam( env_, CPX_PARAM_EPOPT, value ) == 0 ); // ??? OsiDualTolerance == CPLEX Optimality tolerance ???
      break;
    case OsiPrimalTolerance:
      retval = ( CPXsetdblparam( env_, CPX_PARAM_EPRHS, value ) == 0 ); // ??? OsiPrimalTolerance == CPLEX Feasibility tolerance ???
      break;
    case OsiLastDblParam:
      retval = false;
      break;
    }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiCpxSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  bool retval = false;
  switch (key)
    {
    case OsiMaxNumIteration:
      retval = ( CPXgetintparam( env_, CPX_PARAM_ITLIM, &value ) == 0 );  // ??? OsiMaxNumIteration == #Simplex-iterations ???
      break;
    case OsiMaxNumIterationHotStart:
      /* FIXME */
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
OsiCpxSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  bool retval = false;
  switch (key) 
    {
    case OsiDualObjectiveLimit:
      if( getObjSense() == +1 )
	retval = ( CPXgetdblparam( env_, CPX_PARAM_OBJULIM, &value ) == 0 ); // min
      else
	retval = ( CPXgetdblparam( env_, CPX_PARAM_OBJLLIM, &value ) == 0 ); // max
      break;
    case OsiPrimalObjectiveLimit:
      if( getObjSense() == +1 )
	retval = ( CPXgetdblparam( env_, CPX_PARAM_OBJLLIM, &value ) == 0 ); // min
      else
	retval = ( CPXgetdblparam( env_, CPX_PARAM_OBJULIM, &value ) == 0 ); // max
      break;
    case OsiDualTolerance:
      retval = ( CPXgetdblparam( env_, CPX_PARAM_EPOPT, &value ) == 0 ); // ??? OsiDualTolerance == CPLEX Optimality tolerance ???
      break;
    case OsiPrimalTolerance:
      retval = ( CPXgetdblparam( env_, CPX_PARAM_EPRHS, &value ) == 0 ); // ??? OsiPrimalTolerance == CPLEX Feasibility tolerance ???
      break;
    case OsiLastDblParam:
      retval = false;
      break;
    }
  return retval;
}

//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiCpxSolverInterface::isAbandoned() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == 0 || stat == CPX_NUM_BEST_FEAS || stat == CPX_NUM_BEST_INFEAS || 
    stat == CPX_ABORT_FEAS || stat == CPX_ABORT_INFEAS || stat == CPX_ABORT_CROSSOVER;
}

bool OsiCpxSolverInterface::isProvenOptimal() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_OPTIMAL || stat == CPX_OPTIMAL_INFEAS;
}

bool OsiCpxSolverInterface::isProvenPrimalInfeasible() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_INFEASIBLE || stat == CPX_ABORT_PRIM_INFEAS || stat == CPX_ABORT_PRIM_DUAL_INFEAS;
}

bool OsiCpxSolverInterface::isProvenDualInfeasible() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_UNBOUNDED || stat == CPX_ABORT_DUAL_INFEAS || stat == CPX_ABORT_PRIM_DUAL_INFEAS;
}

bool OsiCpxSolverInterface::isPrimalObjectiveLimitReached() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_OBJ_LIM; // ??? CPLEX just have an objective limit - not primal/dual
}

bool OsiCpxSolverInterface::isDualObjectiveLimitReached() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_OBJ_LIM; // ??? CPLEX just have an objective limit - not primal/dual
}

bool OsiCpxSolverInterface::isIterationLimitReached() const
{
  int stat = CPXgetstat( env_, getMutableLpPtr() );

  return stat == CPX_IT_LIM_FEAS || stat == CPX_IT_LIM_INFEAS;
}

//#############################################################################
// WarmStart related methods
//#############################################################################

OsiWarmStart* OsiCpxSolverInterface::getWarmStart() const
{
  OsiWarmStartBasis* ws = NULL;
  int numcols = getNumCols();
  int numrows = getNumRows();
  int *cstat = new int[numcols];
  int *rstat = new int[numrows];
  int restat, i;

  restat = CPXgetbase( env_, getMutableLpPtr(), cstat, rstat );
  if( restat == 0 )
    {
      ws = new OsiWarmStartBasis;
      ws->setSize( numcols, numrows );
      
      for( i = 0; i < numrows; ++i )
	{
	  switch( rstat[i] )
	    {
	    case CPX_BASIC:
	      ws->setArtifStatus( i, OsiWarmStartBasis::basic );
	      break;
	    case CPX_AT_LOWER:
	      ws->setArtifStatus( i, OsiWarmStartBasis::atLowerBound );
	      break;
	    case CPX_AT_UPPER:
	      ws->setArtifStatus( i, OsiWarmStartBasis::atUpperBound );
	      break;
	    default:  // unknown row status
	      delete ws;
	      ws = NULL;
	      goto TERMINATE;
	    }
	}
      for( i = 0; i < numcols; ++i )
	{
	  switch( cstat[i] )
	    {
	    case CPX_BASIC:
	      ws->setStructStatus( i, OsiWarmStartBasis::basic );
	      break;
	    case CPX_AT_LOWER:
	      ws->setStructStatus( i, OsiWarmStartBasis::atLowerBound );
	      break;
	    case CPX_AT_UPPER:
	      ws->setStructStatus( i, OsiWarmStartBasis::atUpperBound );
	      break;
	    case CPX_FREE_SUPER:
	      ws->setStructStatus( i, OsiWarmStartBasis::isFree );
	      break;
	    default:  // unknown column status
	      delete ws;
	      ws = NULL;
	      goto TERMINATE;
	    }
	}
    }
 TERMINATE:
  delete[] cstat;
  delete[] rstat;

  return ws;
}

//-----------------------------------------------------------------------------

bool OsiCpxSolverInterface::setWarmStart(const OsiWarmStart* warmstart)
{
  const OsiWarmStartBasis* ws = dynamic_cast<const OsiWarmStartBasis*>(warmstart);
  int numcols, numrows, i, restat;
  int *cstat, *rstat;
  bool retval = false;

  if( !ws )
    return false;

  numcols = ws->getNumStructural();
  numrows = ws->getNumArtificial();
  
  if( numcols != getNumCols() || numrows != getNumRows() )
    return false;

  cstat = new int[numcols];
  rstat = new int[numrows];
  for( i = 0; i < numrows; ++i )
    {
      switch( ws->getArtifStatus( i ) )
	{
	case OsiWarmStartBasis::basic:
	  rstat[i] = CPX_BASIC;
	  break;
	case OsiWarmStartBasis::atLowerBound:
	  rstat[i] = CPX_AT_LOWER;
	  break;
	case OsiWarmStartBasis::atUpperBound:
	  rstat[i] = CPX_AT_UPPER;
	  break;
	default:  // unknown row status
	  retval = false;
	  goto TERMINATE;
	}
    }
  for( i = 0; i < numcols; ++i )
    {
      switch( ws->getStructStatus( i ) )
	{
	case OsiWarmStartBasis::basic:
	  cstat[i] = CPX_BASIC;
	  break;
	case OsiWarmStartBasis::atLowerBound:
	  cstat[i] = CPX_AT_LOWER;
	  break;
	case OsiWarmStartBasis::atUpperBound:
	  cstat[i] = CPX_AT_UPPER;
	  break;
	case OsiWarmStartBasis::isFree:
	  cstat[i] = CPX_FREE_SUPER;
	  break;
	default:  // unknown row status
	  retval = false;
	  goto TERMINATE;
	}
    }

  // *FIXME* : can this be getMutableLpPtr() ? Does any cached data change by
  // *FIXME* : setting warmstart? Or at least wouldn't it be sufficient to
  // *FIXME* : clear the cached results but not the problem data?
  restat = CPXcopybase( env_, getLpPtr(), cstat, rstat );
  retval = (restat == 0);
 TERMINATE:
  delete[] cstat;
  delete[] rstat;
  return retval;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiCpxSolverInterface::markHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::markHotStart();
}

void OsiCpxSolverInterface::solveFromHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::solveFromHotStart();
}

void OsiCpxSolverInterface::unmarkHotStart()
{
  // *FIXME* : do better... -LL
  OsiSolverInterface::unmarkHotStart();
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows, columns, elements, ...
//------------------------------------------------------------------
int OsiCpxSolverInterface::getNumCols() const
{
  return CPXgetnumcols( env_, getMutableLpPtr() );
}
int OsiCpxSolverInterface::getNumRows() const
{
  return CPXgetnumrows( env_, getMutableLpPtr() );
}
int OsiCpxSolverInterface::getNumElements() const
{
  return CPXgetnumnz( env_, getMutableLpPtr() );
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  

const double * OsiCpxSolverInterface::getColLower() const
{
  if( collower_ == NULL )
    {
      int ncols = CPXgetnumcols( env_, getMutableLpPtr() );
      if( ncols > 0 )
	{
	  collower_ = new double[ncols];
	  CPXgetlb( env_, getMutableLpPtr(), collower_, 0, ncols-1 );
	}
    }
  return collower_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getColUpper() const
{
  if( colupper_ == NULL )
    {
      int ncols = CPXgetnumcols( env_, getMutableLpPtr() );
      if( ncols > 0 )
	{
	  colupper_ = new double[ncols];
	  CPXgetub( env_, getMutableLpPtr(), colupper_, 0, ncols-1 );
	}
    }
  return colupper_;
}
//------------------------------------------------------------------
const char * OsiCpxSolverInterface::getRowSense() const
{
  if ( rowsense_==NULL )
  {
    
    // rowsense is determined with rhs, so invoke rhs
    getRightHandSide();
    assert( rowsense_!=NULL || getNumRows() == 0 );
#if 0
    int err;
    int nrows = getNumRows();
    rowsense_ = new char[nrows];
    err = CPXgetsense( env_, getMutableLpPtr(), rowsense_, 0, nrows-1 );
    checkCPXerror( err, "CPXgetsense", "getRowSense" );
    // search free constraints CPLEX doesn't support: implemented as ranged rows with infinite bounds
    const double* therhs = rhs();
    
    for( int i = 0; i < nrows; ++i )
      if( rowsense_[i] == 'R' && therhs[i] <= -infinity() )
        rowsense_[i] = 'N';
#endif
  }
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRightHandSide() const
{
  if ( rhs_==NULL ) {
     CPXLPptr lp = getMutableLpPtr();
     int nrows = getNumRows();
     if( nrows > 0 ) {
	rhs_ = new double[nrows];
	CPXgetrhs( env_, lp, rhs_, 0, nrows-1 );

	assert( rowrange_ == NULL );
	rowrange_ = new double[nrows];
	CPXgetrngval( env_, lp, rowrange_, 0, nrows-1 );

	assert( rowsense_ == NULL );
	rowsense_ = new char[nrows];
	CPXgetsense( env_, lp, rowsense_, 0, nrows-1 );

	double inf = getInfinity();
	int i;
	for ( i = 0; i < nrows; ++i ) {  
	   if ( rowsense_[i] != 'R' ) {
	      rowrange_[i]=0.0;
	   } else {
	      if ( rhs_[i] <= -inf ) {
		 rowsense_[i] = 'N';
		 rowrange_[i] = 0.0;
		 rhs_[i] = 0.0;
	      } else {
		 rhs_[i] = rhs_[i] + rowrange_[i];
	      }
	   }
	}
     }
  }
  return rhs_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRowRange() const
{
  if ( rowrange_==NULL ) 
    {
      // rowrange is determined with rhs, so invoke rhs
      getRightHandSide();
      assert( rowrange_!=NULL || getNumRows() == 0 );
    }
  return rowrange_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRowLower() const
{
  if ( rowlower_ == NULL )
    {
      int     nrows = getNumRows();
      const   char    *rowsense = getRowSense();
      const   double  *rhs      = getRightHandSide();
      const   double  *rowrange = getRowRange();
    
      if ( nrows > 0 )
	{
	  rowlower_ = new double[nrows];
	  
	  double dum1;
	  for ( int i = 0;  i < nrows;  i++ )
	    convertSenseToBound( rowsense[i], rhs[i], rowrange[i],
				 rowlower_[i], dum1 );
	}
    }
  return rowlower_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRowUpper() const
{  
  if ( rowupper_ == NULL )
    {
      int     nrows = getNumRows();
      const   char    *rowsense = getRowSense();
      const   double  *rhs      = getRightHandSide();
      const   double  *rowrange = getRowRange();
      
      if ( nrows > 0 ) 
	{
	  rowupper_ = new double[nrows];
	  
	  double dum1;
	  for ( int i = 0;  i < nrows;  i++ )
	    convertSenseToBound( rowsense[i], rhs[i], rowrange[i],
				 dum1, rowupper_[i] );
	}
    }
  
  return rowupper_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getObjCoefficients() const
{
  if ( obj_==NULL )
    {
      int ncols = CPXgetnumcols( env_, getMutableLpPtr() );
      if( ncols > 0 )
	{
	  obj_ = new double[ncols];
	  int err = CPXgetobj( env_, getMutableLpPtr(), obj_, 0, ncols-1 );
	  checkCPXerror( err, "CPXgetobj", "getObjCoefficients" );
	}
    }
  return obj_;
}
//------------------------------------------------------------------
double OsiCpxSolverInterface::getObjSense() const
{
  if( CPXgetobjsen( env_, getMutableLpPtr() ) == CPX_MIN )
    return +1.0;
  else
    return -1.0;
}

//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------

bool OsiCpxSolverInterface::isContinuous( int colNumber ) const
{
  return( getCtype()[colNumber] == CPX_CONTINUOUS );
}

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------

const OsiPackedMatrix * OsiCpxSolverInterface::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL ) 
    {
      int nrows = getNumRows();
      int ncols = getNumCols();
      int nelems;
      int *starts   = new int   [nrows + 1];
      int *len      = new int   [nrows];
      
      int requiredSpace;
      int rc = CPXgetrows( env_, getMutableLpPtr(), 
			   &nelems, starts, NULL, NULL, 0, &requiredSpace,
			   0, nrows-1 );
      
      assert( -requiredSpace == getNumElements() );
      int     *indices  = new int   [-requiredSpace];
      double  *elements = new double[-requiredSpace]; 
      
      rc = CPXgetrows( env_, getMutableLpPtr(), 
		       &nelems, starts, indices, elements, -requiredSpace,
		       &requiredSpace, 0, nrows-1 );
      assert( requiredSpace == 0 );
            
      matrixByRow_ = new OsiPackedMatrix();
      
      // Should be able to pass null for length of packed matrix,
      // assignMatrix does not seem to allow (even though documentation
      // say it is possible to do this). 
      // For now compute the length.
      starts[nrows] = nelems;
      for ( int i = 0; i < nrows; ++i )
	len[i]=starts[i+1] - starts[i];
      
      matrixByRow_->assignMatrix( false /* not column ordered */,
				  ncols, nrows, nelems,
				  elements, indices, starts, len /*NULL*/);
      
    }
  return matrixByRow_;
} 

//------------------------------------------------------------------

const OsiPackedMatrix * OsiCpxSolverInterface::getMatrixByCol() const
{
  if ( matrixByCol_ == NULL )
    {
      int nrows = getNumRows();
      int ncols = getNumCols();
      int nelems;
      int *starts = new int   [ncols + 1];
      int *len    = new int   [ncols];
      
      int requiredSpace;
      int rc = CPXgetcols( env_, getMutableLpPtr(), 
			   &nelems, starts, NULL, NULL, 0, &requiredSpace,
			   0, ncols-1 );
      assert( -requiredSpace == getNumElements() );
      
      int     *indices  = new int   [-requiredSpace];
      double  *elements = new double[-requiredSpace]; 
      
      rc = CPXgetcols( env_, getMutableLpPtr(), 
		       &nelems, starts, indices, elements, -requiredSpace,
		       &requiredSpace, 0, ncols-1 );
      assert( requiredSpace == 0);
      
      matrixByCol_ = new OsiPackedMatrix();
      
      // Should be able to pass null for length of packed matrix,
      // assignMatrix does not seem to allow (even though documentation
      // say it is possible to do this). 
      // For now compute the length.
      starts[ncols] = nelems;
      for ( int i = 0; i < ncols; i++ )
	len[i]=starts[i+1] - starts[i];
      
      matrixByCol_->assignMatrix( true /* column ordered */,
				  nrows, ncols, nelems,
				  elements, indices, starts, len /*NULL*/);
      assert( matrixByCol_->getNumCols()==ncols );
      assert( matrixByCol_->getNumRows()==nrows );
    }
  return matrixByCol_;
} 

//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double OsiCpxSolverInterface::getInfinity() const
{
  return CPX_INFBOUND;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

// *FIXME*: what should be done if a certain vector doesn't exist???

const double * OsiCpxSolverInterface::getColSolution() const
{
  if( colsol_==NULL )
    {
      CPXLPptr lp = getMutableLpPtr();
      int ncols = CPXgetnumcols( env_, lp );
      if( ncols > 0 )
	{
	  colsol_ = new double[ncols]; 
	  int probType = CPXgetprobtype(env_,lp);
	  if ( probType == CPXPROB_MIP ) {
	    int err = CPXgetmipx( env_, lp, colsol_, 0, ncols-1 );
	    if ( err == CPXERR_NO_INT_SOLN ) 
	      CoinFillN( colsol_, ncols, 0.0 );
	    else
	      checkCPXerror( err, "CPXgetmipx", "getColSolution" );
	  } else {
	    int err = CPXgetx( env_, lp, colsol_, 0, ncols-1 );
	    if ( err == CPXERR_NO_SOLN ) 
	      CoinFillN( colsol_, ncols, 0.0 );
	    else
	      checkCPXerror( err, "CPXgetx", "getColSolution" );
	  }
	}
    }
  return colsol_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRowPrice() const
{
  if( rowsol_==NULL )
    {
      int nrows = getNumRows();
      if( nrows > 0 )
	{
	  rowsol_ = new double[nrows];
	  int err = CPXgetpi( env_, getMutableLpPtr(), rowsol_, 0, nrows-1 );
	  if ( err == CPXERR_NO_SOLN ) 
	    CoinFillN( rowsol_, nrows, 0.0 );
	  else
	    checkCPXerror( err, "CPXgetpi", "getRowPrice" );
	}
    }
  return rowsol_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getReducedCost() const
{
  if( redcost_==NULL )
    {
      int ncols = CPXgetnumcols( env_, getMutableLpPtr() );
      if( ncols > 0 )
	{
	  redcost_ = new double[ncols]; 
	  int err = CPXgetdj( env_, getMutableLpPtr(), redcost_, 0, ncols-1 );
	  if ( err == CPXERR_NO_SOLN ) 
	    CoinFillN( redcost_, ncols, 0.0 );
	  else
	    checkCPXerror( err, "CPXgetdj", "getReducedCost" );
	}
    }
  return redcost_;
}
//------------------------------------------------------------------
const double * OsiCpxSolverInterface::getRowActivity() const
{
  // *FIXME* : this can be returned for integer programs, just use 
  // *FIXME* : CPXgetmipslack and subtract from the right-hand-side
  if( rowact_==NULL )
    {
      int nrows = getNumRows();
      if( nrows > 0 )
	{
	  rowact_ = new double[nrows];
	  int err = CPXgetax( env_, getMutableLpPtr(), rowact_, 0, nrows-1 );
	  if ( err == CPXERR_NO_SOLN )
	    CoinFillN( rowact_, nrows, 0.0 );
	  else
	    checkCPXerror( err, "CPXgetax", "getRowActivity" );
	}
    }
  return rowact_;
}
//------------------------------------------------------------------
double OsiCpxSolverInterface::getObjValue() const
{
  double objval = 0.0;
  int err;

  CPXLPptr lp = getMutableLpPtr();
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_MIP ) {
    err = CPXgetmipobjval( env_, lp, &objval);
    if ( err == CPXERR_NO_INT_SOLN ) 
      // => return 0.0 as objective value (?? is this the correct behaviour ??)
      objval = 0.0;
    else
      checkCPXerror( err, "CPXgetmipobjval", "getObjValue" );
  } else {
    err = CPXgetobjval( env_, lp, &objval );
    if ( err == CPXERR_NO_SOLN ) 
      objval = 0.0;
    else
      checkCPXerror( err, "CPXgetobjval", "getObjValue" );
  }
  return objval;
}
//------------------------------------------------------------------
int OsiCpxSolverInterface::getIterationCount() const
{
  return CPXgetitcnt( env_, getMutableLpPtr() );
}
//------------------------------------------------------------------
std::vector<double*> OsiCpxSolverInterface::getDualRays(int maxNumRays) const
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getDualRays",
		 "OsiCpxSolverInterface");
  return std::vector<double*>();
}
//------------------------------------------------------------------
std::vector<double*> OsiCpxSolverInterface::getPrimalRays(int maxNumRays) const
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getPrimalRays",
		 "OsiCpxSolverInterface");
  return std::vector<double*>();
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

void setObjCoeff( int elementIndex, double elementValue )
{
  int err = CPXchgobj(env_, getLpPtr(), 1, &elementIndex, &elementValue);
  checkCPXerror(err, "CPXchgobj", "setObjCoeff");

  freeCachedColRim();
}

void OsiCpxSolverInterface::setColLower(int elementIndex, double elementValue)
{
  char c = 'L';
  int err = CPXchgbds( env_, getLpPtr(), 1, &elementIndex, &c, &elementValue );
  checkCPXerror( err, "CPXchgbds", "setColLower" );

  freeCachedColRim();
}
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::setColUpper(int elementIndex, double elementValue)
{  
  char c = 'U';
  int err = CPXchgbds( env_, getLpPtr(), 1, &elementIndex, &c, &elementValue );
  checkCPXerror( err, "CPXchgbds", "setColUpper" );

  freeCachedColRim();
} 
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::setColBounds( int elementIndex, double lower, double upper )
{
  char c[2] = { 'L', 'U' };
  int ind[2];
  double bd[2];
  int err;

  ind[0] = elementIndex;
  ind[1] = elementIndex;
  bd[0] = lower;
  bd[1] = upper;
  err = CPXchgbds( env_, getLpPtr(), 2, ind, c, bd );
  checkCPXerror( err, "CPXchgbds", "setColBounds" );

  freeCachedColRim();
}
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  OsiSolverInterface::setColSetBounds( indexFirst, indexLast, boundList );
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setRowLower( int i, double elementValue )
{
  double rhs   = getRightHandSide()[i];
  double range = getRowRange()[i];
  char   sense = getRowSense()[i];
  double lower, upper;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( lower != elementValue ) {
      convertBoundToSense( elementValue, upper, sense, rhs, range );
      setRowType( i, sense, rhs, range );
      // freeCachedRowRim(); --- invoked in setRowType()
    }
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setRowUpper( int i, double elementValue )
{
  double rhs   = getRightHandSide()[i];
  double range = getRowRange()[i];
  char   sense = getRowSense()[i];
  double lower, upper;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( upper != elementValue ) {
      convertBoundToSense( lower, elementValue, sense, rhs, range );
      setRowType( i, sense, rhs, range );
      // freeCachedRowRim(); --- invoked in setRowType()
    }
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
  double rhs, range;
  char sense;
  
  convertBoundToSense( lower, upper, sense, rhs, range );
  setRowType( elementIndex, sense, rhs, range );
  // freeCachedRowRim(); --- invoked in setRowType()
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  int err;

  if (sense == 'R') {
     rightHandSide -= range;
  }

  //  if( sense == 'N' )
  //    {
  //      sense = 'R';
  //      rightHandSide = -getInfinity();
  //      range = getInfinity();
  //    }
  err = CPXchgsense( env_, getMutableLpPtr(), 1, &i, &sense );
  checkCPXerror( err, "CPXchgsense", "setRowType" );
  err = CPXchgrhs( env_, getMutableLpPtr(), 1, &i, &rightHandSide );
  checkCPXerror( err, "CPXchgrhs", "setRowType" );
  err = CPXchgrngval( env_, getMutableLpPtr(), 1, &i, &range );
  checkCPXerror( err, "CPXchgrngval", "setRowType" );

  freeCachedRowRim();
}
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  OsiSolverInterface::setRowSetBounds( indexFirst, indexLast, boundList );
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  OsiSolverInterface::setRowSetTypes( indexFirst, indexLast, senseList, rhsList, rangeList );
}
//#############################################################################
// *FIXME*: I don't know CPLEX that well, there might be a better way to do
// *FIXME*: this... -LL
void
OsiCpxSolverInterface::setContinuous(int index)
{
  CPXLPptr lp = getMutableLpPtr();
  int probType = CPXgetprobtype(env_,lp);
  if ( probType != CPXPROB_LP ) {
    int err;
    char type = 'C';
    err = CPXchgctype(env_, lp, 1, &index, &type);
    checkCPXerror( err, "CPXchgctype", "setContinuous");
    freeCachedColRim();
  }
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setInteger(int index)
{
  CPXLPptr lp = getMutableLpPtr();
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_LP ) 
    CPXchgprobtype(env_,lp,CPXPROB_MIP);

  int err;
  char type = 'I';
  if (getColLower()[index] == 0.0 && getColUpper()[index] == 1.0)
    type = 'B';
  err = CPXchgctype(env_, lp, 1, &index, &type);
  checkCPXerror( err, "CPXchgctype", "setInteger");

  freeCachedColRim();
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setContinuous(const int* indices, int len)
{
  CPXLPptr lp = getMutableLpPtr();
  int probType = CPXgetprobtype(env_,lp);
  if ( probType != CPXPROB_LP ) {
    int err;
    char* type = new char[len];
    CoinFillN(type, len, 'C');
    err = CPXchgctype(env_, lp, len, const_cast<int*>(indices), type);
    checkCPXerror( err, "CPXchgctype", "setContinuous");
    delete[] type;
    freeCachedColRim();
  }
}
//-----------------------------------------------------------------------------
void
OsiCpxSolverInterface::setInteger(const int* indices, int len)
{
  CPXLPptr lp = getMutableLpPtr();
  int probType = CPXgetprobtype(env_,lp);
  if ( probType == CPXPROB_LP ) 
    CPXchgprobtype(env_,lp,CPXPROB_MIP);

  int err;
  char* type = new char[len];
  CoinFillN(type, len, 'I');
  const double* clb = getColLower();
  const double* cub = getColUpper();
  for (int i = 0; i < len; ++i) {
    if (clb[indices[i]] == 0.0 && cub[indices[i]] == 1.0)
      type[i] = 'B';
  }
  CPXENVptr env = env_;
  err = CPXchgctype(env, lp, len, const_cast<int*>(indices), type);
  checkCPXerror( err, "CPXchgctype", "setInteger");
  delete[] type;

  freeCachedColRim();
}
//#############################################################################

void OsiCpxSolverInterface::setObjSense(double s) 
{
  if( s == +1.0 )
    CPXchgobjsen( env_, getLpPtr(), CPX_MIN );
  else
    CPXchgobjsen( env_, getLpPtr(), CPX_MAX );
}
 
//-----------------------------------------------------------------------------

void OsiCpxSolverInterface::setColSolution(const double * cs) 
{
  int nc = getNumCols();

  if( cs == NULL )
    freeCachedResults();
  else if( nc > 0 )
    {
      // If colsol isn't allocated, then allocate it
      if ( colsol_==NULL ) {
	int nr = getNumRows();
	colsol_ = new double[nc];
	rowsol_ = new double[nr];
	redcost_ = new double[nc];
	rowact_ = new double[nr];
      }
	
      // Copy in new col solution.
      CoinDisjointCopyN( cs, nc, colsol_ );
      
      // CPLEX < 7.0 doesn't support setting a col solution without a row solution
      // -> if a row solution exists or CPLEX version >= 7, then pass into CPLEX
      if ( rowsol_ != NULL || cpxVersionMajor_ >= 7 )
	{
	  int err = CPXcopystart( env_, getMutableLpPtr(), NULL, NULL, 
				  const_cast<double*>( colsol_ ), 
				  const_cast<double*>( rowsol_ ), 
				  NULL, NULL );
	  checkCPXerror( err, "CPXcopystart", "setColSolution" );
	}
    }
}

//-----------------------------------------------------------------------------

void OsiCpxSolverInterface::setRowPrice(const double * rs) 
{
  int nr = getNumRows();

  if( rs == NULL )
    freeCachedResults();
  else if( nr > 0 )
    {
      // If rowsol isn't allocated, then allocate it
      if ( rowsol_==NULL )  {
	int nc = getNumCols();
	colsol_ = new double[nc];
	rowsol_ = new double[nr];
	redcost_ = new double[nc];
	rowact_ = new double[nr];
      }

      // Copy in new row solution.
      CoinDisjointCopyN( rs, nr, rowsol_);
      
      // if a col solution exists, then pass into CPLEX
      if ( colsol_ != NULL )
	{
	  int err = CPXcopystart( env_, getMutableLpPtr(), NULL, NULL, 
				  const_cast<double*>( colsol_ ), 
				  const_cast<double*>( rowsol_ ), 
				  NULL, NULL );
	  checkCPXerror( err, "CPXcopystart", "setRowPrice" );
	}
    }
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiCpxSolverInterface::addCol(const OsiPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  freeCachedColRim();
  freeCachedMatrix();
  
  int err;
  int cmatbeg = 0;

  err = CPXaddcols( env_, getMutableLpPtr(),
		    1, vec.getNumElements(), const_cast<double*>(&obj),
		    &cmatbeg,
		    const_cast<int*>(vec.getIndices()),
		    const_cast<double*>(vec.getElements()),
		    const_cast<double*>(&collb),
		    const_cast<double*>(&colub), NULL );
  checkCPXerror( err, "CPXaddcols", "addCol" );
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::addCols(const int numcols,
			       const OsiPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  freeCachedColRim();
  freeCachedMatrix();
  
  int i;
  for( i = 0; i < numcols; ++i )
    addCol( *(cols[i]), collb[i], colub[i], obj[i] );
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  freeCachedColRim();
  freeCachedMatrix();
  
  int ncols = getNumCols();
  int *delstat = new int[ncols];
  int i, err;

  CoinFillN(delstat, ncols, 0);
  for( i = 0; i < num; ++i )
    delstat[columnIndices[i]] = 1;
  err = CPXdelsetcols( env_, getMutableLpPtr(), delstat );
  checkCPXerror( err, "CPXdelsetcols", "deleteCols" );
  delete[] delstat;
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  freeCachedRowRim();
  freeCachedMatrix();

  char sense;
  double rhs, range;

  convertBoundToSense( rowlb, rowub, sense, rhs, range );
  addRow( vec, sense, rhs, range );
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  freeCachedRowRim();
  freeCachedMatrix();

  int err;
  int rmatbeg = 0;

  double rhs;
  double range;
  char sense = rowsen;
  if (rowsen == 'R') {
     rhs = rowrhs - rowrng;
     range = rowrng;
  } else {
     rhs = rowrhs;
     range = 0.0;
  }

  err = CPXaddrows( env_, getMutableLpPtr(), 0, 1, vec.getNumElements(), 
		    &rhs,
		    &sense,
		    &rmatbeg,
		    const_cast<int*>(vec.getIndices()),
		    const_cast<double*>(vec.getElements()),
		    NULL, NULL );
  checkCPXerror( err, "CPXaddrows", "addRow" );
  if( rowsen == 'R' )
    {
      int row = getNumRows() - 1;
      err = CPXchgrngval( env_, getMutableLpPtr(), 1, &row, &range );
      checkCPXerror( err, "CPXchgrngval", "addRow" );
    }
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  freeCachedRowRim();
  freeCachedMatrix();

  int i;

  for( i = 0; i < numrows; ++i )
    addRow( *(rows[i]), rowlb[i], rowub[i] );
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  freeCachedRowRim();
  freeCachedMatrix();

  int i;

  for( i = 0; i < numrows; ++i )
    addRow( *(rows[i]), rowsen[i], rowrhs[i], rowrng[i] );
}
//-----------------------------------------------------------------------------
void 
OsiCpxSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  freeCachedRowRim();
  freeCachedMatrix();

  int nrows = getNumRows();
  int *delstat = new int[nrows];
  int i, err;

  CoinFillN(delstat, nrows, 0);
  for( i = 0; i < num; ++i )
    delstat[rowIndices[i]] = 1;
  err = CPXdelsetrows( env_, getMutableLpPtr(), delstat );
  checkCPXerror( err, "CPXdelsetrows", "deleteRows" );
  delete[] delstat;
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiCpxSolverInterface::loadProblem( const OsiPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const double* rowlb, const double* rowub )
{
  freeCachedData();

  const double inf = getInfinity();
  
  int nrows = matrix.getNumRows();
  char   * rowSense = new char  [nrows];
  double * rowRhs   = new double[nrows];
  double * rowRange = new double[nrows];
  
  int i;
  for ( i = nrows - 1; i >= 0; --i ) {
     const double lower = rowlb ? rowlb[i] : -inf;
     const double upper = rowub ? rowub[i] : inf;
     rowRange[i] = 0.0;
     if (lower > -inf) {
	rowRhs[i] = lower;
	if (upper < inf) {
	   if (upper==lower) {
	      rowSense[i] = 'E';
	   } else {
	      rowSense[i] = 'R';
	      rowRange[i] = upper - lower;
	   }
	} else {
	   rowSense[i] = 'G';
	}
     } else {
	if (upper < inf) {
	   rowSense[i] = 'L';
	   rowRhs[i] = upper;
	} else {
	   rowSense[i] = 'N';
	   rowRhs[i] = 0.0;
	}
     }
  }

  loadProblem(matrix, collb, colub, obj, rowSense, rowRhs, rowRange ); 
  delete [] rowSense;
  delete [] rowRhs;
  delete [] rowRange;
}
			    
//-----------------------------------------------------------------------------

void
OsiCpxSolverInterface::assignProblem( OsiPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      double*& rowlb, double*& rowub )
{
  freeCachedData();

  loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
  delete matrix;   matrix = 0;
  delete[] collb;  collb = 0;
  delete[] colub;  colub = 0;
  delete[] obj;    obj = 0;
  delete[] rowlb;  rowlb = 0;
  delete[] rowub;  rowub = 0;
}

//-----------------------------------------------------------------------------

void
OsiCpxSolverInterface::loadProblem( const OsiPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const char* rowsen, const double* rowrhs,
				    const double* rowrng )
{
  freeCachedData();

  int nc=matrix.getNumCols();
  int nr=matrix.getNumRows();

  if( nr == 0 || nc == 0 )   // empty LP?
    gutsOfDestructor();      // -> kill old LP
  else
    {
      assert( rowsen != NULL );
      assert( rowrhs != NULL );
      
      int i;
      
      // Set column values to defaults if NULL pointer passed
      double * clb;  
      double * cub;
      double * ob;
      double * rr = NULL;
      double * rhs;
      if ( collb!=NULL ) {
	clb=const_cast<double*>(collb);
      } else {
	clb = new double[nc];
	OsiFillN(clb, nc, 0.0);
      }
      if ( colub!=NULL ) {
	cub=const_cast<double*>(colub);
      } else {
	cub = new double[nc];
	OsiFillN(cub, nc, getInfinity());
      }
      if ( obj!=NULL ) {
	ob=const_cast<double*>(obj);
      } else {
	ob = new double[nc];
	OsiFillN(ob, nc, 0.0);
      } 
      if ( rowrng != NULL ) {
	rhs = new double[nr];
	rr = new double[nr];
	for ( i=0; i<nr; i++ ) {
	  if (rowsen[i] == 'R') {
	    rhs[i] = rowrhs[i] - rowrng[i];
	    rr[i] = rowrng[i];
	  } else {
	    rhs[i] = rowrhs[i];
	    rr[i] = 0.0;
	  }
	}
      } else {
	rhs = const_cast<double*>(rowrhs);
      }

      bool freeMatrixRequired = false;
      OsiPackedMatrix * m = NULL;
      if ( !matrix.isColOrdered() ) 
	{
	  m = new OsiPackedMatrix();
	  m->reverseOrderedCopyOf(matrix);
	  freeMatrixRequired = true;
	} 
      else 
	m = const_cast<OsiPackedMatrix *>(&matrix);
      
      assert( nc == m->getNumCols() );
      assert( nr == m->getNumRows() );
      assert( m->isColOrdered() ); 
      
      int objDirection = CPXgetobjsen( env_, getMutableLpPtr() );
      
      int err = CPXcopylp( env_, getLpPtr(), 
			   nc, nr,
			   // Leave ObjSense alone(set to current value).
			   objDirection,
			   ob, 
			   rhs,
			   const_cast<char *>(rowsen),
			   const_cast<int *>(m->getVectorStarts()),
			   const_cast<int *>(m->getVectorLengths()),
			   const_cast<int *>(m->getIndices()),
			   const_cast<double *>(m->getElements()),
			   const_cast<double *>(clb), 
			   const_cast<double *>(cub), 
			   rr );
      checkCPXerror( err, "CPXcopylp", "loadProblem" );
            
      if ( collb == NULL )
	delete[] clb;
      if ( colub == NULL ) 
	delete[] cub;
      if ( obj   == NULL )
	delete[] ob;
      if ( rowrng != NULL ) {
	delete[] rr;
	delete[] rhs;
      }
      
      if ( freeMatrixRequired ) 
	delete m;
    }
}
   
//-----------------------------------------------------------------------------

void
OsiCpxSolverInterface::assignProblem( OsiPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      char*& rowsen, double*& rowrhs,
				      double*& rowrng )
{
  freeCachedData();

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
// Read mps files
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::readMps( const char * filename,
				     const char * extension )
{
  freeCachedData();

  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + "." + e;
  int err = CPXreadcopyprob( env_, getLpPtr(), const_cast<char*>( fullname.c_str() ), NULL );
  checkCPXerror( err, "CPXreadcopyprob", "readMps" );
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------
void OsiCpxSolverInterface::writeMps( const char * filename,
				      const char * extension ) const
{
  char filetype[4] = "MPS";
  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + "." + e;
  int err = CPXwriteprob( env_, getMutableLpPtr(), const_cast<char*>( fullname.c_str() ), filetype );
  checkCPXerror( err, "CPXwriteprob", "writeMps" );
}

//#############################################################################
// CPX specific public interfaces
//#############################################################################

CPXENVptr OsiCpxSolverInterface::getEnvironmentPtr()
{
  assert( env_ != NULL );
  return env_;
}

CPXLPptr OsiCpxSolverInterface::getLpPtr()
{
  freeCachedData();
  return getMutableLpPtr();
}

//-----------------------------------------------------------------------------

const char * OsiCpxSolverInterface::getCtype() const
{
  if ( ctype_==NULL )
    {
      int ncols = CPXgetnumcols( env_, getMutableLpPtr() );
      if( ncols > 0 )
	{
	  ctype_ = new char[ncols];
	  int err = CPXgetctype( env_, getMutableLpPtr(), ctype_, 0, ncols-1 );
	  if ( err == CPXERR_NOT_MIP ) 
	    CoinFillN( ctype_, ncols, 'C' );
	  else
	    checkCPXerror( err, "CPXgetctype", "getCtype" );
	}
    }
  return ctype_;
}

//#############################################################################
// Static instance counter methods
//#############################################################################

void OsiCpxSolverInterface::incrementInstanceCounter()
{
  if ( numInstances_ == 0 )
    {
      int err;
      env_ = CPXopenCPLEXdevelop( &err );
      checkCPXerror( err, "CPXopenCPLEXdevelop", "incrementInstanceCounter" );
      assert( env_ != NULL );
      CPXsetintparam( env_, CPX_PARAM_SCRIND, CPX_ON ); // for testing purposes
      char logfileName[]="cplex.log";
      char filemode[]="a+";
      CPXFILEptr fp = CPXfopen( logfileName, filemode );
      CPXsetlogfile( env_, fp );
      err = sscanf( CPXversion( env_ ), "%d.%d.%d", &cpxVersionMajor_, &cpxVersionMinor_, &cpxVersionMinorMinor_ );
      assert( err == 3 );
    }
  numInstances_++;
}

//-----------------------------------------------------------------------------

void OsiCpxSolverInterface::decrementInstanceCounter()
{
  assert( numInstances_ != 0 );
  numInstances_--;
  if ( numInstances_ == 0 )
    {
      int err = CPXcloseCPLEX( &env_ );
      checkCPXerror( err, "CPXcloseCPLEX", "decrementInstanceCounter" );
      env_ = NULL;
    }
}

//-----------------------------------------------------------------------------

unsigned int OsiCpxSolverInterface::getNumInstances()
{
  return numInstances_;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiCpxSolverInterface::OsiCpxSolverInterface ()
  : OsiSolverInterface(),
    lp_(NULL),
    obj_(NULL),
    collower_(NULL),
    colupper_(NULL),
    ctype_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    rowlower_(NULL),
    rowupper_(NULL),
    colsol_(NULL),
    rowsol_(NULL),
    redcost_(NULL),
    rowact_(NULL),
    matrixByRow_(NULL),
    matrixByCol_(NULL)
{
  incrementInstanceCounter();
  gutsOfConstructor();
}


//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiSolverInterface * OsiCpxSolverInterface::clone() const
{
  return( new OsiCpxSolverInterface( *this ) );
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiCpxSolverInterface::OsiCpxSolverInterface( const OsiCpxSolverInterface & source )
  : OsiSolverInterface(source),
    lp_(NULL),
    obj_(NULL),
    collower_(NULL),
    colupper_(NULL),
    ctype_(NULL),
    rowsense_(NULL),
    rhs_(NULL),
    rowrange_(NULL),
    rowlower_(NULL),
    rowupper_(NULL),
    colsol_(NULL),
    rowsol_(NULL),
    redcost_(NULL),
    rowact_(NULL),
    matrixByRow_(NULL),
    matrixByCol_(NULL)
{
  incrementInstanceCounter();  
  gutsOfConstructor();
  gutsOfCopy( source );
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiCpxSolverInterface::~OsiCpxSolverInterface ()
{
  gutsOfDestructor();
  decrementInstanceCounter();
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiCpxSolverInterface& OsiCpxSolverInterface::operator=( const OsiCpxSolverInterface& rhs )
{
  if (this != &rhs)
    {    
      OsiSolverInterface::operator=( rhs );
      gutsOfDestructor();
      gutsOfConstructor();
      if ( rhs.lp_ !=NULL )
	gutsOfCopy( rhs );
    }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiCpxSolverInterface::applyColCut( const OsiColCut & cc )
{
  const double * cplexColLB = getColLower();
  const double * cplexColUB = getColUpper();
  const OsiPackedVector & lbs = cc.lbs();
  const OsiPackedVector & ubs = cc.ubs();
  int i;

  for( i = 0; i < lbs.getNumElements(); ++i ) 
    if ( lbs.getElements()[i] > cplexColLB[lbs.getIndices()[i]] )
      setColLower( lbs.getIndices()[i], lbs.getElements()[i] );
  for( i = 0; i < ubs.getNumElements(); ++i )
    if ( ubs.getElements()[i] < cplexColUB[ubs.getIndices()[i]] )
      setColUpper( ubs.getIndices()[i], ubs.getElements()[i] );
}

//-----------------------------------------------------------------------------

void OsiCpxSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  int err = 0;
  double rhs = 0.0;
  double rng = 0.0;
  char sns;
  double lb = rowCut.lb();
  double ub = rowCut.ub();
  if( lb <= -getInfinity() && ub >= getInfinity() )   // free constraint
    {
      rhs = -getInfinity();
      rng = 2*getInfinity();  // CPLEX doesn't support free constraints
      sns = 'R';           // -> implement them as ranged rows with infinite bounds
    }
  else if( lb <= -getInfinity() )  // <= constraint
    {
      rhs = ub;
      sns = 'L';
    }
  else if( ub >= getInfinity() )  // >= constraint
    {
      rhs = lb;
      sns = 'G';
    }
  else if( ub == lb )  // = constraint
    {
      rhs = ub;
      sns = 'E';
    }
  else  // range constraint
    {
      rhs = lb;
      rng = ub - lb;
      sns = 'R';
    }
  int rmatbeg = 0;
  err = CPXaddrows( env_, getLpPtr(), 0, 1, rowCut.row().getNumElements(),
		    &rhs, &sns, &rmatbeg, 
		    const_cast<int*>( rowCut.row().getIndices() ), 
		    const_cast<double*>( rowCut.row().getElements() ),
		    NULL, NULL );
  checkCPXerror( err, "CPXaddrows", "applyRowCut" );
  if( sns == 'R' )
    {
      err = CPXchgcoef( env_, getLpPtr(), CPXgetnumrows(env_, getLpPtr())-1,
			-2, rng );
      checkCPXerror( err, "CPXchgcoef", "applyRowCut" );
    }
}

//#############################################################################
// Private methods (non-static and static) and static data
//#############################################################################

//------------------------------------------------------------------
// Static data
//------------------------------------------------------------------      
CPXENVptr OsiCpxSolverInterface::env_ = NULL;

int OsiCpxSolverInterface::cpxVersionMajor_ = 0;
int OsiCpxSolverInterface::cpxVersionMinor_ = 0;
int OsiCpxSolverInterface::cpxVersionMinorMinor_ = 0;
unsigned int OsiCpxSolverInterface::numInstances_ = 0;
 
//-------------------------------------------------------------------
// Get pointer to CPXLPptr.
// const methods should use getMutableLpPtr().
// non-const methods should use lpPtr().
//------------------------------------------------------------------- 
CPXLPptr OsiCpxSolverInterface::getMutableLpPtr() const
{
  if ( lp_ == NULL )
    {
      int err;
      char pn[] = "OSI_CPLEX";
      assert(env_ != NULL);
      lp_ = CPXcreateprob( env_, &err, pn );
      checkCPXerror( err, "CPXcreateprob", "getMutableLpPtr" );
      assert( lp_ != NULL ); 
    }
  return lp_;
}

//-------------------------------------------------------------------

void OsiCpxSolverInterface::gutsOfCopy( const OsiCpxSolverInterface & source )
{
  // Set Objective Sense
  setObjSense(source.getObjSense());

  // Set Rim and constraints
  const double* obj = source.getObjCoefficients();
  const double* rhs = source.getRightHandSide();
  const char* sense = source.getRowSense();
  const OsiPackedMatrix * cols = source.getMatrixByCol();
  const double* lb = source.getColLower();
  const double* ub = source.getColUpper();
  loadProblem(*cols,lb,ub,obj,sense,rhs,source.getRowRange());

  // Set MIP information
  const char * colType = source.getCtype();
  if ( colType != NULL )
    CPXcopyctype( env_, getLpPtr(), const_cast<char *>(colType) );
  
  // Set Solution
  setColSolution(source.getColSolution());
  setRowPrice(source.getRowPrice());

  // Should also copy row and col names.
#if 0
  char** cname = new char*[numcols];
  char* cnamestore = NULL;
  int surplus;
  err = CPXgetcolname( env_, source.lp_, cname, NULL, 0, &surplus, 0, numcols-1 );
  if( err != CPXERR_NO_NAMES )
    {
      cnamestore = new char[-surplus];
      err = CPXgetcolname( env_, source.lp_, cname, cnamestore, -surplus, &surplus, 0, numcols-1 );
      checkCPXerror( err, "CPXgetcolname", "gutsOfCopy" );
      assert( surplus == 0 );
    }
  else
    {
      delete [] cname;
      cname = NULL;
    }
  
  char** rname = new char*[numrows];
  char* rnamestore = NULL;
  err = CPXgetrowname( env_, source.lp_, rname, NULL, 0, &surplus, 0, numrows-1 );
  if( err != CPXERR_NO_NAMES )
    {
      rnamestore = new char[-surplus];
      err = CPXgetrowname( env_, source.lp_, rname, rnamestore, -surplus, &surplus, 0, numrows-1 );
      checkCPXerror( err, "CPXgetrowname", "gutsOfCopy" );
      assert( surplus == 0 );
    }
  else
    {
      delete [] rname;
      rname = NULL;
    }

  err = CPXcopylpwnames( env_, getLpPtr(), 
			 numcols, numrows, objsen, 
			 const_cast<double *>(obj), 
			 const_cast<double *>(rhs), 
			 const_cast<char *>(sense),
			 const_cast<int *>(cols->vectorStarts()),
			 const_cast<int *>(cols->vectorLengths()),
			 const_cast<int *>(cols->indices()),
			 const_cast<double *>(cols->elements()),
			 const_cast<double *>(lb), 
			 const_cast<double *>(ub), 
			 rng, 
			 cname, rname);
  checkCPXerror( err, "CPXcopylpwnames", "gutsOfCopy" );
  
  if( rname != NULL )
    {
      delete [] rnamestore;
      delete [] rname;
    }
  if( cname != NULL )
    {
      delete [] cnamestore;
      delete [] cname;
    }
  delete [] rng;
#endif
 
}

//-------------------------------------------------------------------
void OsiCpxSolverInterface::gutsOfConstructor()
{  
#if 0
  // CPXcreateprob was moved to getLpPtr() method.
  int err;
  lp_ = CPXcreateprob( env_, &err, "OSI_CPLEX" );
  checkCPXerror( err, "CPXcreateprob", "gutsOfConstructor" );
  assert( lp_ != NULL );
#endif
}

//-------------------------------------------------------------------
void OsiCpxSolverInterface::gutsOfDestructor()
{  
  if ( lp_ != NULL )
    {
      int err = CPXfreeprob( env_, &lp_ );
      checkCPXerror( err, "CPXfreeprob", "gutsOfDestructor" );
      lp_=NULL;
      freeCachedData();
    }
  assert( lp_==NULL );
  assert( obj_==NULL );
  assert( collower_==NULL );
  assert( colupper_==NULL );
  assert( ctype_==NULL );
  assert( rowsense_==NULL );
  assert( rhs_==NULL );
  assert( rowrange_==NULL );
  assert( rowlower_==NULL );
  assert( rowupper_==NULL );
  assert( colsol_==NULL );
  assert( rowsol_==NULL );
  assert( redcost_==NULL );
  assert( rowact_==NULL );
  assert( matrixByRow_==NULL );
  assert( matrixByCol_==NULL );
}

//-------------------------------------------------------------------
/// free cached vectors

void OsiCpxSolverInterface::freeCachedColRim()
{
  freeCacheChar( ctype_ );
  freeCacheDouble( obj_ );  
  freeCacheDouble( collower_ ); 
  freeCacheDouble( colupper_ ); 
  assert( obj_==NULL );
  assert( collower_==NULL );
  assert( colupper_==NULL );
  assert( ctype_==NULL );
}

void OsiCpxSolverInterface::freeCachedRowRim()
{
  freeCacheChar( rowsense_ );
  freeCacheDouble( rhs_ );
  freeCacheDouble( rowrange_ );
  freeCacheDouble( rowlower_ );
  freeCacheDouble( rowupper_ );
  assert( rowsense_==NULL ); 
  assert( rhs_==NULL ); 
  assert( rowrange_==NULL ); 
  assert( rowlower_==NULL ); 
  assert( rowupper_==NULL );
 }

void OsiCpxSolverInterface::freeCachedResults()
{
  freeCacheDouble( colsol_ ); 
  freeCacheDouble( rowsol_ );
  freeCacheDouble( redcost_ );
  freeCacheDouble( rowact_ );
  assert( colsol_==NULL );
  assert( rowsol_==NULL );
  assert( redcost_==NULL );
  assert( rowact_==NULL );
}

void OsiCpxSolverInterface::freeCachedMatrix()
{
  freeCacheMatrix( matrixByRow_ );
  freeCacheMatrix( matrixByCol_ );
  assert( matrixByRow_==NULL ); 
  assert( matrixByCol_==NULL ); 
}


void OsiCpxSolverInterface::freeCachedData()
{
  freeCachedColRim();
  freeCachedRowRim();
  freeCachedResults();
  freeCachedMatrix();
}

#endif
