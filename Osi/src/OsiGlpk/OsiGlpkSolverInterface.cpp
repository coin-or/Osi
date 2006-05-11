//-----------------------------------------------------------------------------
// name:     OSI Interface for GLPK
//-----------------------------------------------------------------------------
// Copyright (C) 2001, 2002 Vivian De Smedt
// Copyright (C) 2002, 2003 Braden Hunsaker 
// Copyright (C) 2003, 2004 University of Pittsburgh
// Copyright (C) 2004 Joseph Young
//    University of Pittsburgh coding done by Brady Hunsaker
// All Rights Reserved.
//
// Comments:
//   
//    Areas that may need work in the code are marked with '???'.
//
//    As of version 4.7, GLPK problems can be of class LPX_LP or LPX_MIP.
// The difference is that LPX_MIP problems have a separate MIP data block, 
// and certain functions are only allowed for LPX_MIP problems.
//
//    In (much) earlier versions of GLPK, if an LPX_MIP problem was
// changed back into a LPX_LP problem, then the MIP data was lost,
// including which columns are integer.  However, LPX_MIP problems
// still had access to LP information (like lpx_get_status).
//
//    It appears that this behavior is no longer true in version 4.7.
// Therefore it may be worthwhile to adjust the interface to change
// the class of the problem as needed.  For the time being, however,
// all OSI problems are set to type LPX_MIP.  The only trick is
// differentiating when the user calls status routines like
// isProvenOptimal().  For now, we assume that the user is referring
// to the most recent solution call.  We add an extra variable,
// bbWasLast_, to the class to record which solution call was most
// recent (lp or mip).
//
//
// Possible areas of improvement
// -----------------------------
//
// Methods that are not implemented:
//
//  getPrimalRays, getDualRays
//
// Methods that are implemented, but do not do what you probably expect:
//
//  setColSolution, setRowPrice  
//
// Many methods have not been attempted to be improved (for speed) 
// to take advantage of how GLPK works.  The emphasis has been on
// making things functional and not too messy.  There's plenty of room
// for more efficient implementations.
//

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <string>
#include <iostream>
#include <stdio.h>

#include "CoinError.hpp"

#include "OsiGlpkSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

//-----------------------------------------------------------------------------

inline void checkGLPKerror( int err, std::string glpkfuncname, std::string osimethod )
{
	if( err != 0 )
	{
		char s[100];
		sprintf( s, "%s returned error %d", glpkfuncname.c_str(), err );
		std::cout << "ERROR: " << s << " (" << osimethod << " in OsiGlpkSolverInterface)" << std::endl;
		throw CoinError( s, osimethod.c_str(), "OsiGlpkSolverInterface" );
	}
}

//#############################################################################
// Solve methods
//#############################################################################

void OsiGlpkSolverInterface::initialSolve()
{
	// This method should solve the lp relaxation of a mip problem.
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );

	lpx_set_int_parm(model, LPX_K_MSGLEV, 1);  // suppress most output 
	lpx_set_int_parm(model, LPX_K_PRESOL, 1);  // turn on presolver
	int err = lpx_simplex( model );
	iter_used_ = lpx_get_int_parm(model, LPX_K_ITCNT);

	isIterationLimitReached_ = false;
	isAbandoned_ = false;
	isPrimInfeasible_ = false;
	isDualInfeasible_ = false;

	/* When the presolver is turned on, lpx_simplex() will not be able
	to tell whether the objective function has hit it's upper or lower
	limit. */

	isObjLowerLimitReached_ = false;
	isObjUpperLimitReached_ = false;

	switch ( err )
	{
	case LPX_E_ITLIM:
		isIterationLimitReached_ = true;
		break;

	// maybe more exit codes should count as abandoned
	case LPX_E_TMLIM:
	case LPX_E_FAULT:
	case LPX_E_SING:
		isAbandoned_ = true;
		break;

	case LPX_E_NOPFS:
		isPrimInfeasible_ = true;
		break;

	case LPX_E_NODFS:
		isDualInfeasible_ = true;
		break;
		
	}
	// Record that simplex was most recent
	bbWasLast_ = 0;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::resolve()
{
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );

	lpx_set_int_parm(model, LPX_K_MSGLEV, 1);  // suppress most output 
	lpx_set_int_parm(model, LPX_K_DUAL, 1); // Use dual simplex if dual feasible
	lpx_set_int_parm(model, LPX_K_PRESOL, 0);  // turn off presolver

	// lpx_simplex will use the current basis if possible
	int err = lpx_simplex( model );
	iter_used_ = lpx_get_int_parm(model, LPX_K_ITCNT);

	isIterationLimitReached_ = false;
	isAbandoned_ = false;
	isObjLowerLimitReached_ = false;
	isObjUpperLimitReached_ = false;
	isPrimInfeasible_ = false;
	isDualInfeasible_ = false;
	switch ( err )
	{
	case LPX_E_ITLIM:
		isIterationLimitReached_ = true;
		break;
	case LPX_E_OBJLL:
	  isObjLowerLimitReached_ = true;
	  break;
	case LPX_E_OBJUL:
	  isObjUpperLimitReached_ = true;
	  break;

		// maybe more exit codes should count as abandoned
	case LPX_E_FAULT:
	case LPX_E_SING:
		isAbandoned_ = true;
		break;
	case LPX_E_NOPFS:
		isPrimInfeasible_ = true;
		break;
	case LPX_E_NODFS:
		isDualInfeasible_ = true;
		break;
	}
	// Record that simplex was most recent
	bbWasLast_ = 0;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::branchAndBound()
{
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );
	if( lpx_get_num_int( model ) ) {

	        // Must have an LP solution before running lpx_integer
 	        if (lpx_get_status(model) != LPX_OPT)
	          initialSolve();
		// What if there's an error there?
		int err = lpx_integer( model );
		iter_used_ = lpx_get_int_parm(model, LPX_K_ITCNT);
		// Uncertain whether GLPK 4.7 keeps iteration count correctly
		// for MIPs ???

		isIterationLimitReached_ = false;
		isAbandoned_ = false;
		isPrimInfeasible_ = false;
		isDualInfeasible_ = false;
		switch( err )
		{
		case LPX_E_ITLIM:
			isIterationLimitReached_ = true;
			break;

		case LPX_E_FAULT:
		case LPX_E_SING:
			isAbandoned_ = true;
			break;
		}
		// Record that b&b was most recent
		bbWasLast_ = 1;
	}
	else {
		resolve();
	}
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiGlpkSolverInterface::setIntParam( OsiIntParam key, int value )
{
	bool retval = false;
	switch( key )
    {
    case OsiMaxNumIteration:
		if( value >= 0 )
		{
			maxIteration_ = value;
			lpx_set_int_parm(getMutableModelPtr(), LPX_K_ITLIM, 
					 value); 
			retval = true;
		}
		else
			retval = false;
		break;

    case OsiMaxNumIterationHotStart:
		if( value >= 0 )
		{
			hotStartMaxIteration_ = value;
			retval = true;
		}
		else
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
OsiGlpkSolverInterface::setDblParam( OsiDblParam key, double value )
{
  bool retval = false;
  switch ( key )
    {
    case OsiDualObjectiveLimit:
      // as of 4.7, GLPK only uses this if it does dual simplex
			dualObjectiveLimit_ = value;
			if (getObjSense()==1)  // minimization
			  lpx_set_real_parm( getMutableModelPtr(), 
					       LPX_K_OBJUL, value);
			else                   // maximization
			  lpx_set_real_parm( getMutableModelPtr(), 
					       LPX_K_OBJLL, value);
			retval = true; 
		break;

    case OsiPrimalObjectiveLimit:
      // as of 4.7, GLPK only uses this if it does dual simplex
			primalObjectiveLimit_ = value;
			if (getObjSense()==1) // minimization
			  lpx_set_real_parm( getMutableModelPtr(), 
					       LPX_K_OBJLL, value);
			else
			  lpx_set_real_parm( getMutableModelPtr(), 
					       LPX_K_OBJUL, value);
			retval = true; 
		break;

    case OsiDualTolerance:
		if( value >= 0 )
		{
		  dualTolerance_ = value;
		  lpx_set_real_parm( getMutableModelPtr(), 
				     LPX_K_TOLDJ, value);
		  retval = true;
		}
		else
		  retval = false;
		break;

    case OsiPrimalTolerance:
		if( value >= 0 )
		{
		  primalTolerance_ = value;
		  lpx_set_real_parm( getMutableModelPtr(), 
				     LPX_K_TOLBND, value);
		  retval = true; 
		}
		else
		  retval = false;
		break;

    case OsiObjOffset:
                lpx_set_obj_coef( getMutableModelPtr(), 0, value );
                retval = true;
                break;

    case OsiLastDblParam:
		retval = false;
		break;
    }
	return retval;
}

//-----------------------------------------------------------------------------

bool
OsiGlpkSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  switch (key) {
  case OsiProbName:
    lpx_set_prob_name( getMutableModelPtr(), const_cast<char *>(value.c_str()));
    return true;
  case OsiSolverName:
    return false;
  case OsiLastStrParam:
    return false;
  }
  return false;
}

//-----------------------------------------------------------------------------

bool
OsiGlpkSolverInterface::getIntParam( OsiIntParam key, int& value ) const
{
	bool retval = false;
	switch( key )
    {
    case OsiMaxNumIteration:
		value = maxIteration_;
		retval = true;
		break;

    case OsiMaxNumIterationHotStart:
		value = hotStartMaxIteration_;
		retval = true;
		break;

    case OsiLastIntParam:
		retval = false;
		break;
    }
	return retval;
}

//-----------------------------------------------------------------------------

bool
OsiGlpkSolverInterface::getDblParam( OsiDblParam key, double& value ) const
{
	bool retval = false;
	switch( key )
    {
    case OsiDualObjectiveLimit:
		value = dualObjectiveLimit_;
		retval = true;
		break;

    case OsiPrimalObjectiveLimit:
		value = primalObjectiveLimit_;
		retval = true;
		break;

    case OsiDualTolerance:
		value = dualTolerance_;
		retval = true;
		break;

    case OsiPrimalTolerance:
		value = primalTolerance_;
		retval = true;
		break;

    case OsiObjOffset:
                value = lpx_get_obj_coef(getMutableModelPtr(),0);
		retval = true;
		break;

    case OsiLastDblParam:
		retval = false;
		break;
    }
	return retval;
}


bool
OsiGlpkSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  //	bool retval = false;
  switch (key) {
  case OsiProbName:
    value = lpx_get_prob_name( getMutableModelPtr() );
    break;
  case OsiSolverName:
    value = "glpk";
    break;
  case OsiLastStrParam:
    return false;
  }
  return true;
}
//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiGlpkSolverInterface::isAbandoned() const
{
	return isAbandoned_;
}

bool OsiGlpkSolverInterface::isProvenOptimal() const
{
        LPX *model = getMutableModelPtr();

	if( bbWasLast_ == 0 )
	  {
	    int stat = lpx_get_status( model );
	    return stat == LPX_OPT;
	  }
	else
	  {
	    int stat = lpx_mip_status( model );
	    return stat == LPX_I_OPT;
	  }
}

bool OsiGlpkSolverInterface::isProvenPrimalInfeasible() const
{
        LPX *model = getMutableModelPtr();

	if(isPrimInfeasible_==true)
		return true;

	if( bbWasLast_ == 0 )
		return lpx_get_prim_stat( model ) == LPX_P_NOFEAS;
	else
		return lpx_mip_status( model ) == LPX_NOFEAS;
}

bool OsiGlpkSolverInterface::isProvenDualInfeasible() const
{
        LPX *model = getMutableModelPtr();
	
	if(isDualInfeasible_==true)
		return true;

	if( bbWasLast_ == 0 )
		return lpx_get_dual_stat( model ) == LPX_D_NOFEAS;
	else
	  // Not sure what to do for MIPs;  does it just mean unbounded?
	  // ??? for now, return false
	        return false;
}

bool OsiGlpkSolverInterface::isPrimalObjectiveLimitReached() const
{
  if (getObjSense()==1)  // minimization
    return isObjLowerLimitReached_;
  else  // maximization
    return isObjUpperLimitReached_;
}

bool OsiGlpkSolverInterface::isDualObjectiveLimitReached() const
{
  if (getObjSense()==1)  // minimization
    return isObjUpperLimitReached_;
  else  // maximization
    return isObjLowerLimitReached_;
}

bool OsiGlpkSolverInterface::isIterationLimitReached() const
{
	return isIterationLimitReached_;
}

//#############################################################################
// WarmStart related methods
//#############################################################################

CoinWarmStart* OsiGlpkSolverInterface::getWarmStart() const
{
	CoinWarmStartBasis *ws = NULL;
	LPX *model = getMutableModelPtr();

	ws = new CoinWarmStartBasis();

	int numcols = getNumCols();
	int numrows = getNumRows();
	ws->setSize( numcols, numrows );

	int i;
	for( i = 0; i < numrows; i++ )
	{
		int stat;
		double val;
		double dualVal;
		stat=lpx_get_row_stat(model,i+1);
		val=lpx_get_row_prim(model,i+1);
		dualVal=lpx_get_row_dual(model,i+1);
		switch( stat ) {
		case LPX_BS:
			ws->setArtifStatus( i, CoinWarmStartBasis::basic );
			break;

		case LPX_NS: // ??? I'am not completly sure it's the best interpretation.
		case LPX_NL:
			ws->setArtifStatus( i, CoinWarmStartBasis::atLowerBound );
			break;

		case LPX_NU:
			ws->setArtifStatus( i, CoinWarmStartBasis::atUpperBound );
			break;

		case LPX_NF:
			ws->setArtifStatus( i, CoinWarmStartBasis::isFree );
			break;

		default:
			assert( false );
			break;
		}
		//		ws->setArtifValue( i, val);
		//		ws->setArtifDualValue( i, dualVal );
	}

	int j;
	for( j = 0; j < numcols; j++ )
	{
		int stat;
		//		double val;
		//		double dualVal;
	        stat=lpx_get_col_stat(model,j+1);
		switch( stat ) {
		case LPX_BS:
			ws->setStructStatus( j, CoinWarmStartBasis::basic );
			break;

		case LPX_NS: // ??? I'am not completly sure it's the best interpretation.
		case LPX_NL:
			ws->setStructStatus( j, CoinWarmStartBasis::atLowerBound );
			break;

		case LPX_NU:
			ws->setStructStatus( j, CoinWarmStartBasis::atUpperBound );
			break;

		case LPX_NF:
			ws->setStructStatus( j, CoinWarmStartBasis::isFree );
			break;

		default:
			assert( false );
			break;
		}
		//		ws->setStructValue( j, val );
		//		ws->setStructDualValue( j, dualVal );
	}

	return ws;
}

//-----------------------------------------------------------------------------

bool OsiGlpkSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
	
	const CoinWarmStartBasis *ws = 
	  dynamic_cast<const CoinWarmStartBasis *>(warmstart);
	LPX *model = getMutableModelPtr();
	
	if( !ws )
		return false;
	
	int numcols = ws->getNumStructural();
	int numrows = ws->getNumArtificial();
	
	if( numcols != getNumCols() || numrows != getNumRows() )
		return false;
	
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );
	
	int i;
	for( i = 0; i < numrows; i++)
	{
		int stat;
		switch( ws->getArtifStatus(i) )
		{
		case CoinWarmStartBasis::basic:
			stat = LPX_BS;
			break;

		case CoinWarmStartBasis::atLowerBound:
			stat = LPX_NL;
			break;

		case CoinWarmStartBasis::atUpperBound:
			stat = LPX_NU;
			break;

		case CoinWarmStartBasis::isFree:
			stat = LPX_NF;
			break;

		default:  // unknown row status
			assert( false );
			return false;
		}		
		lpx_set_row_stat( model, i+1, stat );
	}
	
	int j;
	for( j = 0; j < numcols; j++)
	{
		int stat;
		switch( ws->getStructStatus( j ) )
		{
		case CoinWarmStartBasis::basic:
			stat = LPX_BS;
			break;

		case CoinWarmStartBasis::atLowerBound:
			stat = LPX_NL;
			break;

		case CoinWarmStartBasis::atUpperBound:
			stat = LPX_NU;
			break;

		case CoinWarmStartBasis::isFree:
			stat = LPX_NF;
			break;

		default:  // unknown col status
			assert( false );
			return false;
		}
		lpx_set_col_stat( model, j+1, stat );
	}
	
	return true;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiGlpkSolverInterface::markHotStart()
{
        LPX *model = getMutableModelPtr();
	int numcols, numrows;

	numcols = getNumCols();
	numrows = getNumRows();
	if( numcols > hotStartCStatSize_ )
    {
		delete[] hotStartCStat_;
		delete [] hotStartCVal_;
		delete [] hotStartCDualVal_;
		hotStartCStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numcols ) ); // get some extra space for future hot starts
		hotStartCStat_ = new int[hotStartCStatSize_];
		hotStartCVal_ = new double[hotStartCStatSize_];
		hotStartCDualVal_ = new double[hotStartCStatSize_];
    }
	int j;
	for( j = 0; j < numcols; j++ ) {
		int stat;
		double val;
		double dualVal;
		stat=lpx_get_col_stat(model,j);
		val=lpx_get_col_prim(model,j);
		dualVal=lpx_get_col_dual(model,j);
		hotStartCStat_[j] = stat;
		hotStartCVal_[j] = val;
		hotStartCDualVal_[j] = dualVal;
	}

	if( numrows > hotStartRStatSize_ )
    {
		delete [] hotStartRStat_;
		delete [] hotStartRVal_;
		delete [] hotStartRDualVal_;
		hotStartRStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numrows ) ); // get some extra space for future hot starts
		hotStartRStat_ = new int[hotStartRStatSize_];
		hotStartRVal_ = new double[hotStartRStatSize_];
		hotStartRDualVal_ = new double[hotStartRStatSize_];
    }
	int i;
	for( i = 0; i < numrows; i++ ) {
		int stat;
		double val;
		double dualVal;
		stat=lpx_get_row_stat(model,i+1);
		val=lpx_get_row_prim(model,i+1);
		dualVal=lpx_get_row_dual(model,i+1);
		hotStartRStat_[i] = stat;
		hotStartRVal_[i] = val;
		hotStartRDualVal_[i] = dualVal;
	}
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::solveFromHotStart()
{
        LPX *model = getMutableModelPtr();
	int numcols, numrows;

	numcols = getNumCols();
	numrows = getNumRows();

	assert( numcols <= hotStartCStatSize_ );
	assert( numrows <= hotStartRStatSize_ );

	int j;
	for( j = 0; j < numcols; j++ ) {
	  //	glp_put_col_soln( lp_, j+1, hotStartCStat_[j], hotStartCVal_[j], hotStartCDualVal_[j] );
	  // GLPK 4.7 doesn't have a way to set the values--just the status
	  lpx_set_col_stat( model, j+1, hotStartCStat_[j]);
	}
	int i;
	for( i = 0; i < numrows; i++ ) {
	  //	glp_put_row_soln( lp_, i+1, hotStartRStat_[i], hotStartRVal_[i], hotStartRDualVal_[i] );
	  lpx_set_row_stat( model, i+1, hotStartRStat_[i]);
	}

	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );

	int maxIteration = maxIteration_;
	maxIteration_ = hotStartMaxIteration_;
	resolve();
	maxIteration_ = maxIteration;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::unmarkHotStart()
{
	// ??? be lazy with deallocating memory and do nothing here, deallocate memory in the destructor.
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//-----------------------------------------------------------------------------
// Get number of rows, columns, elements, ...
//-----------------------------------------------------------------------------
int OsiGlpkSolverInterface::getNumCols() const
{
	return lpx_get_num_cols( getMutableModelPtr() );
}

int OsiGlpkSolverInterface::getNumRows() const
{
	return lpx_get_num_rows( getMutableModelPtr() );
}

int OsiGlpkSolverInterface::getNumElements() const
{
	return lpx_get_num_nz( getMutableModelPtr() );
}

//-----------------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getColLower() const
{
        LPX *model = getMutableModelPtr();

	if( collower_ == NULL )
	{
		assert( colupper_ == NULL );
		int numcols = getNumCols();
		if( numcols > 0 )
		{
			collower_ = new double[numcols];
			colupper_ = new double[numcols];
		}

		double inf = getInfinity();

		int i;
		for( i = 0; i < numcols; i++)
		{
			int type;
			double lb;
			double ub;
			type=lpx_get_col_type(model,i+1);
			lb=lpx_get_col_lb(model,i+1);
			ub=lpx_get_col_ub(model,i+1);
			switch ( type )
			{
			case LPX_FR:
				lb = -inf;
				ub = inf;
				break;

			case LPX_LO:
				ub = inf;
				break;

			case LPX_UP:
				lb = -inf;
				break;

			case LPX_FX:
			case LPX_DB:
				break;

			default:
				assert( false );
				break;
			}
			collower_[i] = lb;
			colupper_[i] = ub;
		}
	}
	return collower_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getColUpper() const
{
	if( colupper_ == NULL )
	{
		getColLower();
		//assert( colupper_ != NULL ); // Could be null if no columns.
	}
	return colupper_;
}

//-----------------------------------------------------------------------------

const char * OsiGlpkSolverInterface::getRowSense() const
{
	// Could be in OsiSolverInterfaceImpl.
	if( rowsense_ == NULL )
	{
		assert( rhs_ == NULL && rowrange_ == NULL );

		int numrows = getNumRows();
		if( numrows > 0 )
		{
			rowsense_ = new char[numrows];
			rhs_ = new double[numrows];
			rowrange_ = new double[numrows];
		}

		const double *rowlower = getRowLower();
		const double *rowupper = getRowUpper();
		int i;
		for( i = 0; i < numrows; i++ )
		{
			char sense;
			double right;
			double range;
			convertBoundToSense( rowlower[i], rowupper[i], sense, right, range );
			rowsense_[i] = sense;
			rhs_[i] = right;
			rowrange_[i] = range;
		}
	}
	return rowsense_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRightHandSide() const
{
	if( rhs_ == NULL )
	{
		getRowSense();
		//assert( rhs_ != NULL ); // Could be null if no rows.
	}
	return rhs_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRowRange() const
{
	if( rowrange_ == NULL )
	{
		getRowSense();
		//assert( rowrange_ != NULL ); // Could be null if no rows.
	}
	return rowrange_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRowLower() const
{
        LPX *model = getMutableModelPtr();

	if( rowlower_ == NULL )
	{
		assert( rowupper_ == NULL );
		int numrows = getNumRows();
		if( numrows > 0 )
		{
			rowlower_ = new double[numrows];
			rowupper_ = new double[numrows];
		}
		int i;
		for( i = 0; i < numrows; i++ )
		{
			double inf = getInfinity();
			int type;
			double lb;
			double ub;
			type=lpx_get_row_type(model,i+1);
			lb=lpx_get_row_lb(model,i+1);
			ub=lpx_get_row_ub(model,i+1);
			switch( type )
			{
			case LPX_FR:
				lb = -inf;
				ub = inf;
				break;

			case LPX_LO:
				ub = inf;
				break;

			case LPX_UP:
				lb = -inf;
				break;

			case LPX_DB:
			case LPX_FX:
				break;

			default:
				assert( false );
				break;
			}
			rowlower_[i] = lb;
			rowupper_[i] = ub;
		}
	}
	return rowlower_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRowUpper() const
{
	if( rowupper_ == NULL )
	{
		getRowLower();
		//assert( rowupper_ != NULL ); // Could be null if no rows.
	}
	return rowupper_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getObjCoefficients() const
{
	if( obj_ == NULL )
	{
	        LPX *model = getMutableModelPtr();

		int numcols = getNumCols();
		if( numcols > 0 )
		{
			obj_ = new double[numcols];
		}
		int i;
		for( i = 0; i < numcols; i++ )
		{
			obj_[i] = lpx_get_obj_coef( model, i + 1);
		}
	}
	return obj_;
}

//-----------------------------------------------------------------------------

double OsiGlpkSolverInterface::getObjSense() const
{
        LPX *model = getMutableModelPtr();

	if( lpx_get_obj_dir( model ) == LPX_MIN )
		// Minimization.
		return +1.0;
	else if( lpx_get_obj_dir( model ) == LPX_MAX )
		// Maximization.
		return -1.0;
	else
		assert( false );
	return 0;
}

//-----------------------------------------------------------------------------
// Return information on integrality
//-----------------------------------------------------------------------------

bool OsiGlpkSolverInterface::isContinuous( int colNumber ) const
{
  return lpx_get_col_kind( getMutableModelPtr(), colNumber+1 ) == LPX_CV;
}

//-----------------------------------------------------------------------------
// Row and column copies of the matrix ...
//-----------------------------------------------------------------------------

const CoinPackedMatrix * OsiGlpkSolverInterface::getMatrixByRow() const
{
	if( matrixByRow_ == NULL )
	{
	        LPX *model = getMutableModelPtr();

		matrixByRow_ = new CoinPackedMatrix();
		matrixByRow_->transpose();  // converts to row-order
		matrixByRow_->setDimensions( 0, getNumCols() );

		int numcols = getNumCols();
		int *colind = new int[numcols+1];
		double *colelem = new double[numcols+1];
		int i;
		for( i = 0; i < getNumRows(); i++ )
		{
			int colsize = lpx_get_mat_row( model, i+1, colind, colelem);
			int j;
			for( j = 1; j <= colsize; j++ )
			{
				--colind[j];
			}

			// Note:  lpx_get_mat_row apparently may return the
			// elements in decreasing order.  This differs from
			// people's standard expectations but is not an error.

			matrixByRow_->appendRow( colsize, colind+1, colelem+1 );
		}
		delete [] colind;
		delete [] colelem;
		if( numcols )
			matrixByRow_->removeGaps();
	}
	return matrixByRow_;
}

//-----------------------------------------------------------------------------

const CoinPackedMatrix * OsiGlpkSolverInterface::getMatrixByCol() const
{
	if( matrixByCol_ == NULL )
	{
   	        LPX *model = getMutableModelPtr();

		matrixByCol_ = new CoinPackedMatrix();
		matrixByCol_->setDimensions( getNumRows(), 0 );

		int numrows = getNumRows();
		int *rowind = new int[numrows+1];
		double *rowelem = new double[numrows+1];
		int j;
		for( j = 0; j < getNumCols(); j++ )
		{
			int rowsize = lpx_get_mat_col( model, j+1, rowind, rowelem);
			int i;
			for( i = 1; i <= rowsize; i++ )
			{
				--rowind[i];
			}
			matrixByCol_->appendCol( rowsize, rowind+1, rowelem+1 );
		}
		delete [] rowind;
		delete [] rowelem;
		if( numrows )
			matrixByCol_->removeGaps();
	}
	return matrixByCol_;
}

//-----------------------------------------------------------------------------
// Get solver's value for infinity
//-----------------------------------------------------------------------------
double OsiGlpkSolverInterface::getInfinity() const
{
  	return 1E+300;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

const double * OsiGlpkSolverInterface::getColSolution() const
{
	if( !colsol_ )
	{
 	        LPX *model = getMutableModelPtr();

		int numcols = getNumCols();
		if( numcols > 0 )
		{
			colsol_ = new double[numcols];
			redcost_ = new double[numcols];
		}

		// Function calls differ if last solve was b&b or simplex
		if (bbWasLast_ == 0)  // simplex/continuous solution
		  {
		    // Note: if the problem has not been solved, GLPK
		    //  returns 0, but OSI requires something
		    //  within bounds
		    if (lpx_get_status( model ) == LPX_UNDEF)
		      {
			// Make sure bounds are in cache.
			// this really handles upper bds too
			getColLower();  
			// now collower_ and colupper_ hold the bounds
			int j;
			for( j = 0; j < numcols; j++ )
			  {
			    colsol_[j] = (collower_[j] > 0) ? collower_[j]:0;
			    if (colupper_[j] < colsol_[j]) 
			      colsol_[j] = colupper_[j];
			  }
		      }
		    else  // this is what will usually be executed
		      {
			int j;
			for( j = 0; j < numcols; j++ )
			  {
			    // int status;
			    // status=lpx_get_col_stat(model,j+1);
			    colsol_[j] = lpx_get_col_prim(model,j+1);
			    redcost_[j] = lpx_get_col_dual(model,j+1);
			  }
		      }
		  }
		else   // b&b/integer solution
		  {
		    int j;
		    for( j = 0; j < numcols; j++ )
		      colsol_[j] = lpx_mip_col_val( model, j+1);
		  }
	}
	return colsol_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRowPrice() const
{
	if ( !rowsol_ )
		getRowActivity();
	return rowsol_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getReducedCost() const
{
	if ( !redcost_ )
		getColSolution();
	return redcost_;
}

//-----------------------------------------------------------------------------

const double * OsiGlpkSolverInterface::getRowActivity() const
{
#if 1
	if( rowact_ == NULL )
	{

 	        LPX *model = getMutableModelPtr();

		int numrows = getNumRows();
		if( numrows > 0 )
		{
			rowact_ = new double[numrows];
			rowsol_ = new double[numrows];
		}

		// Function calls differ if last solve was b&b or simplex
		if (bbWasLast_ == 0)  // simplex/continuous solution
		  {
		    int i;
		    for( i = 0; i < numrows; i++ )
		      {
			//int status;
			//status=lpx_get_row_stat(model,i+1);
			rowact_[i] = lpx_get_row_prim(model,i+1);
			rowsol_[i] = lpx_get_row_dual(model,i+1);
		      }
		  }
		else  // last call was b&b mip solution
		  {
		    int i;
		    for( i = 0; i < numrows; i++ )
		      {
			rowact_[i] = lpx_mip_row_val( model, i+1);
			rowsol_[i] = 0;  // no dual values for mip
		      }
		  }
	}
	return rowact_;
#else
	// ??? Consider removing this unused code.
	// Could be in OsiSolverInterfaceImpl.
	if( !rowact_ )
	{
		const double *colsol = getColSolution();
		const CoinPackedMatrix * matrix = getMatrixByRow();

		int numrows = getNumRows();
		rowact_ = new double[numrows];

		assert( numrows == matrix->getNumRows() );

		int i;
		for( i = 0; i < matrix->getNumRows(); i++ )
		{
			const CoinShallowPackedVector row = matrix->getVector(i);
			const double *elements = row.getElements();
			const int *indices = row.getIndices();

			rowact_[i] = 0;
			int j;
			for( j = 0; j < row.getNumElements(); j++ )
			{
				rowact_[i] += colsol[indices[j]] * elements[j];
			}
		}
	}
	return rowact_;
#endif
}

//-----------------------------------------------------------------------------

double OsiGlpkSolverInterface::getObjValue() const
{
  return OsiSolverInterface::getObjValue();
}

//-----------------------------------------------------------------------------

int OsiGlpkSolverInterface::getIterationCount() const
{
	return iter_used_;
}

//-----------------------------------------------------------------------------

std::vector<double*> OsiGlpkSolverInterface::getDualRays(int maxNumRays) const
{
	// ??? not yet implemented.
	throw CoinError("method is not yet implemented", "getDualRays", "OsiGlpkSolverInterface");
	return std::vector<double*>();
}

//-----------------------------------------------------------------------------

std::vector<double*> OsiGlpkSolverInterface::getPrimalRays(int maxNumRays) const
{
	// ??? not yet implemented.
	throw CoinError("method is not yet implemented", "getPrimalRays", "OsiGlpkSolverInterface");
	return std::vector<double*>();
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

void OsiGlpkSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );
	lpx_set_obj_coef( getMutableModelPtr(), elementIndex+1, elementValue );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setColLower(int elementIndex, double elementValue)
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );

	double inf = getInfinity();

	int type;
	double lb;
	double ub;

	type=lpx_get_col_type(getMutableModelPtr(),elementIndex+1);
	ub=lpx_get_col_ub(getMutableModelPtr(),elementIndex+1);
	lb = elementValue;
	switch( type )
	{
	case LPX_UP:
	case LPX_DB:
	case LPX_FX:
		break;

	case LPX_FR:
	case LPX_LO:
		ub = inf;
		break;

	default:
		assert( false );
	}
	setColBounds( elementIndex, lb, ub );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setColUpper(int elementIndex, double elementValue)
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );

	double inf = getInfinity();

	int type;
	double lb;
	double ub;

	type=lpx_get_col_type(getMutableModelPtr(),elementIndex+1);
	lb=lpx_get_col_lb(getMutableModelPtr(),elementIndex+1);
	ub = elementValue;
	switch( type )
	{
	case LPX_LO:
	case LPX_DB:
	case LPX_FX:
		break;

	case LPX_FR:
	case LPX_UP:
		lb = -inf;
		break;

	default:
		assert( false );
	}
	setColBounds( elementIndex, lb, ub );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setColBounds( int elementIndex, double lower, double upper )
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );

	double inf = getInfinity();

	int type;
	if( lower == upper )
		type = LPX_FX;
	else if( lower > -inf && upper < inf )
		type = LPX_DB;
	else if( lower > -inf )
		type = LPX_LO;
	else if( upper < inf)
		type = LPX_UP;
	else
		type = LPX_FR;

	lpx_set_col_bnds( getMutableModelPtr(), elementIndex+1, type, lower, upper );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setColSetBounds(const int* indexFirst,
					     const int* indexLast,
					     const double* boundList)
{
	OsiSolverInterface::setColSetBounds( indexFirst, indexLast, boundList );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setRowLower( int elementIndex, double elementValue )
{
	// Could be in OsiSolverInterfaceImpl.
	double inf = getInfinity();

	int type;
	double lb;
	double ub;

	type=lpx_get_row_type(getMutableModelPtr(),elementIndex+1);
	ub=lpx_get_row_ub(getMutableModelPtr(),elementIndex+1);
	lb = elementValue;
	switch( type )
	{
	case LPX_UP:
	case LPX_DB:
	case LPX_FX:
		break;

	case LPX_FR:
	case LPX_LO:
		ub = inf;
		break;

	default:
		assert( false );
	}
	setRowBounds( elementIndex, lb, ub );
}

//-----------------------------------------------------------------------------
void
OsiGlpkSolverInterface::setRowUpper( int elementIndex, double elementValue )
{
	// Could be in OsiSolverInterfaceImpl.
	double inf = getInfinity();

	int type;
	double lb;
	double ub;

	type=lpx_get_row_type(getMutableModelPtr(),elementIndex+1);
	lb=lpx_get_row_lb(getMutableModelPtr(),elementIndex+1);
	ub = elementValue;
	switch( type )
	{
	case LPX_LO:
	case LPX_DB:
	case LPX_FX:
		break;

	case LPX_FR:
	case LPX_UP:
		lb = -inf;
		break;

	default:
		assert( false );
	}
	setRowBounds( elementIndex, lb, ub );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_ROW );

	double inf = getInfinity();

	int type;
	if( lower == upper )
		type = LPX_FX;
	else if( lower > -inf && upper < inf )
		type = LPX_DB;
	else if( lower > -inf )
		type = LPX_LO;
	else if( upper < inf)
		type = LPX_UP;
	else
		type = LPX_FR;

	lpx_set_row_bnds( getMutableModelPtr(), elementIndex+1, type, lower, upper );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setRowType(int elementIndex, char sense, double rightHandSide,
								   double range)
{
	// Could be in OsiSolverInterfaceImpl.
	double lower;
	double upper;
    convertSenseToBound( sense, rightHandSide, range, lower, upper );
	setRowBounds( elementIndex, lower, upper );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setRowSetBounds(const int* indexFirst,
					     const int* indexLast,
					     const double* boundList)
{
	// Could be in OsiSolverInterface (should'nt be implemeted at here).
	OsiSolverInterface::setRowSetBounds( indexFirst, indexLast, boundList );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setRowSetTypes(const int* indexFirst,
				       const int* indexLast,
				       const char* senseList,
				       const double* rhsList,
				       const double* rangeList)
{
	// Could be in OsiSolverInterface (should'nt be implemeted at here).
	OsiSolverInterface::setRowSetTypes( indexFirst, indexLast, senseList, rhsList, rangeList );
}

//#############################################################################

void
OsiGlpkSolverInterface::setContinuous(int index)
{
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );
	lpx_set_col_kind( model, index+1, LPX_CV );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setInteger(int index)
{
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_COLUMN );
	lpx_set_col_kind( model, index+1, LPX_IV );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setContinuous(const int* indices, int len)
{
	// Could be in OsiSolverInterfaceImpl.
	int i;
	for( i = 0; i < len; i++ )
	{
		setContinuous( indices[i] );
	}
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::setInteger(const int* indices, int len)
{
	// Could be in OsiSolverInterfaceImpl.
	int i;
	for( i = 0; i < len; i++ )
	{
		setInteger( indices[i] );
	}
}

//#############################################################################

void OsiGlpkSolverInterface::setObjSense(double s)
{
	freeCachedData( OsiGlpkSolverInterface::FREECACHED_RESULTS );
	if( s == +1.0 )
		lpx_set_obj_dir( getMutableModelPtr(), LPX_MIN );
	else
		lpx_set_obj_dir( getMutableModelPtr(), LPX_MAX );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setColSolution(const double * cs)
{
  // You probably don't want to use this function.  You probably want
  // setWarmStart instead.
  // This implementation changes the cached information, 
  // BUT DOES NOT TELL GLPK about the changes.  In that sense, it's not
  // really useful.  It is added to conform to current OSI expectations.

  // Other results (such as row prices) might not make sense with this 
  // new solution, but we can't free all the results we have, since the 
  // row prices may have already been set with setRowPrice.
  if (cs == 0)
    delete [] colsol_;
  else
    {
      int nc = getNumCols();

      if (colsol_ == 0)
	colsol_ = new double[nc];

      // Copy in new col solution.
      CoinDisjointCopyN( cs, nc, colsol_ );
    }
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::setRowPrice(const double * rs)
{
  // You probably don't want to use this function.  You probably want
  // setWarmStart instead.
  // This implementation changes the cached information, 
  // BUT DOES NOT TELL GLPK about the changes.  In that sense, it's not
  // really useful.  It is added to conform to current OSI expectations.

  // Other results (such as column solutions) might not make sense with this 
  // new solution, but we can't free all the results we have, since the 
  // column solutions may have already been set with setColSolution.
  if (rs == 0)
    delete [] rowsol_;
  else
    {
      int nr = getNumRows();

      if (rowsol_ == 0)
	rowsol_ = new double[nr];

      // Copy in new col solution.
      CoinDisjointCopyN( rs, nr, rowsol_ );
    }
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################

void
OsiGlpkSolverInterface::addCol(const CoinPackedVectorBase& vec,
			       const double collb, const double colub,
			       const double obj)
{
  // Note: GLPK expects only non-zero coefficients will be given in 
  //   lpx_set_mat_col and will abort if there are any zeros.  So any
  //   zeros must be removed prior to calling lpx_set_mat_col.
        LPX *model = getMutableModelPtr();
	freeCachedData(OsiGlpkSolverInterface::KEEPCACHED_ROW);

	lpx_add_cols( model, 1 );
	int numcols = getNumCols();
	setColBounds( numcols-1, collb, colub );
	setObjCoeff( numcols-1, obj );
	int i;
	// For GLPK, we don't want the arrays to start at 0
	// We also need to weed out any 0.0 coefficients
	const int *indices = vec.getIndices();
	const double *elements = vec.getElements();
	int numrows = getNumRows();

	int *indices_adj = new int[1+vec.getNumElements()];
	double *elements_adj = new double[1+vec.getNumElements()];

	int count = 0;
	for( i = 0; i < vec.getNumElements(); i++ ) {
	  if (elements[i] != 0.0)
	    {
		if ( indices[i]+1 > numrows ) {
			lpx_add_rows( model, indices[i]+1 - numrows );
			numrows = indices[i]+1;
			// ??? could do this more efficiently with a single call based on the max
		}
		count++;
		// GLPK arrays start at 1
		indices_adj[count] = indices[i] + 1;
		elements_adj[count] = elements[i];
	    }
	}
	lpx_set_mat_col( model, numcols, count, indices_adj, elements_adj );
	delete [] indices_adj;
	delete [] elements_adj;
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::addCols(const int numcols,
				const CoinPackedVectorBase * const * cols,
				const double* collb, const double* colub,
				const double* obj)
{
  // ??? We could do this more efficiently now 
	// Could be in OsiSolverInterfaceImpl.
	int i;
	for( i = 0; i < numcols; ++i )
		addCol( *(cols[i]), collb[i], colub[i], obj[i] );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::deleteCols(const int num, const int * columnIndices)
{
	int *columnIndicesPlus1 = new int[num+1];
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::KEEPCACHED_ROW );

	for( int i = 0; i < num; i++ )
	{
		columnIndicesPlus1[i+1]=columnIndices[i]+1;
	}
	lpx_del_cols(model,num,columnIndicesPlus1);
	delete [] columnIndicesPlus1;
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::addRow(const CoinPackedVectorBase& vec,
							   const double rowlb, const double rowub)
{
  // Note: GLPK expects only non-zero coefficients will be given in 
  //   lpx_set_mat_row and will abort if there are any zeros.  So any
  //   zeros must be removed prior to calling lpx_set_mat_row.

        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::KEEPCACHED_COLUMN );

	lpx_add_rows( model, 1 );
	int numrows = getNumRows();
	setRowBounds( numrows-1, rowlb, rowub );
	int i;
	const int *indices = vec.getIndices();
	const double *elements = vec.getElements();
	int numcols = getNumCols();

	// For GLPK, we don't want the arrays to start at 0
	// Also, we need to weed out any 0.0 elements
	int *indices_adj = new int[1+vec.getNumElements()];
	double *elements_adj = new double[1+vec.getNumElements()];

	int count = 0;
	for( i = 0; i < vec.getNumElements(); i++ ) {
	  if ( elements[i] != 0.0 )
	    {
		if ( indices[i]+1 > numcols ) {
		  // ??? Could do this more efficiently with a single call
			lpx_add_cols( model, indices[i]+1 - numcols );
			numcols = indices[i]+1;
		}
		count++;
		elements_adj[count] = elements[i];
		indices_adj[count] = indices[i] + 1;
	    }
	}
	lpx_set_mat_row( model, numrows, count, indices_adj, elements_adj );
	delete [] indices_adj;
	delete [] elements_adj;
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::addRow(const CoinPackedVectorBase& vec,
			       const char rowsen, const double rowrhs,
			       const double rowrng)
{
	// Could be in OsiSolverInterfaceImpl.
	double lb;
	double ub;
	convertSenseToBound( rowsen, rowrhs, rowrng, lb, ub);
	addRow( vec, lb, ub );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::addRows(const int numrows, 
				const CoinPackedVectorBase * const * rows,
				const double* rowlb, const double* rowub)
{
  // ??? Could do this more efficiently now
	// Could be in OsiSolverInterfaceImpl.
	int i;
	for( i = 0; i < numrows; ++i )
		addRow( *(rows[i]), rowlb[i], rowub[i] );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::addRows(const int numrows,
				const CoinPackedVectorBase * const * rows,
				const char* rowsen, const double* rowrhs,
				const double* rowrng)
{
	// Could be in OsiSolverInterfaceImpl.
	int i;
	for( i = 0; i < numrows; ++i )
		addRow( *(rows[i]), rowsen[i], rowrhs[i], rowrng[i] );
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::deleteRows(const int num, const int * rowIndices)
{

	int *rowIndicesPlus1 = new int[num+1];
        LPX *model = getMutableModelPtr();
	freeCachedData( OsiGlpkSolverInterface::KEEPCACHED_COLUMN );

	for( int i = 0; i < num; i++ )
	{
		rowIndicesPlus1[i+1]=rowIndices[i]+1;
	}
	lpx_del_rows(model, num, rowIndicesPlus1);
	delete [] rowIndicesPlus1;
}

//#############################################################################
// Methods to input a problem
//#############################################################################

// Currently support for the default values (when passed a 0 pointer) is
// inefficient.  Should improve this. ???

void
OsiGlpkSolverInterface::loadProblem( const CoinPackedMatrix& matrix,
				     const double* collb, const double* colub,
				     const double* obj,
				     const double* rowlb, const double* rowub )
{
	// Could be in OsiSolverInterfaceImpl.
  // Actually, this could not.  I needed to add a glpk call below in case 
  // rows are present but have no nonzero coefficients.
  // JJF - If model exists then this adds! so I have added guts
  gutsOfDestructor();
  gutsOfConstructor();
  //freeCachedData( OsiGlpkSolverInterface::KEEPCACHED_NONE );

	double inf = getInfinity();

	if( matrix.isColOrdered() )
	{
		int i;
		for( i = 0; i < matrix.getNumCols(); i++ )
		{
			addCol( matrix.getVector(i), collb ? collb[i]:0.0, 
				colub ? colub[i]:inf, obj ? obj[i]:0.0 );
		}

		// Also make sure there are enough columns 
		if( matrix.getNumRows() > getNumRows())
		  lpx_add_rows( getMutableModelPtr(), matrix.getNumRows() - getNumRows() );
		int j;
		for( j = 0; j < matrix.getNumRows(); j++ )
		{
			setRowBounds( j, rowlb ? rowlb[j]:-inf, 
				      rowub ? rowub[j]:inf );
		}
	}
	else
	{
		int j;
		for( j = 0; j < matrix.getNumRows(); j++ )
		{
			addRow( matrix.getVector(j), rowlb ? rowlb[j]:-inf, 
				rowub ? rowub[j]:inf );
		}
		// Make sure there are enough columns
		if( matrix.getNumCols() > getNumCols())
		  lpx_add_cols( getMutableModelPtr(), matrix.getNumCols() - getNumCols() );

		int i;
		for( i = 0; i < matrix.getNumCols(); i++ )
		{
			setColBounds( i, collb ? collb[i]:0.0, 
				      colub ? colub[i]:inf );
			setObjCoeff( i, obj ? obj[i]:0.0 );
		}
	}
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::assignProblem( CoinPackedMatrix*& matrix,
				       double*& collb, double*& colub,
				       double*& obj,
				       double*& rowlb, double*& rowub )
{
	// Could be in OsiSolverInterfaceImpl.
	loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
	delete matrix;
	matrix = NULL;
	delete[] collb;
	collb = NULL;
	delete[] colub;
	colub = NULL;
	delete[] obj;
	obj = NULL;
	delete[] rowlb;
	rowlb = NULL;
	delete[] rowub;
	rowub = NULL;
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::loadProblem( const CoinPackedMatrix& matrix,
				     const double* collb, const double* colub,
				     const double* obj,
				     const char* rowsen, const double* rowrhs,
				     const double* rowrng )
{
	// Could be in OsiSolverInterfaceImpl.
	int numrows = matrix.getNumRows();
	double * rowlb = new double[numrows];
	double * rowub = new double[numrows];

	int i;
	for( i = numrows - 1; i >= 0; --i )
	{
		convertSenseToBound( rowsen ? rowsen[i]:'G', 
				     rowrhs ? rowrhs[i]:0.0, 
				     rowrng ? rowrng[i]:0.0, 
				     rowlb[i], rowub[i] );
	}

	loadProblem( matrix, collb, colub, obj, rowlb, rowub );
	delete [] rowlb;
	delete [] rowub;
}

//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::assignProblem( CoinPackedMatrix*& matrix,
				       double*& collb, double*& colub,
				       double*& obj,
				       char*& rowsen, double*& rowrhs,
				       double*& rowrng )
{
	// Could be in OsiSolverInterfaceImpl.
	loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );
	delete matrix;
	matrix = NULL;
	delete[] collb;
	collb = NULL;
	delete[] colub;
	colub = NULL;
	delete[] obj;
	obj = NULL;
	delete[] rowsen;
	rowsen = NULL;
	delete[] rowrhs;
	rowrhs = NULL;
	delete[] rowrng;
	rowrng = NULL;
}
//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  freeCachedData( OsiGlpkSolverInterface::KEEPCACHED_NONE );
  LPX *model = getMutableModelPtr();
  double inf = getInfinity();

  // Can't send 0 to lpx_add_xxx
  if (numcols > 0)
    lpx_add_cols( model, numcols );
  if (numrows > 0)
    lpx_add_rows( model, numrows );

  // How many elements?  Column-major, so indices of start are columns
  int numelem = start[ numcols ];
  //  int numelem = 0;
  //  while ( index[numelem] != 0 )
  //    numelem++;
  int *index_adj = new int[1+numelem];
  double *value_adj = new double[1+numelem];

  int i;
  for ( i=1; i <= numelem; i++)
    {
      index_adj[i] = index[i-1] + 1;
      value_adj[i] = value[i-1];
    }

  for( i = 0; i < numcols; i++ )
  {
	setColBounds( i, collb ? collb[i]:0.0, 
		    colub ? colub[i]:inf );
	lpx_set_mat_col( model, i+1, start[i+1]-start[i], &(index_adj[start[i]]), &(value_adj[start[i]]) );
    setObjCoeff( i, obj ? obj[i]:0.0 );
  }
  int j;
  for( j = 0; j < numrows; j++ )
  {
      setRowBounds( j, rowlb ? rowlb[j]:-inf, rowub ? rowub[j]:inf );
  }

  delete [] index_adj;
  delete [] value_adj;
  
}
//-----------------------------------------------------------------------------

void
OsiGlpkSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
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

   loadProblem( numcols, numrows, start, index, value, collb, colub, obj,
		rowlb, rowub );

   delete[] rowlb;
   delete[] rowub;
}

//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------

int OsiGlpkSolverInterface::readMps( const char * filename,
									 const char * extension )
{
#if 0
	std::string f( filename );
	std::string e( extension );
	std::string fullname = f + "." + e;
	lp_ = lpx_read_mps( const_cast<char*>( fullname.c_str() ) );
	assert( lp_ );
	return 0;  // Should actually check for errors here
#endif
  // just call base class method
  return OsiSolverInterface::readMps(filename,extension);
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::writeMps( const char * filename,
				       const char * extension,
				       double objSense ) const
{
	// Could be in OsiSolverInterfaceImpl.
#if 1
	std::string f( filename );
	std::string e( extension );
	std::string fullname = f + "." + e;
	lpx_write_mps( getMutableModelPtr(), const_cast<char*>( fullname.c_str() ));
#else
	// Fall back on native MPS writer.  
	// These few lines of code haven't been tested. 2004/10/15
	std::string f( filename );
	std::string e( extension );
	std::string fullname = f + "." + e;

	OsiSolverInterface::writeMpsNative(fullname.c_str(), 
					   NULL, NULL, 0, 2, objSense); 
#endif
}

//############################################################################
// GLPK-specific methods
//############################################################################

// Get a pointer to the instance
LPX * OsiGlpkSolverInterface::getModelPtr ()
{
  freeCachedResults();
  return lp_;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-----------------------------------------------------------------------------
// Default Constructor
//-----------------------------------------------------------------------------

OsiGlpkSolverInterface::OsiGlpkSolverInterface ()
: 
OsiSolverInterface()
{
	gutsOfConstructor();
}

//-----------------------------------------------------------------------------
// Clone
//-----------------------------------------------------------------------------

OsiSolverInterface * OsiGlpkSolverInterface::clone(bool copyData) const
{
  if (copyData)
    return( new OsiGlpkSolverInterface( *this ) );
  else
    return( new OsiGlpkSolverInterface() );
}

//-----------------------------------------------------------------------------
// Copy constructor
//-----------------------------------------------------------------------------

OsiGlpkSolverInterface::OsiGlpkSolverInterface( const OsiGlpkSolverInterface & source )
: 
OsiSolverInterface(source)
{
	gutsOfConstructor();
	gutsOfCopy( source );
}

//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

OsiGlpkSolverInterface::~OsiGlpkSolverInterface ()
{
	gutsOfDestructor();
}

// Resets 
// ??? look over this carefully to be sure it is correct
void OsiGlpkSolverInterface::reset()
{
   setInitialData();  // this is from the base class OsiSolverInterface
   gutsOfDestructor();

   bbWasLast_ = 0;
   
   maxIteration_ = INT_MAX;
   hotStartMaxIteration_ = 0;
   
   lp_ = lpx_create_prob();
   char name[] = "OSI_GLPK";
   lpx_set_prob_name( lp_, name );
   lpx_set_class( lp_, LPX_MIP );  
   // See note at top of file.   Use LPX_MIP even for LPs.
   assert( lp_ != NULL );

}


//-----------------------------------------------------------------------------
// Assignment operator
//-----------------------------------------------------------------------------

OsiGlpkSolverInterface& OsiGlpkSolverInterface::operator=( const OsiGlpkSolverInterface& rhs )
{
	if( this != &rhs )
	{
		OsiSolverInterface::operator=( rhs );
		gutsOfDestructor();
		gutsOfConstructor();
		if( rhs.lp_ != NULL )
			gutsOfCopy( rhs );
	}
	return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiGlpkSolverInterface::applyColCut( const OsiColCut & cc )
{
	// Could be in OsiSolverInterfaceImpl.
	const double * colLb = getColLower();
	const double * colUb = getColUpper();
	const CoinPackedVector & lbs = cc.lbs();
	const CoinPackedVector & ubs = cc.ubs();
	int i;
#if 0
        // replaced (JJF) because colLb and colUb are invalidated by sets
	for( i = 0; i < lbs.getNumElements(); ++i )
		if( lbs.getElements()[i] > colLb[lbs.getIndices()[i]] )
			setColLower( lbs.getIndices()[i], lbs.getElements()[i] );
		for( i = 0; i < ubs.getNumElements(); ++i )
			if( ubs.getElements()[i] < colUb[ubs.getIndices()[i]] )
				setColUpper( ubs.getIndices()[i], ubs.getElements()[i] );
#else
	double inf = getInfinity();

	int type;
        // lower bounds
	for( i = 0; i < lbs.getNumElements(); ++i ) {
          int column = lbs.getIndices()[i];
          double lower = lbs.getElements()[i];
          double upper = colUb[column];
          if( lower > colLb[column] ) {
            // update cached version as well
            collower_[column] = lower;
            if( lower == upper )
              type = LPX_FX;
            else if( lower > -inf && upper < inf )
              type = LPX_DB;
            else if( lower > -inf )
              type = LPX_LO;
            else if( upper < inf)
              type = LPX_UP;
            else
              type = LPX_FR;
            
            lpx_set_col_bnds( getMutableModelPtr(), column+1, type, lower, upper );
          }
        }
        // lower bounds
	for( i = 0; i < ubs.getNumElements(); ++i ) {
          int column = ubs.getIndices()[i];
          double upper = ubs.getElements()[i];
          double lower = colLb[column];
          if( upper < colUb[column] ) {
            // update cached version as well
            colupper_[column] = upper;
            if( lower == upper )
              type = LPX_FX;
            else if( lower > -inf && upper < inf )
              type = LPX_DB;
            else if( lower > -inf )
              type = LPX_LO;
            else if( upper < inf)
              type = LPX_UP;
            else
              type = LPX_FR;
            
            lpx_set_col_bnds( getMutableModelPtr(), column+1, type, lower, upper );
          }
        }
#endif
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
	// Could be in OsiSolverInterfaceImpl.
	addRow( rowCut.row(), rowCut.lb(), rowCut.ub() );
}

//#############################################################################
// Private methods (non-static and static) and static data
//#############################################################################

LPX * OsiGlpkSolverInterface::getMutableModelPtr( void ) const
{
  return lp_;
}

void OsiGlpkSolverInterface::gutsOfCopy( const OsiGlpkSolverInterface & source )
{
	// Set Objective Sense
	setObjSense( source.getObjSense() );

	// Set Rim and constraints
	const double* obj = source.getObjCoefficients();
	const double* rhs = source.getRightHandSide();
	const char* sense = source.getRowSense();
	const CoinPackedMatrix* cols = source.getMatrixByCol();
	const double* lb = source.getColLower();
	const double* ub = source.getColUpper();
	loadProblem( *cols, lb, ub, obj, sense, rhs, source.getRowRange() );

	bbWasLast_ = source.bbWasLast_;
	iter_used_ = source.iter_used_;

	// ??? Need to copy parameters somehow?
	

	int numcols = getNumCols();
	// Set MIP information
	int j;
	for( j = 0; j < numcols; j++ )
	{
		if( !source.isContinuous(j) )
		{
			setInteger(j);
		}
	}

	// Set Solution

	// In this case, it's easier to use GLPK's own functions than
	// to go through COIN/OSI interface.
	LPX *srcmodel = source.getMutableModelPtr();
	LPX *model = getMutableModelPtr();
	int tagx;
	for ( j = 1; j <= numcols; j++ )
	  {
	    tagx=lpx_get_col_stat(model,j);
	    lpx_set_col_stat( model, j, tagx );
	  }
	int numrows = getNumRows();
	for ( j = 1; j <= numrows; j++ )
	  {
	    tagx=lpx_get_row_stat(srcmodel,j);
	    lpx_set_row_stat( model, j, tagx );
	  }

	// In case the cache is different, we'll use setColSolution
	setColSolution(source.getColSolution());
	setRowPrice(source.getRowPrice());

	// Now we have GLPK construct the basis so it has solution values
	lpx_warm_up( model );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::gutsOfConstructor()
{
        bbWasLast_ = 0;
        iter_used_ = 0;
	obj_ = NULL;
	collower_ = NULL;
	colupper_ = NULL;
	ctype_ = NULL;
	rowsense_ = NULL;
	rhs_ = NULL;
	rowrange_ = NULL;
	rowlower_ = NULL;
	rowupper_ = NULL;
	colsol_ = NULL;
	rowsol_ = NULL;
	redcost_ = NULL;
	rowact_ = NULL;
	matrixByRow_ = NULL;
	matrixByCol_ = NULL;

	maxIteration_ = INT_MAX;
	hotStartMaxIteration_ = 0;

	hotStartCStat_ = NULL;
	hotStartCStatSize_ = 0;
	hotStartRStat_ = NULL;
	hotStartRStatSize_ = 0;

	isIterationLimitReached_ = false;
	isAbandoned_ = false;
	isPrimInfeasible_ = false;
	isDualInfeasible_ = false;

	lp_ = lpx_create_prob();
	char name[] = "OSI_GLPK";
	lpx_set_prob_name( lp_, name );
	// Make all problems MIPs.  See note at top of file.
	lpx_set_class( lp_, LPX_MIP );
	assert( lp_ != NULL );
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::gutsOfDestructor()
{
	if( lp_ != NULL )
	{
		lpx_delete_prob( lp_ );
		lp_=NULL;
		freeAllMemory();
	}
	assert( lp_ == NULL );
	assert( obj_ == NULL );
	assert( collower_ == NULL );
	assert( colupper_ == NULL );
	assert( ctype_ == NULL );
	assert( rowsense_ == NULL );
	assert( rhs_ == NULL );
	assert( rowrange_ == NULL );
	assert( rowlower_ == NULL );
	assert( rowupper_ == NULL );
	assert( colsol_ == NULL );
	assert( rowsol_ == NULL );
	assert( redcost_ == NULL );
	assert( rowact_ == NULL );
	assert( matrixByRow_ == NULL );
	assert( matrixByCol_ == NULL );
}

//-----------------------------------------------------------------------------
// free cached vectors
//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeCachedColRim()
{
	delete [] ctype_;
	delete [] obj_;
	delete [] collower_;
	delete [] colupper_;
	ctype_ = NULL;
	obj_ = NULL;
	collower_ = NULL;
	colupper_ = NULL;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeCachedRowRim()
{
	delete [] rowsense_;
	delete [] rhs_;
	delete [] rowrange_;
	delete [] rowlower_;
	delete [] rowupper_;
	rowsense_ = NULL;
	rhs_ = NULL;
	rowrange_ = NULL;
	rowlower_ = NULL;
	rowupper_ = NULL;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeCachedMatrix()
{
	delete matrixByRow_;
	delete matrixByCol_;
	matrixByRow_ = NULL;
	matrixByCol_ = NULL;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeCachedResults()
{
        iter_used_ = 0;
	isAbandoned_ = false;
	isIterationLimitReached_ = false;
	isPrimInfeasible_ = false;
	isDualInfeasible_ = false;
	delete [] colsol_;
	delete [] rowsol_;
	delete [] redcost_;
	delete [] rowact_;
	colsol_ = NULL;
	rowsol_ = NULL;
	redcost_ = NULL;
	rowact_ = NULL;
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeCachedData( int keepCached )
{
	if( !(keepCached & OsiGlpkSolverInterface::KEEPCACHED_COLUMN) )
		freeCachedColRim();
	if( !(keepCached & OsiGlpkSolverInterface::KEEPCACHED_ROW) )
		freeCachedRowRim();
	if( !(keepCached & OsiGlpkSolverInterface::KEEPCACHED_MATRIX) )
		freeCachedMatrix();
	if( !(keepCached & OsiGlpkSolverInterface::KEEPCACHED_RESULTS) )
		freeCachedResults();
}

//-----------------------------------------------------------------------------

void OsiGlpkSolverInterface::freeAllMemory()
{
	freeCachedData();
	delete[] hotStartCStat_;
	delete[] hotStartRStat_;
	hotStartCStat_ = NULL;
	hotStartCStatSize_ = 0;
	hotStartRStat_ = NULL;
	hotStartRStatSize_ = 0;
}
