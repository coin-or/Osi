//-----------------------------------------------------------------------------
// name:     OSI Interface for Gurobi
// template: OSI Cplex Interface written by T. Achterberg
// author:   Stefan Vigerske
//           Humboldt University Berlin
// comments: please scan this file for '???' and 'TODO' and read the comments
//-----------------------------------------------------------------------------
// Copyright (C) 2009 Humboldt University Berlin and others.
// All Rights Reserved.

// $Id$

#include <iostream>
#include <cassert>
#include <string>
#include <numeric>

#include "CoinPragma.hpp"
#include "CoinError.hpp"

#include "OsiGrbSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

extern "C" {
#include "gurobi_c.h"
}

#define GUROBI_CALL(m, x) do \
{ \
  int _retcode; \
  if( (_retcode = (x)) != 0 ) \
  { \
    char s[1001]; \
    sprintf( s, "%s:%d: Error <%d> from GUROBI function call: ", __FILE__, __LINE__, _retcode ); \
    if (OsiGrbSolverInterface::globalenv_) \
      strncat(s, GRBgeterrormsg(OsiGrbSolverInterface::globalenv_), 1000); \
    throw CoinError( s, m, "OsiGrbSolverInterface", __FILE__, __LINE__ ); \
  } \
} \
while( false )

// #define DEBUG 1

#ifdef DEBUG
#define debugMessage printf
#else
#define debugMessage if( false ) printf
#endif


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

inline void freeCacheMatrix( CoinPackedMatrix*& ptr )
{
  if( ptr != NULL )
    {
      delete ptr;
      ptr = NULL;
    }
}

void
OsiGrbSolverInterface::switchToLP( void )
{
  debugMessage("OsiGrbSolverInterface::switchToLP()\n");

  if( probtypemip_ )
  {
     GRBmodel* lp = getMutableLpPtr();
     int nc = getNumCols();
     
     char* contattr = new char[nc];
     CoinFillN(contattr, nc, 'C');
     
     GUROBI_CALL( "switchToLP", GRBsetcharattrarray(lp, GRB_CHAR_ATTR_VTYPE, 0, nc, contattr) );
     
     delete[] contattr;

     probtypemip_ = false;
  }
}

void
OsiGrbSolverInterface::switchToMIP( void )
{
  debugMessage("OsiGrbSolverInterface::switchToMIP()\n");

  if( !probtypemip_ )
  {
     GRBmodel* lp = getMutableLpPtr();
     int nc = getNumCols();

     assert(coltype_ != NULL);

     GUROBI_CALL( "switchToMIP", GRBsetcharattrarray(lp, GRB_CHAR_ATTR_VTYPE, 0, nc, coltype_) );

     probtypemip_ = true;
  }
}

void
OsiGrbSolverInterface::resizeColType( int minsize )
{
  debugMessage("OsiGrbSolverInterface::resizeColType()\n");

  if( minsize > coltypesize_ )
  {
     int newcoltypesize = 2*coltypesize_;
     if( minsize > newcoltypesize )
        newcoltypesize = minsize;
     char *newcoltype = new char[newcoltypesize];

     if( coltype_ != NULL )
     {
        CoinDisjointCopyN( coltype_, coltypesize_, newcoltype );
        delete[] coltype_;
     }
     coltype_ = newcoltype;
     coltypesize_ = newcoltypesize;
  }
  assert(minsize == 0 || coltype_ != NULL);
  assert(coltypesize_ >= minsize);
}

void
OsiGrbSolverInterface::freeColType()
{
  debugMessage("OsiGrbSolverInterface::freeColType()\n");

   if( coltypesize_ > 0 )
   {
      delete[] coltype_;
      coltype_ = NULL;
      coltypesize_ = 0;
   }
   assert(coltype_ == NULL);
}


//#############################################################################
// Solve methods
//#############################################################################

void OsiGrbSolverInterface::initialSolve()
{
  debugMessage("OsiGrbSolverInterface::initialSolve()\n");
  bool takeHint, gotHint;
  OsiHintStrength strength;

  switchToLP();
  
  GRBmodel* lp = getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS );

  /* set whether dual or primal */
  int algorithm = GRB_LPMETHOD_PRIMAL;
  gotHint = getHintParam(OsiDoDualInInitial,takeHint,strength);
  assert (gotHint);
  if (strength!=OsiHintIgnore)
  	algorithm = takeHint ? GRB_LPMETHOD_DUAL : GRB_LPMETHOD_PRIMAL;

	GUROBI_CALL( "initialSolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_LPMETHOD, algorithm) );

	/* set whether presolve or not */
  int presolve = GRB_PRESOLVE_AUTO;
  gotHint = (getHintParam(OsiDoPresolveInInitial,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
  	presolve = takeHint ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF;

  GUROBI_CALL( "initialSolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve) );

	/* set whether output or not */
  GUROBI_CALL( "initialSolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0)) );

	/* optimize */
  GUROBI_CALL( "initialSolve", GRBoptimize(lp) );
  
  /* reoptimize without presolve if status unclear */
  int stat;
  GUROBI_CALL( "initialSolve", GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &stat) );

  if (stat == GRB_INF_OR_UNBD && presolve != GRB_PRESOLVE_OFF) {
    GUROBI_CALL( "initialSolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
    GUROBI_CALL( "initialSolve", GRBoptimize(lp) );
    GUROBI_CALL( "initialSolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve) );
  }
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::resolve()
{
  debugMessage("OsiGrbSolverInterface::resolve()\n");
  bool takeHint, gotHint;
  OsiHintStrength strength;

  switchToLP();
  
  GRBmodel* lp = getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS );

  /* set whether primal or dual */
  int algorithm = GRB_LPMETHOD_DUAL;
  gotHint = getHintParam(OsiDoDualInResolve,takeHint,strength);
  assert (gotHint);
  if (strength != OsiHintIgnore)
  	algorithm = takeHint ? GRB_LPMETHOD_DUAL : GRB_LPMETHOD_PRIMAL;

  GUROBI_CALL( "resolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_LPMETHOD, algorithm) );

	/* set whether presolve or not */
  int presolve = GRB_PRESOLVE_OFF;
  gotHint = getHintParam(OsiDoPresolveInResolve,takeHint,strength);
  assert (gotHint);
  if (strength != OsiHintIgnore)
  	presolve = takeHint ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF;

  GUROBI_CALL( "resolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve) );

	/* set whether output or not */
  GUROBI_CALL( "resolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0)) );

	/* optimize */
  GUROBI_CALL( "resolve", GRBoptimize(lp) );

  /* reoptimize if status unclear */
  int stat;
  GUROBI_CALL( "resolve", GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &stat) );

  if (stat == GRB_INF_OR_UNBD && presolve != GRB_PRESOLVE_OFF) {
    GUROBI_CALL( "resolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF) );
    GUROBI_CALL( "resolve", GRBoptimize(lp) );
    GUROBI_CALL( "resolve", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve) );
  }
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::branchAndBound()
{
  debugMessage("OsiGrbSolverInterface::branchAndBound()\n");

  switchToMIP();

  GRBmodel* lp = getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS );

  GUROBI_CALL( "branchAndBound", GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0)) );

  GUROBI_CALL( "branchAndBound", GRBoptimize( lp ) );
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiGrbSolverInterface::setIntParam(OsiIntParam key, int value)
{
  debugMessage("OsiGrbSolverInterface::setIntParam(%d, %d)\n", key, value);

  bool retval = false;
  switch (key)
  {
    case OsiMaxNumIteration:
      GUROBI_CALL( "setIntParam", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, (double)value) );
      retval = true;
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
    default:
      retval = false ;
      break ;
  }
  return retval;
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::setDblParam(OsiDblParam key, double value)
{
  debugMessage("OsiGrbSolverInterface::setDblParam(%d, %g)\n", key, value);

  switch (key)
  {
//TODO
//  	case OsiDualObjectiveLimit:
//  		break;
  	case OsiPrimalObjectiveLimit:
  	  GUROBI_CALL( "setDblParam", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_CUTOFF, value) );
      return true;
  	case OsiDualTolerance:
  	  GUROBI_CALL( "setDblParam", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_OPTIMALITYTOL, value) );
  	  return true;
  	case OsiPrimalTolerance:
  	  GUROBI_CALL( "setDblParam", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_FEASIBILITYTOL, value) );
  	  return true;
  	case OsiObjOffset:
  		return OsiSolverInterface::setDblParam(key,value);
  	case OsiLastDblParam:
  	default:
  		return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  debugMessage("OsiGrbSolverInterface::setStrParam(%d, %s)\n", key, value.c_str());

  switch (key) {
    case OsiProbName:
      return OsiSolverInterface::setStrParam(key,value);

    case OsiSolverName:
    case OsiLastStrParam:
    default:
      return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  debugMessage("OsiGrbSolverInterface::getIntParam(%d)\n", key);

  switch (key)
  {
    case OsiMaxNumIteration:
    	double dblval;
      GUROBI_CALL( "getIntParam", GRBupdatemodel(getMutableLpPtr()) );
      GUROBI_CALL( "getIntParam", GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, &dblval) );
    	value = (int) dblval;
      return true;

    case OsiMaxNumIterationHotStart:
      value = hotStartMaxIteration_;
      return true;

    case OsiLastIntParam:
    default:
      return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  debugMessage("OsiGrbSolverInterface::getDblParam(%d)\n", key);

  GUROBI_CALL( "getDblParam", GRBupdatemodel(getMutableLpPtr()) );

  switch (key) 
  {
//TODO
//    	case OsiDualObjectiveLimit:
//    		break;
   	case OsiPrimalObjectiveLimit:
   	  GUROBI_CALL( "getDblParam", GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_CUTOFF, &value) );
   	  return true;
  	case OsiDualTolerance:
  	  GUROBI_CALL( "getDblParam", GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_OPTIMALITYTOL, &value) );
  		return true;
  	case OsiPrimalTolerance:
  	  GUROBI_CALL( "getDblParam", GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_FEASIBILITYTOL, &value) );
  	  return true;
    case OsiObjOffset:
      return OsiSolverInterface::getDblParam(key, value);
    case OsiLastDblParam:
    default:
      return false;
  }
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  debugMessage("OsiGrbSolverInterface::getStrParam(%d)\n", key);

  switch (key) {
  	case OsiProbName:
  		return OsiSolverInterface::getStrParam(key, value);

  	case OsiSolverName:
  		value = "gurobi";
  		return true;

  	case OsiLastStrParam:
  	default:
  		return false;
  }
}

//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiGrbSolverInterface::isAbandoned() const
{
  debugMessage("OsiGrbSolverInterface::isAbandoned()\n");

  GUROBI_CALL( "isAbandoned", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isAbandoned", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );
	
	return (
			stat == GRB_LOADED ||
			stat == GRB_NUMERIC || 
			stat == GRB_INTERRUPTED
		);
}

bool OsiGrbSolverInterface::isProvenOptimal() const
{
  debugMessage("OsiGrbSolverInterface::isProvenOptimal()\n");

  GUROBI_CALL( "isProvenOptimal", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isProvenOptimal", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );

	return (stat == GRB_OPTIMAL);
}

bool OsiGrbSolverInterface::isProvenPrimalInfeasible() const
{
  debugMessage("OsiGrbSolverInterface::isProvenPrimalInfeasible()\n");

  GUROBI_CALL( "isProvenPrimalInfeasible", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isProvenPrimalInfeasible", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );

	return (stat == GRB_INFEASIBLE);
}

bool OsiGrbSolverInterface::isProvenDualInfeasible() const
{
  debugMessage("OsiGrbSolverInterface::isProvenDualInfeasible()\n");

  GUROBI_CALL( "isProvenDualInfeasible", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isProvenDualInfeasible", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );

	return (stat == GRB_UNBOUNDED);
}

bool OsiGrbSolverInterface::isPrimalObjectiveLimitReached() const
{
  debugMessage("OsiGrbSolverInterface::isPrimalObjectiveLimitReached()\n");

  GUROBI_CALL( "isPrimalObjectiveLimitReached", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isPrimalObjectiveLimitReached", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );

	return (stat == GRB_CUTOFF);
}

bool OsiGrbSolverInterface::isDualObjectiveLimitReached() const
{
  debugMessage("OsiGrbSolverInterface::isDualObjectiveLimitReached()\n");

  return false;
}

bool OsiGrbSolverInterface::isIterationLimitReached() const
{
  debugMessage("OsiGrbSolverInterface::isIterationLimitReached()\n");

  GUROBI_CALL( "isIterationLimitReached", GRBupdatemodel(getMutableLpPtr()) );

  int stat;
  GUROBI_CALL( "isIterationLimitReached", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat) );

	return (stat == GRB_ITERATION_LIMIT);
}

//#############################################################################
// WarmStart related methods
//#############################################################################

CoinWarmStart* OsiGrbSolverInterface::getEmptyWarmStart () const
{ return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartBasis())) ; }

CoinWarmStart* OsiGrbSolverInterface::getWarmStart() const
{
  debugMessage("OsiGrbSolverInterface::getWarmStart()\n");

  CoinWarmStartBasis* ws = NULL;
  int numcols = getNumCols();
  int numrows = getNumRows();
  int *cstat = new int[numcols];
  int *rstat = new int[numrows];
  int i;

  assert(!probtypemip_);
  
  GUROBI_CALL( "getWarmStart", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getWarmStart", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, cstat) );
  GUROBI_CALL( "getWarmStart", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, rstat) );
  
	ws = new CoinWarmStartBasis;
	ws->setSize( numcols, numrows );

	char sense;
	for( i = 0; i < numrows; ++i )
	{
	  switch( rstat[i] )
	  {
	  	case GRB_BASIC:
	  		ws->setArtifStatus( i, CoinWarmStartBasis::basic );
	  		break;
	  	case GRB_NONBASIC_LOWER:
	  	case GRB_NONBASIC_UPPER:
	  	  GUROBI_CALL( "getWarmStart", GRBgetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, i, &sense) );
	  		ws->setArtifStatus( i, (sense == '>' ? CoinWarmStartBasis::atUpperBound : CoinWarmStartBasis::atLowerBound) );
	  		break;
	  	default:  // unknown row status
	  		delete ws;
	  		delete[] rstat;
	  		delete[] cstat;
	  		return NULL;
	  }
	}
	
	for( i = 0; i < numcols; ++i )
	{
		switch( cstat[i] )
		{
			case GRB_BASIC:
				ws->setStructStatus( i, CoinWarmStartBasis::basic );
				break;
			case GRB_NONBASIC_LOWER:
				ws->setStructStatus( i, CoinWarmStartBasis::atLowerBound );
				break;
			case GRB_NONBASIC_UPPER:
				ws->setStructStatus( i, CoinWarmStartBasis::atUpperBound );
				break;
			case GRB_SUPERBASIC:
				ws->setStructStatus( i, CoinWarmStartBasis::isFree );
				break;
			default:  // unknown column status
				delete[] rstat;
				delete[] cstat;
				delete ws;
				return NULL;
		}
	}

	delete[] cstat;
  delete[] rstat;

  return ws;
}

//-----------------------------------------------------------------------------

bool OsiGrbSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  debugMessage("OsiGrbSolverInterface::setWarmStart(%p)\n", (void*)warmstart);

  const CoinWarmStartBasis* ws = dynamic_cast<const CoinWarmStartBasis*>(warmstart);
  int numcols, numrows, i;
  int *stat;
  bool retval = true;

  if( !ws )
    return false;

  numcols = ws->getNumStructural();
  numrows = ws->getNumArtificial();
  
  if( numcols != getNumCols() || numrows != getNumRows() )
    return false;

  switchToLP();

  stat = new int[numcols];
  for( i = 0; i < numrows; ++i )
  {
  	switch( ws->getArtifStatus( i ) )
  	{
  		case CoinWarmStartBasis::basic:
  			stat[i] = GRB_BASIC;
  			break;
  		case CoinWarmStartBasis::atLowerBound:
  		case CoinWarmStartBasis::atUpperBound:
  			stat[i] = GRB_NONBASIC_LOWER;
  			break;
  		default:  // unknown row status
  			delete[] stat;
  			return false;
  	}
  }
  
  GUROBI_CALL( "setWarmStart", GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_CBASIS, 0, numrows, stat) );

  if (numcols > numrows)
  {
  	delete[] stat; 
    stat = new int[numcols];
  }  
  
  for( i = 0; i < numcols; ++i )
  {
  	switch( ws->getStructStatus( i ) )
  	{
  		case CoinWarmStartBasis::basic:
  			stat[i] = GRB_BASIC;
  			break;
  		case CoinWarmStartBasis::atLowerBound:
  			stat[i] = GRB_NONBASIC_LOWER;
  			break;
  		case CoinWarmStartBasis::atUpperBound:
  			stat[i] = GRB_NONBASIC_UPPER;
  			break;
  		case CoinWarmStartBasis::isFree:
  		default:  // unknown row status
  			delete[] stat;
  			return false;
  	}
  }
  
  GUROBI_CALL( "setWarmStart", GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_VBASIS, 0, numcols, stat) );

  delete[] stat;
  return retval;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiGrbSolverInterface::markHotStart()
{
  debugMessage("OsiGrbSolverInterface::markHotStart()\n");

  int numcols, numrows;

  assert(!probtypemip_);

  numcols = getNumCols();
  numrows = getNumRows();
  if( numcols > hotStartCStatSize_ )
  {
  	delete[] hotStartCStat_;
  	hotStartCStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numcols ) ); // get some extra space for future hot starts
  	hotStartCStat_ = new int[hotStartCStatSize_];
  }
  if( numrows > hotStartRStatSize_ )
  {
  	delete[] hotStartRStat_;
  	hotStartRStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numrows ) ); // get some extra space for future hot starts
  	hotStartRStat_ = new int[hotStartRStatSize_];
  }

  GUROBI_CALL( "markHotStart", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "markHotStart", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, hotStartCStat_) );
  GUROBI_CALL( "markHotStart", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, hotStartRStat_) );
}

void OsiGrbSolverInterface::solveFromHotStart()
{
  debugMessage("OsiGrbSolverInterface::solveFromHotStart()\n");

  double maxiter;

  switchToLP();

  assert( getNumCols() <= hotStartCStatSize_ );
  assert( getNumRows() <= hotStartRStatSize_ );

  GUROBI_CALL( "solveFromHotStart", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "solveFromHotStart", GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_CBASIS, 0, getNumRows(), hotStartRStat_ ) );
  GUROBI_CALL( "solveFromHotStart", GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_VBASIS, 0, getNumCols(), hotStartCStat_ ) );

  GUROBI_CALL( "solveFromHotStart", GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, &maxiter) );
  GUROBI_CALL( "solveFromHotStart", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, (double)hotStartMaxIteration_) );
  
  resolve();

  GUROBI_CALL( "solveFromHotStart", GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, maxiter) );
}

void OsiGrbSolverInterface::unmarkHotStart()
{
  debugMessage("OsiGrbSolverInterface::unmarkHotStart()\n");

  // ??? be lazy with deallocating memory and do nothing here, deallocate memory in the destructor
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows, columns, elements, ...
//------------------------------------------------------------------
int OsiGrbSolverInterface::getNumCols() const
{
  debugMessage("OsiGrbSolverInterface::getNumCols()\n");
  
  int numcols;

  GUROBI_CALL( "getNumCols", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getNumCols", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMVARS, &numcols) );

  return numcols;
}

int OsiGrbSolverInterface::getNumRows() const
{
  debugMessage("OsiGrbSolverInterface::getNumRows()\n");

  int numrows;

  GUROBI_CALL( "getNumRows", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getNumRows", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMCONSTRS, &numrows) );

  return numrows;
}

int OsiGrbSolverInterface::getNumElements() const
{
  debugMessage("OsiGrbSolverInterface::getNumElements()\n");

  int numnz;

  GUROBI_CALL( "getNumElements", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getNumElements", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMNZS, &numnz) );

  return numnz;
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  

const double * OsiGrbSolverInterface::getColLower() const
{
  debugMessage("OsiGrbSolverInterface::getColLower()\n");

  if( collower_ == NULL )
  {
  	int ncols = getNumCols();
  	if( ncols > 0 )
  	{
  		collower_ = new double[ncols];
  		GUROBI_CALL( "getColLower", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getColLower", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_LB, 0, ncols, collower_) );
  	}
  }
  
  return collower_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getColUpper() const
{
  debugMessage("OsiGrbSolverInterface::getColUpper()\n");

  if( colupper_ == NULL )
  {
  	int ncols = getNumCols();
  	if( ncols > 0 )
  	{
  		colupper_ = new double[ncols];
  		GUROBI_CALL( "getColUpper", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getColUpper", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_UB, 0, ncols, colupper_) );
  	}
  }
  
  return colupper_;
}

//------------------------------------------------------------------
const char * OsiGrbSolverInterface::getRowSense() const
{
  debugMessage("OsiGrbSolverInterface::getRowSense()\n");

  if ( rowsense_==NULL )
  {
  	int nrows = getNumRows();
  	if( nrows > 0 )
  	{
  		rowsense_ = new char[nrows];
  		GUROBI_CALL( "getRowSense", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getRowSense", GRBgetcharattrarray(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, 0, nrows, rowsense_) );
  	  
  	  for (int i = 0; i < nrows; ++i)
  	  {
  	  	switch (rowsense_[i])
  	  	{
  	  		case GRB_LESS_EQUAL:
  	  			rowsense_[i] = 'L';
  	  			break;
  	  		case GRB_GREATER_EQUAL:
  	  			rowsense_[i] = 'G';
  	  			break;
  	  		case GRB_EQUAL:
  	  			rowsense_[i] = 'E';
  	  			break;
  	  	}
  	  }
  	}
  }
  
  return rowsense_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRightHandSide() const
{
  debugMessage("OsiGrbSolverInterface::getRightHandSide()\n");

  if ( rhs_==NULL )
  {
  	int nrows = getNumRows();
  	if( nrows > 0 )
  	{
  		rhs_ = new double[nrows];
  		GUROBI_CALL( "getRightHandSide", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getRightHandSide", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_RHS, 0, nrows, rhs_) );
  	}
  }
  
  return rhs_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRowRange() const
{
  debugMessage("OsiGrbSolverInterface::getRowRange()\n");

  if ( rowrange_==NULL ) 
  {
  	// Gurobi does not have row ranges, or in other words, they are always 0.0
  	int nrows = getNumRows();
  	if (nrows > 0)
    	rowrange_ = CoinCopyOfArrayOrZero((double*)NULL, nrows);
  }
  
  return rowrange_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRowLower() const
{
  debugMessage("OsiGrbSolverInterface::getRowLower()\n");

  if ( rowlower_ == NULL )
  {
  	int nrows = getNumRows();
  	if ( nrows > 0 )
  	{
    	const char*   rowsense = getRowSense();
    	const double* rhs      = getRightHandSide();
    	
  		rowlower_ = new double[nrows];

  		double dum1;
  		for ( int i = 0;  i < nrows;  i++ )
  			convertSenseToBound( rowsense[i], rhs[i], 0.0, rowlower_[i], dum1 );
  	}
  }
  
  return rowlower_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRowUpper() const
{  
  debugMessage("OsiGrbSolverInterface::getRowUpper()\n");

  if ( rowupper_ == NULL )
  {
  	int nrows = getNumRows();
  	if ( nrows > 0 ) 
  	{
    	const   char    *rowsense = getRowSense();
    	const   double  *rhs      = getRightHandSide();
  		rowupper_ = new double[nrows];

  		double dum1;
  		for ( int i = 0;  i < nrows;  i++ )
  			convertSenseToBound( rowsense[i], rhs[i], 0.0, dum1, rowupper_[i] );
  	}
  }
  
  return rowupper_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getObjCoefficients() const
{
  debugMessage("OsiGrbSolverInterface::getObjCoefficients()\n");

  if ( obj_==NULL )
  {
  	int ncols = getNumCols();
  	if( ncols > 0 )
  	{
  		obj_ = new double[ncols];
  		GUROBI_CALL( "getObjCoefficients", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getObjCoefficients", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_OBJ, 0, ncols, obj_) );
  	}
  }
  
  return obj_;
}

//------------------------------------------------------------------
double OsiGrbSolverInterface::getObjSense() const
{
  debugMessage("OsiGrbSolverInterface::getObjSense()\n");

  int sense;
  GUROBI_CALL( "getObjSense", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getObjSense", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &sense) );
 
  return (double)sense;
}

//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------

bool OsiGrbSolverInterface::isContinuous( int colNumber ) const
{
  debugMessage("OsiGrbSolverInterface::isContinuous(%d)\n", colNumber);

  return getCtype()[colNumber] == 'C';
}

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------

const CoinPackedMatrix * OsiGrbSolverInterface::getMatrixByRow() const
{
  debugMessage("OsiGrbSolverInterface::getMatrixByRow()\n");
  
  if ( matrixByRow_ == NULL ) 
  {
  	int nrows = getNumRows();
  	int ncols = getNumCols();
  	
    if ( nrows == 0 ) {
    	matrixByRow_ = new CoinPackedMatrix();
    	matrixByRow_->setDimensions(0, ncols);
    	return matrixByRow_;
    }

  	int nelems;
  	int *starts   = new int   [nrows + 1];
  	int *len      = new int   [nrows];

  	GUROBI_CALL( "getMatrixByRow", GRBupdatemodel(getMutableLpPtr()) );

  	GUROBI_CALL( "getMatrixByRow", GRBgetconstrs(getMutableLpPtr(), &nelems, NULL, NULL, NULL, 0, nrows) );

  	assert( nelems == getNumElements() );
  	int     *indices  = new int   [nelems];
  	double  *elements = new double[nelems]; 

  	GUROBI_CALL( "getMatrixByRow", GRBgetconstrs(getMutableLpPtr(), &nelems, starts, indices, elements, 0, nrows) );

  	matrixByRow_ = new CoinPackedMatrix();

  	// Should be able to pass null for length of packed matrix,
  	// assignMatrix does not seem to allow (even though documentation
  	// say it is possible to do this). 
  	// For now compute the length.
  	starts[nrows] = nelems;
  	for ( int i = 0; i < nrows; ++i )
  		len[i] = starts[i+1] - starts[i];

  	matrixByRow_->assignMatrix( false /* not column ordered */,
  			ncols, nrows, nelems,
  			elements, indices, starts, len);
  }
  
  return matrixByRow_;
}

//------------------------------------------------------------------

const CoinPackedMatrix * OsiGrbSolverInterface::getMatrixByCol() const
{
	debugMessage("OsiGrbSolverInterface::getMatrixByCol()\n");

	if ( matrixByCol_ == NULL )
	{
		int nrows = getNumRows();
		int ncols = getNumCols();

    matrixByCol_ = new CoinPackedMatrix();

		if ( ncols > 0 )
		{
		  int nelems;
		  int *starts = new int   [ncols + 1];
		  int *len    = new int   [ncols];

		  GUROBI_CALL( "getMatrixByCol", GRBupdatemodel(getMutableLpPtr()) );

		  GUROBI_CALL( "getMatrixByCol", GRBgetvars(getMutableLpPtr(), &nelems, NULL, NULL, NULL, 0, ncols) );

		  int     *indices  = new int   [nelems];
		  double  *elements = new double[nelems];

		  GUROBI_CALL( "getMatrixByCol", GRBgetvars(getMutableLpPtr(), &nelems, starts, indices, elements, 0, ncols) );

		  // Should be able to pass null for length of packed matrix,
		  // assignMatrix does not seem to allow (even though documentation
		  // say it is possible to do this).
		  // For now compute the length.
		  starts[ncols] = nelems;
		  for ( int i = 0; i < ncols; i++ )
		    len[i] = starts[i+1] - starts[i];

		  matrixByCol_->assignMatrix( true /* column ordered */,
		      nrows, ncols, nelems,
		      elements, indices, starts, len);

		  assert( matrixByCol_->getNumCols()==ncols );
		  assert( matrixByCol_->getNumRows()==nrows );
		}
		else
		  matrixByCol_->setDimensions(nrows, ncols);
	}
	
	return matrixByCol_;
} 

//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double OsiGrbSolverInterface::getInfinity() const
{
  debugMessage("OsiGrbSolverInterface::getInfinity()\n");

  return GRB_INFINITY;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

// *FIXME*: what should be done if a certain vector doesn't exist???

const double * OsiGrbSolverInterface::getColSolution() const
{
	debugMessage("OsiGrbSolverInterface::getColSolution()\n");

	if( colsol_==NULL )
	{
		int ncols = getNumCols();
		if( ncols > 0 )
		{
			colsol_ = new double[ncols];

			GUROBI_CALL( "getColSolution", GRBupdatemodel(getMutableLpPtr()) );

			GUROBI_CALL( "getColSolution", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_X, 0, ncols, colsol_) );
		}
	}
	
	return colsol_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRowPrice() const
{
  debugMessage("OsiGrbSolverInterface::getRowPrice()\n");

  if( rowsol_==NULL )
  {
  	int nrows = getNumRows();
  	if( nrows > 0 )
  	{
  		rowsol_ = new double[nrows];
  		
  		GUROBI_CALL( "getRowPrice", GRBupdatemodel(getMutableLpPtr()) );
  	  
  		GUROBI_CALL( "getRowPrice", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_PI, 0, nrows, rowsol_) );
  	}
  }
  
  return rowsol_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getReducedCost() const
{
  debugMessage("OsiGrbSolverInterface::getReducedCost()\n");

  if( redcost_==NULL )
  {
  	int ncols = getNumCols();
  	if( ncols > 0 )
  	{
  		redcost_ = new double[ncols];
  		
  		GUROBI_CALL( "getReducedCost", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getReducedCost", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_RC, 0, ncols, redcost_) );
  	}
  }
  
  return redcost_;
}

//------------------------------------------------------------------
const double * OsiGrbSolverInterface::getRowActivity() const
{
  debugMessage("OsiGrbSolverInterface::getRowActivity()\n");

  if( rowact_==NULL )
  {
  	int nrows = getNumRows();
  	if( nrows > 0 )
  	{
  		rowact_ = new double[nrows];
  		
  		GUROBI_CALL( "getRowActivity", GRBupdatemodel(getMutableLpPtr()) );

  		GUROBI_CALL( "getRowActivity", GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_SLACK, 0, nrows, rowact_) );
	  	
  		for( int r = 0; r < nrows; ++r )
  		  rowact_[r] += getRightHandSide()[r];
  	}
  }
  return rowact_;
}

//------------------------------------------------------------------
double OsiGrbSolverInterface::getObjValue() const
{
  debugMessage("OsiGrbSolverInterface::getObjValue()\n");

  double objval = 0.0;
 
  GUROBI_CALL( "getObjValue", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getObjValue", GRBgetdblattr(getMutableLpPtr(), GRB_DBL_ATTR_OBJVAL, &objval) );

  // Adjust objective function value by constant term in objective function
  double objOffset;
  getDblParam(OsiObjOffset,objOffset);
  objval = objval - objOffset;

  return objval;
}

//------------------------------------------------------------------
int OsiGrbSolverInterface::getIterationCount() const
{
  debugMessage("OsiGrbSolverInterface::getIterationCount()\n");

  double itercnt;
  
  GUROBI_CALL( "getIterationCount", GRBupdatemodel(getMutableLpPtr()) );

  GUROBI_CALL( "getIterationCount", GRBgetdblattr(getMutableLpPtr(), GRB_DBL_ATTR_ITERCOUNT, &itercnt) );

  return (int)itercnt;
}

//------------------------------------------------------------------
std::vector<double*> OsiGrbSolverInterface::getDualRays(int maxNumRays) const
{
  debugMessage("OsiGrbSolverInterface::getDualRays(%d)\n", maxNumRays);

  OsiGrbSolverInterface solver(*this);

  const int numcols = getNumCols();
  const int numrows = getNumRows();
  int* index = new int[CoinMax(numcols,numrows)];
  int i;
  for ( i = CoinMax(numcols,numrows)-1; i >= 0; --i) {
  	index[i] = i;
  }
  double* obj = new double[CoinMax(numcols,2*numrows)];
  CoinFillN(obj, numcols, 0.0);
  solver.setObjCoeffSet(index, index+numcols, obj);

  double* clb = new double[2*numrows];
  double* cub = new double[2*numrows];

  const double plusone = 1.0;
  const double minusone = -1.0;
  const char* sense = getRowSense();

  const CoinPackedVectorBase** cols = new const CoinPackedVectorBase*[numrows];
  int newcols = 0;
  for (i = 0; i < numrows; ++i) {
  	switch (sense[i]) {
  		case 'L':
  			cols[newcols++] = new CoinShallowPackedVector(1, &index[i], &minusone, false);
  			break;
  		case 'G':
  			cols[newcols++] = new CoinShallowPackedVector(1, &index[i], &plusone, false);
  			break;
  		case 'N':
  		case 'E':
  			break;
  	}
  }

  CoinFillN(obj, newcols, 1.0);
  CoinFillN(clb, newcols, 0.0);
  CoinFillN(cub, newcols, getInfinity());

  solver.addCols(newcols, cols, clb, cub, obj);
  delete[] index;
  delete[] cols;
  delete[] clb;
  delete[] cub;
  delete[] obj;

  solver.setObjSense(1.0); // minimize
  solver.initialSolve();

  const double* solverpi = solver.getRowPrice();
  double* pi = new double[numrows];
  for ( i = numrows - 1; i >= 0; --i) {
  	pi[i] = -solverpi[i];
  }
  return std::vector<double*>(1, pi);
}

//------------------------------------------------------------------
std::vector<double*> OsiGrbSolverInterface::getPrimalRays(int maxNumRays) const
{
  debugMessage("OsiGrbSolverInterface::getPrimalRays(%d)\n", maxNumRays);

  return std::vector<double*>();
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

void OsiGrbSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
  debugMessage("OsiGrbSolverInterface::setObjCoeff(%d, %g)\n", elementIndex, elementValue);
  
  GUROBI_CALL( "setObjCoeff", GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_OBJ, elementIndex, elementValue) );

  if(obj_ != NULL)
    obj_[elementIndex] = elementValue;
}

//-----------------------------------------------------------------------------

void OsiGrbSolverInterface::setObjCoeffSet(const int* indexFirst,
					   const int* indexLast,
					   const double* coeffList)
{
  debugMessage("OsiGrbSolverInterface::setObjCoeffSet(%p, %p, %p)\n", (void*)indexFirst, (void*)indexLast, (void*)coeffList);

  const int cnt = (int)(indexLast - indexFirst);
  
  GUROBI_CALL( "setObjCoeffSet", GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_OBJ, cnt, const_cast<int*>(indexFirst), const_cast<double*>(coeffList)) );

	if (obj_ != NULL)
		for (int i = 0; i < cnt; ++i)
			obj_[indexFirst[i]] = coeffList[i];
	
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColLower(int elementIndex, double elementValue)
{
  debugMessage("OsiGrbSolverInterface::setColLower(%d, %g)\n", elementIndex, elementValue);

  GUROBI_CALL( "setColLower", GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_LB, elementIndex, elementValue) );

  if(collower_ != NULL)
    collower_[elementIndex] = elementValue;
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColUpper(int elementIndex, double elementValue)
{  
  debugMessage("OsiGrbSolverInterface::setColUpper(%d, %g)\n", elementIndex, elementValue);

  GUROBI_CALL( "setColUpper", GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_UB, elementIndex, elementValue) );
	
  if(colupper_ != NULL)
    colupper_[elementIndex] = elementValue;

}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColBounds( int elementIndex, double lower, double upper )
{
  debugMessage("OsiGrbSolverInterface::setColBounds(%d, %g, %g)\n", elementIndex, lower, upper);

  setColLower(elementIndex, lower);
  setColUpper(elementIndex, upper);
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  debugMessage("OsiGrbSolverInterface::setColSetBounds(%p, %p, %p)\n", (void*)indexFirst, (void*)indexLast, (void*)boundList);

  const int cnt = (int)(indexLast - indexFirst);
  if (cnt <= 0)
  	return;
  
  double* lbList = new double[cnt];
  double* ubList = new double[cnt];
  
  for (int i = 0; i < cnt; ++i)
  {
  	lbList[i] = boundList[2*i];
  	ubList[i] = boundList[2*i+1];

    if(collower_ != NULL)
    	collower_[indexFirst[i]] = boundList[2*i];
    if(colupper_ != NULL)
    	colupper_[indexFirst[i]] = boundList[2*i+1];
  }
  
  GUROBI_CALL( "setColSetBounds", GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_LB, cnt, const_cast<int*>(indexFirst), lbList) );

  GUROBI_CALL( "setColSetBounds", GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_UB, cnt, const_cast<int*>(indexFirst), ubList) );

	delete[] lbList;
	delete[] ubList;
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowLower( int i, double elementValue )
{
  debugMessage("OsiGrbSolverInterface::setRowLower(%d, %g)\n", i, elementValue);

  double rhs   = getRightHandSide()[i];
  double range = 0.0;
  char   sense = getRowSense()[i];
  double lower = 0, upper = 0;

  convertSenseToBound( sense, rhs, 0.0, lower, upper );
  if( lower != elementValue ) {
  	convertBoundToSense( elementValue, upper, sense, rhs, range );
  	if (sense == 'R')
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowLower", "OsiGrbSolverInterface", __FILE__, __LINE__);
  	setRowType( i, sense, rhs, range );
  }
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowUpper( int i, double elementValue )
{
  debugMessage("OsiGrbSolverInterface::setRowUpper(%d, %g)\n", i, elementValue);

  double rhs   = getRightHandSide()[i];
  double range = 0.0;
  char   sense = getRowSense()[i];
  double lower = 0, upper = 0;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( upper != elementValue ) {
  	convertBoundToSense( lower, elementValue, sense, rhs, range );
  	if (sense == 'R')
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowUpper", "OsiGrbSolverInterface", __FILE__, __LINE__);
  	setRowType( i, sense, rhs, range );
  }
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
  debugMessage("OsiGrbSolverInterface::setRowBounds(%d, %g, %g)\n", elementIndex, lower, upper);

  double rhs, range;
  char sense;
  
  convertBoundToSense( lower, upper, sense, rhs, range );
  if (sense == 'R')
		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface", __FILE__, __LINE__);
  
  setRowType( elementIndex, sense, rhs, range );
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowType(int i, char sense, double rightHandSide, double range)
{
  debugMessage("OsiGrbSolverInterface::setRowType(%d, %c, %g, %g)\n", i, sense, rightHandSide, range);
  
  if (range)
		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface", __FILE__, __LINE__);

  char grbsense;
  switch (sense)
  {
  	case 'R':
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface", __FILE__, __LINE__);
  		
  	case 'N':
      grbsense = GRB_LESS_EQUAL;
      rightHandSide = getInfinity();
      break;
      
  	case 'L':
  		grbsense = GRB_LESS_EQUAL;
  		break;
  		
  	case 'G':
  		grbsense = GRB_GREATER_EQUAL;
  		break;
  		
  	case 'E':
  		grbsense = GRB_EQUAL;
  		break;
  		
  	default:
  	  std::cerr << "Unknown row sense: " << sense << std::endl;
  	  exit(-1);
  }
  
  GUROBI_CALL( "setRowType", GRBsetcharattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_CHAR_ATTR_SENSE, i, grbsense ) );

  if(rowsense_ != NULL)
    rowsense_[i] = sense;
  
  GUROBI_CALL( "setRowType", GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_RHS, i, rightHandSide) );

  if(rhs_ != NULL)
    rhs_[i] = rightHandSide;
  
  if (rowlower_ != NULL || rowupper_ != NULL)
  {
  	double dummy;
  	convertSenseToBound(sense, rightHandSide, 0.0, 
  			rowlower_ ? rowlower_[i] : dummy,
  			rowupper_ ? rowupper_[i] : dummy);
  }
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  debugMessage("OsiGrbSolverInterface::setRowSetBounds(%p, %p, %p)\n", (void*)indexFirst, (void*)indexLast, (void*)boundList);

  const long int cnt = indexLast - indexFirst;
  if (cnt <= 0)
    return;

  char* sense = new char[cnt];
  double* rhs = new double[cnt];
  double* range = new double[cnt];
  for (int i = 0; i < cnt; ++i) {
    convertBoundToSense(boundList[2*i], boundList[2*i+1], sense[i], rhs[i], range[i]);
    if (range[i])
      throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface", __FILE__, __LINE__);
  }
  setRowSetTypes(indexFirst, indexLast, sense, rhs, range);

  delete[] range;
  delete[] rhs;
  delete[] sense;
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  debugMessage("OsiGrbSolverInterface::setRowSetTypes(%p, %p, %p, %p, %p)\n", 
  		(void*)indexFirst, (void*)indexLast, (void*)senseList, (void*)rhsList, (void*)rangeList);

   const int cnt = (int)(indexLast - indexFirst);
   if (cnt <= 0)
  	 return;

   char* grbsense = new char[cnt];
   double* rhs = new double[cnt];
   for (int i = 0; i < cnt; ++i) {
  	 rhs[i] = rhsList[i];
  	 switch (senseList[i])
  	 {
  		 case 'R':
    		 throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface", __FILE__, __LINE__);
    		 
  		 case 'N':
  			 grbsense[i] = GRB_LESS_EQUAL;
  			 rhs[i] = getInfinity();
  			 break;
  			 
  		 case 'L':
  			 grbsense[i] = GRB_LESS_EQUAL;
  			 break;
  			 
  		 case 'G':
  			 grbsense[i] = GRB_GREATER_EQUAL;
  			 break;
  			 
  		 case 'E':
  			 grbsense[i] = GRB_EQUAL;
  			 break;
  	 }
  	 
     if(rowsense_ != NULL)
  		 rowsense_[indexFirst[i]] = senseList[i];

     if(rhs_ != NULL)
    	 rhs_[indexFirst[i]] = rhs[i];
   }
   
   GUROBI_CALL( "setRowSetTypes", GRBsetcharattrlist(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW), GRB_CHAR_ATTR_SENSE, cnt, const_cast<int*>(indexFirst), grbsense) );

   GUROBI_CALL( "setRowSetTypes", GRBsetdblattrlist(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW), GRB_DBL_ATTR_RHS, cnt, const_cast<int*>(indexFirst), rhs) );

   delete[] rhs;
   delete[] grbsense;
   delete[] grbsense;
}

//#############################################################################
void
OsiGrbSolverInterface::setContinuous(int index)
{
  debugMessage("OsiGrbSolverInterface::setContinuous(%d)\n", index);

  assert(coltype_ != NULL);
  assert(coltypesize_ >= getNumCols());

  coltype_[index] = GRB_CONTINUOUS;

  if ( probtypemip_ )
  {
    GUROBI_CALL( "setContinuous", GRBsetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_VTYPE, index, GRB_CONTINUOUS) );
  }
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setInteger(int index)
{
  debugMessage("OsiGrbSolverInterface::setInteger(%d)\n", index);

  assert(coltype_ != NULL);
  assert(coltypesize_ >= getNumCols());

  if( getColLower()[index] == 0.0 && getColUpper()[index] == 1.0 )
     coltype_[index] = GRB_BINARY;
  else
     coltype_[index] = GRB_INTEGER;

  if ( probtypemip_ )
  {
    GUROBI_CALL( "setInteger", GRBsetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_VTYPE, index, coltype_[index]) );
  }
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setContinuous(const int* indices, int len)
{
  debugMessage("OsiGrbSolverInterface::setContinuous(%p, %d)\n", (void*)indices, len);

  for( int i = 0; i < len; ++i )
     setContinuous(indices[i]);
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setInteger(const int* indices, int len)
{
  debugMessage("OsiGrbSolverInterface::setInteger(%p, %d)\n", (void*)indices, len);

  for( int i = 0; i < len; ++i )
     setInteger(indices[i]);
}

//#############################################################################

void OsiGrbSolverInterface::setObjSense(double s) 
{
  debugMessage("OsiGrbSolverInterface::setObjSense(%g)\n", s);
  
  GUROBI_CALL( "setObjSense", GRBsetintattr(getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS ), GRB_INT_ATTR_MODELSENSE, (int)s) );
}

//-----------------------------------------------------------------------------

void OsiGrbSolverInterface::setColSolution(const double * cs) 
{
  debugMessage("OsiGrbSolverInterface::setColSolution(%p)\n", (void*)cs);

  int nc = getNumCols();

  if( cs == NULL )
  	freeCachedResults();
  else if( nc > 0 )
  {
  	// If colsol isn't allocated, then allocate it
  	if ( colsol_ == NULL )
  		colsol_ = new double[nc];

  	// Copy in new col solution.
  	CoinDisjointCopyN( cs, nc, colsol_ );

  	*messageHandler() << "OsiGrb::setColSolution: Gurobi does not allow setting the column solution. Command is ignored." << CoinMessageEol;
  }
}

//-----------------------------------------------------------------------------

void OsiGrbSolverInterface::setRowPrice(const double * rs) 
{
  debugMessage("OsiGrbSolverInterface::setRowPrice(%p)\n", (void*)rs);

  int nr = getNumRows();

  if( rs == NULL )
    freeCachedResults();
  else if( nr > 0 )
  {
  	// If rowsol isn't allocated, then allocate it
  	if ( rowsol_ == NULL )
  		rowsol_ = new double[nr];

  	// Copy in new row solution.
  	CoinDisjointCopyN( rs, nr, rowsol_ );

  	*messageHandler() << "OsiGrb::setRowPrice: Gurobi does not allow setting the row price. Command is ignored." << CoinMessageEol;
  }
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiGrbSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  debugMessage("OsiGrbSolverInterface::addCol(%p, %g, %g, %g)\n", (void*)&vec, collb, colub, obj);

  int nc = getNumCols();
  assert(coltypesize_ >= nc);

  resizeColType(nc + 1);
  coltype_[nc] = GRB_CONTINUOUS;
  
  GUROBI_CALL( "addCol", GRBaddvar(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_ROW ),
  		vec.getNumElements(),
  		const_cast<int*>(vec.getIndices()),
  		const_cast<double*>(vec.getElements()),
  		obj, collb, colub, coltype_[nc], NULL) );
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  debugMessage("OsiGrbSolverInterface::addCols(%d, %p, %p, %p, %p)\n", numcols, (void*)cols, (void*)collb, (void*)colub, (void*)obj);

  int nc = getNumCols();
  assert(coltypesize_ >= nc);

  resizeColType(nc + numcols);
  CoinFillN(&coltype_[nc], numcols, GRB_CONTINUOUS);

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
  
  GUROBI_CALL( "addCols", GRBaddvars(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numcols, nz,
  		start, index, elem,
  		const_cast<double*>(obj),
  		const_cast<double*>(collb), const_cast<double*>(colub),
  		NULL, NULL) );

  delete[] start;
  delete[] elem;
  delete[] index;
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  debugMessage("OsiGrbSolverInterface::deleteCols(%d, %p)\n", num, (void*)columnIndices);

  GUROBI_CALL( "deleteCols", GRBdelvars(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_ROW ), num, const_cast<int*>(columnIndices)) );
  
  if (coltype_)
  {
  	delete[] coltype_;
  	coltype_ = NULL;
  }
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  debugMessage("OsiGrbSolverInterface::addRow(%p, %g, %g)\n", (void*)&vec, rowlb, rowub);

  char sense;
  double rhs, range;

  convertBoundToSense( rowlb, rowub, sense, rhs, range );
  if (sense == 'R' || range)
  	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface", __FILE__, __LINE__);
  
  addRow( vec, sense, rhs, 0.0 );
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  debugMessage("OsiGrbSolverInterface::addRow(%p, %c, %g, %g)\n", (void*)&vec, rowsen, rowrhs, rowrng);

  if (rowsen == 'R' || rowrng)
  	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface", __FILE__, __LINE__);

  char grbsense;
  double rhs = rowrhs;

  switch( rowsen )
  {
  	case 'N':
  		grbsense = GRB_LESS_EQUAL;
  		rhs   = getInfinity();
  		break;
  		
  	case 'L':
  		grbsense = GRB_LESS_EQUAL;
  		break;

  	case 'G':
  		grbsense = GRB_GREATER_EQUAL;
  		break;
  		
  	case 'E':
  		grbsense = GRB_EQUAL;
  		break;
	
  	default:
  	  std::cerr << "Unknown row sense: " << rowsen << std::endl;
  	  exit(-1);
  }
    
  GUROBI_CALL( "addRow", GRBaddconstr(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ),
  		vec.getNumElements(),
  		const_cast<int*>(vec.getIndices()),
  		const_cast<double*>(vec.getElements()),
  		grbsense, rhs, NULL) );
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  debugMessage("OsiGrbSolverInterface::addRows(%d, %p, %p, %p)\n", numrows, (void*)rows, (void*)rowlb, (void*)rowub);

  int i;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

  int* index = new int[nz];
  double* elem = new double[nz];
  int* start = new int[numrows+1];
  char* grbsense = new char[numrows];
  double* rhs = new double[numrows];
  double range;

  nz = 0;
  start[0] = 0;
  for (i = 0; i < numrows; ++i)
  {
    const CoinPackedVectorBase* row = rows[i];
    const int len = row->getNumElements();
    CoinDisjointCopyN(row->getIndices(), len, index+nz);
    CoinDisjointCopyN(row->getElements(), len, elem+nz);
    nz += len;
    start[i+1] = nz;
    
    convertBoundToSense( rowlb[i], rowub[i], grbsense[i], rhs[i], range );
    if (range || grbsense[i] == 'R')
    {
      delete[] start;
      delete[] elem;
      delete[] index;
      delete[] grbsense;
      delete[] rhs;
    	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface", __FILE__, __LINE__);
    }
    
    switch (grbsense[i])
    {
    	case 'N':
    		grbsense[i] = GRB_LESS_EQUAL;
    		rhs[i] = getInfinity();
    		break;
    	case 'L':
    		grbsense[i] = GRB_LESS_EQUAL;
    		break;
    	case 'G':
    		grbsense[i] = GRB_GREATER_EQUAL;
    		break;
    	case 'E':
    		grbsense[i] = GRB_EQUAL;
    		break;
    }
  }
  
  GUROBI_CALL( "addRows", GRBaddconstrs(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numrows, nz,
  		start, index, elem,
  		grbsense, rhs, NULL) );

  delete[] start;
  delete[] elem;
  delete[] index;
  delete[] grbsense;
  delete[] rhs;
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,
			       const double* rowrng)
{
  debugMessage("OsiGrbSolverInterface::addRows(%d, %p, %p, %p, %p)\n", numrows, (void*)rows, (void*)rowsen, (void*)rowrhs, (void*)rowrng);

  int i;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

  int* index = new int[nz];
  double* elem = new double[nz];
  int* start = new int[numrows+1];
  char* grbsense = new char[numrows];
  double* rhs = new double[numrows];

  nz = 0;
  start[0] = 0;
  for (i = 0; i < numrows; ++i)
  {
    const CoinPackedVectorBase* row = rows[i];
    const int len = row->getNumElements();
    CoinDisjointCopyN(row->getIndices(), len, index+nz);
    CoinDisjointCopyN(row->getElements(), len, elem+nz);
    nz += len;
    start[i+1] = nz;
    
    rhs[i] = rowrhs[i];
    if (rowrng[i] || rowsen[i] == 'R')
    {
      delete[] start;
      delete[] elem;
      delete[] index;
      delete[] grbsense;
      delete[] rhs;
    	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface", __FILE__, __LINE__);
    }
    
    switch (rowsen[i])
    {
    	case 'N':
    		grbsense[i] = GRB_LESS_EQUAL;
    		rhs[i] = getInfinity();
    		break;
    	case 'L':
    		grbsense[i] = GRB_LESS_EQUAL;
    		break;
    	case 'G':
    		grbsense[i] = GRB_GREATER_EQUAL;
    		break;
    	case 'E':
    		grbsense[i] = GRB_EQUAL;
    		break;
    }
  }
  
  GUROBI_CALL( "addRows", GRBaddconstrs(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numrows, nz,
  		start, index, elem,
  		grbsense, rhs, NULL) );

  delete[] start;
  delete[] elem;
  delete[] index;
  delete[] grbsense;
  delete[] rhs;
}

//-----------------------------------------------------------------------------
void 
OsiGrbSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  debugMessage("OsiGrbSolverInterface::deleteRows(%d, %p)\n", num, (void*)rowIndices);
  
  GUROBI_CALL( "deleteRows", GRBdelconstrs(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ), num, const_cast<int*>(rowIndices)) );
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiGrbSolverInterface::loadProblem( const CoinPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const double* rowlb, const double* rowub )
{
  debugMessage("OsiGrbSolverInterface::loadProblem(1)(%p, %p, %p, %p, %p, %p)\n", (void*)&matrix, (void*)collb, (void*)colub, (void*)obj, (void*)rowlb, (void*)rowub);

  const double inf = getInfinity();

  int nrows = matrix.getNumRows();
  char   * rowSense = new char  [nrows];
  double * rowRhs   = new double[nrows];
  double * rowRange = new double[nrows];

  int i;
  for ( i = nrows - 1; i >= 0; --i )
  {
  	const double lower = rowlb ? rowlb[i] : -inf;
  	const double upper = rowub ? rowub[i] : inf;
  	convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], rowRange[i] );
  }

  loadProblem( matrix, collb, colub, obj, rowSense, rowRhs, rowRange ); 
  delete [] rowSense;
  delete [] rowRhs;
  delete [] rowRange;
}
			    
//-----------------------------------------------------------------------------

void
OsiGrbSolverInterface::assignProblem( CoinPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      double*& rowlb, double*& rowub )
{
  debugMessage("OsiGrbSolverInterface::assignProblem()\n");

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
OsiGrbSolverInterface::loadProblem( const CoinPackedMatrix& matrix,
				    const double* collb, const double* colub,
				    const double* obj,
				    const char* rowsen, const double* rowrhs,
				    const double* rowrng )
{
	debugMessage("OsiGrbSolverInterface::loadProblem(2)(%p, %p, %p, %p, %p, %p, %p)\n",
			(void*)&matrix, (void*)collb, (void*)colub, (void*)obj, (void*)rowsen, (void*)rowrhs, (void*)rowrng);

	int nc=matrix.getNumCols();
	int nr=matrix.getNumRows();
	int i;

	char* myrowsen = new char[nr];
	double* myrowrhs = new double[nr];
	
	for ( i=0; i<nr; i++ )
	{
		if (rowrng && rowrng[i])
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);
		
		if (rowrhs)
			myrowrhs[i] = rowrhs[i];
		else
			myrowrhs[i] = 0.0;

		if (rowsen)
			switch (rowsen[i])
			{
				case 'R':
					throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);

				case 'N':
					myrowsen[i] = GRB_LESS_EQUAL;
					myrowrhs[i] = getInfinity();
					break;

				case 'L':
					myrowsen[i] = GRB_LESS_EQUAL;
					break;

				case 'G':
					myrowsen[i] = GRB_GREATER_EQUAL;
					break;

				case 'E':
					myrowsen[i] = GRB_EQUAL;
					break;
			}
		else
			myrowsen[i] = 'G';
	}

	// Set column values to defaults if NULL pointer passed
	double * clb;  
	double * cub;
	double * ob;
	if ( collb != NULL )
		clb = const_cast<double*>(collb);
	else
	{
		clb = new double[nc];
		CoinFillN(clb, nc, 0.0);
	}

	if ( colub!=NULL )
		cub = const_cast<double*>(colub);
	else
	{
		cub = new double[nc];
		CoinFillN(cub, nc, getInfinity());
	}

	if ( obj!=NULL )
		ob = const_cast<double*>(obj);
	else
	{
		ob = new double[nc];
		CoinFillN(ob, nc, 0.0);
	}

	bool freeMatrixRequired = false;
	CoinPackedMatrix * m = NULL;
	if ( !matrix.isColOrdered() ) 
	{
		m = new CoinPackedMatrix();
		m->reverseOrderedCopyOf(matrix);
		freeMatrixRequired = true;
	}
	else 
		m = const_cast<CoinPackedMatrix *>(&matrix);

	assert( nc == m->getNumCols() );
	assert( nr == m->getNumRows() );
	assert( m->isColOrdered() ); 
	
	int modelsense;
	GUROBI_CALL( "loadProblem", GRBupdatemodel(getMutableLpPtr()) );
  
	GUROBI_CALL( "loadProblem", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &modelsense) );

	gutsOfDestructor(); // kill old LP, if any

	std::string pn;
	getStrParam(OsiProbName, pn);
	GUROBI_CALL( "loadProblem", GRBloadmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()),
			nc, nr,
			modelsense,
			0.0,
			ob, 
			const_cast<char *>(myrowsen),
			myrowrhs,
			const_cast<int *>(m->getVectorStarts()),
			const_cast<int *>(m->getVectorLengths()),
			const_cast<int *>(m->getIndices()),
			const_cast<double *>(m->getElements()),
			const_cast<double *>(clb), 
			const_cast<double *>(cub), 
			NULL, NULL, NULL) );

	delete[] myrowsen;
	delete[] myrowrhs;
	if ( collb == NULL )
		delete[] clb;
	if ( colub == NULL ) 
		delete[] cub;
	if ( obj   == NULL )
		delete[] ob;
	if ( freeMatrixRequired ) 
		delete m;

	resizeColType(nc);
	CoinFillN(coltype_, nc, GRB_CONTINUOUS);
}
   
//-----------------------------------------------------------------------------

void
OsiGrbSolverInterface::assignProblem( CoinPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      char*& rowsen, double*& rowrhs,
				      double*& rowrng )
{
  debugMessage("OsiGrbSolverInterface::assignProblem()\n");

  loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );
  delete matrix;   matrix = NULL;
  delete[] collb;  collb = NULL;
  delete[] colub;  colub = NULL;
  delete[] obj;    obj = NULL;
  delete[] rowsen; rowsen = NULL;
  delete[] rowrhs; rowrhs = NULL;
  delete[] rowrng; rowrng = NULL;
}

//-----------------------------------------------------------------------------

void
OsiGrbSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub )
{
  debugMessage("OsiGrbSolverInterface::loadProblem(3)()\n");

  const double inf = getInfinity();
  
  char   * rowSense = new char  [numrows];
  double * rowRhs   = new double[numrows];
  double rowRange;;
  
  for ( int i = numrows - 1; i >= 0; --i ) {
    const double lower = rowlb ? rowlb[i] : -inf;
    const double upper = rowub ? rowub[i] : inf;
    convertBoundToSense( lower, upper, rowSense[i], rowRhs[i], rowRange );
    if (rowSense[i] == 'R' || rowRange)
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);
  }

  loadProblem(numcols, numrows, start, index, value, collb, colub, obj, rowSense, rowRhs, NULL);
  
  delete [] rowSense;
  delete [] rowRhs;
}

//-----------------------------------------------------------------------------

void
OsiGrbSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const char* rowsen, const double* rowrhs,
				   const double* rowrng )
{
  debugMessage("OsiGrbSolverInterface::loadProblem(4)(%d, %d, %p, %p, %p, %p, %p, %p, %p, %p, %p)\n",
  		numcols, numrows, (void*)start, (void*)index, (void*)value, (void*)collb, (void*)colub, (void*)obj, (void*)rowsen, (void*)rowrhs, (void*)rowrng);

	int nc = numcols;
	int nr = numrows;
	int i;

	char* myrowsen = new char[nr];
	double* myrowrhs = new double[nr];
	
	for ( i=0; i<nr; i++ )
	{
		if (rowrng && rowrng[i])
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);
		
		if (rowrhs)
			myrowrhs[i] = rowrhs[i];
		else
			myrowrhs[i] = 0.0;

		if (rowsen)
			switch (rowsen[i])
			{
				case 'R':
					throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);

				case 'N':
					myrowsen[i] = GRB_LESS_EQUAL;
					myrowrhs[i] = getInfinity();
					break;

				case 'L':
					myrowsen[i] = GRB_LESS_EQUAL;
					break;

				case 'G':
					myrowsen[i] = GRB_GREATER_EQUAL;
					break;

				case 'E':
					myrowsen[i] = GRB_EQUAL;
					break;
			}
		else
			myrowsen[i] = 'G';
	}

	// Set column values to defaults if NULL pointer passed
	double * clb;  
	double * cub;
	double * ob;
	if ( collb != NULL )
		clb = const_cast<double*>(collb);
	else
	{
		clb = new double[nc];
		CoinFillN(clb, nc, 0.0);
	}

	if ( colub!=NULL )
		cub = const_cast<double*>(colub);
	else
	{
		cub = new double[nc];
		CoinFillN(cub, nc, getInfinity());
	}

	if ( obj!=NULL )
		ob = const_cast<double*>(obj);
	else
	{
		ob = new double[nc];
		CoinFillN(ob, nc, 0.0);
	}
	
	int* len = new int[nc];
	for (i = 0; i < nc; ++i)
		len[i] = start[i+1] - start[i];

	GUROBI_CALL( "loadProblem", GRBupdatemodel(getMutableLpPtr()) );

	int modelsense;
	GUROBI_CALL( "loadProblem", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &modelsense) );

	gutsOfDestructor(); // kill old LP, if any

	std::string pn;
	getStrParam(OsiProbName, pn);
	GUROBI_CALL( "loadProblem", GRBloadmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()),
			nc, nr,
			modelsense,
			0.0,
			ob, 
			myrowsen,
			myrowrhs,
			const_cast<int*>(start), len,
			const_cast<int*>(index),
			const_cast<double*>(value),
			const_cast<double *>(clb), 
			const_cast<double *>(cub), 
			NULL, NULL, NULL) );

	delete[] myrowsen;
	delete[] myrowrhs;
	if ( collb == NULL )
		delete[] clb;
	if ( colub == NULL ) 
		delete[] cub;
	if ( obj   == NULL )
		delete[] ob;
	delete[] len;

	resizeColType(nc);
	CoinFillN(coltype_, nc, GRB_CONTINUOUS);
}
 
//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------
int OsiGrbSolverInterface::readMps( const char * filename,
				     const char * extension )
{
  debugMessage("OsiGrbSolverInterface::readMps(%s, %s)\n", filename, extension);

  // just call base class method
  return OsiSolverInterface::readMps(filename,extension);
}


//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::writeMps( const char * filename,
				      const char * extension,
				      double objSense ) const
{
  debugMessage("OsiGrbSolverInterface::writeMps(%s, %s, %g)\n", filename, extension, objSense);
 
  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + "." + e;
  //TODO make sure that fullname ends with mps
  GUROBI_CALL( "writeMps", GRBwrite(getMutableLpPtr(), const_cast<char*>(fullname.c_str())) );
}

//#############################################################################
// CPX specific public interfaces
//#############################################################################

GRBenv* OsiGrbSolverInterface::getEnvironmentPtr() const
{
  assert( localenv_ != NULL || globalenv_ != NULL );
  return localenv_ ? localenv_ : globalenv_;
}

GRBmodel* OsiGrbSolverInterface::getLpPtr( int keepCached )
{
  freeCachedData( keepCached );
  return getMutableLpPtr();
}

//-----------------------------------------------------------------------------

const char * OsiGrbSolverInterface::getCtype() const
{
  debugMessage("OsiGrbSolverInterface::getCtype()\n");

  return coltype_;
}

//#############################################################################
// Static instance counter methods
//#############################################################################

void OsiGrbSolverInterface::incrementInstanceCounter()
{
	if ( numInstances_ == 0 && !globalenv_)
	{
	  GUROBI_CALL( "incrementInstanceCounter", GRBloadenv( &globalenv_, NULL ) );
		assert( globalenv_ != NULL );
		globalenv_is_ours = true;
	}
  numInstances_++;
}

//-----------------------------------------------------------------------------

void OsiGrbSolverInterface::decrementInstanceCounter()
{
	assert( numInstances_ != 0 );
	assert( globalenv_ != NULL );
	numInstances_--;
	if ( numInstances_ == 0 && globalenv_is_ours)
	{
		GRBfreeenv( globalenv_ );
		globalenv_ = NULL;
	}
}

//-----------------------------------------------------------------------------

unsigned int OsiGrbSolverInterface::getNumInstances()
{
  return numInstances_;
}

void OsiGrbSolverInterface::setEnvironment(GRBenv* globalenv)
{
	if (numInstances_)
	{
		assert(globalenv_);
		throw CoinError("Cannot set global GUROBI environment, since some OsiGrb instance is still using it.", "setEnvironment", "OsiGrbSolverInterface", __FILE__, __LINE__);
	}
	
	assert(!globalenv_ || !globalenv_is_ours);
	
	globalenv_ = globalenv;
	globalenv_is_ours = false;
}


//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiGrbSolverInterface::OsiGrbSolverInterface(bool use_local_env)
  : OsiSolverInterface(),
    localenv_(NULL),
    lp_(NULL),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(1000000), // ??? default iteration limit for strong branching is large
    obj_(NULL),
    collower_(NULL),
    colupper_(NULL),
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
    matrixByCol_(NULL),
    coltype_(NULL),
    coltypesize_(0),
    probtypemip_(false)
{
  debugMessage("OsiGrbSolverInterface::OsiGrbSolverInterface()\n");
  
  if (use_local_env)
  {
    GUROBI_CALL( "OsiGrbSolverInterface", GRBloadenv( &localenv_, NULL ) );
		assert( localenv_ != NULL );
  }
  else
  	incrementInstanceCounter();
    
  gutsOfConstructor();
  
  // change Osi default to Gurobi default
  setHintParam(OsiDoDualInInitial,true,OsiHintTry);
}

OsiGrbSolverInterface::OsiGrbSolverInterface(GRBenv* localgrbenv)
  : OsiSolverInterface(),
    localenv_(localgrbenv),
    lp_(NULL),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(1000000), // ??? default iteration limit for strong branching is large
    obj_(NULL),
    collower_(NULL),
    colupper_(NULL),
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
    matrixByCol_(NULL),
    coltype_(NULL),
    coltypesize_(0),
    probtypemip_(false)
{
  debugMessage("OsiGrbSolverInterface::OsiGrbSolverInterface()\n");

  if (localenv_ == NULL)
  { // if user called this constructor with NULL pointer, we assume that he meant that a local environment should be created
    GUROBI_CALL( "OsiGrbSolverInterface", GRBloadenv( &localenv_, NULL ) );
		assert( localenv_ != NULL );
  }
    
  gutsOfConstructor();
  
  // change Osi default to Gurobi default
  setHintParam(OsiDoDualInInitial,true,OsiHintTry);
}


//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiSolverInterface * OsiGrbSolverInterface::clone(bool copyData) const
{
  debugMessage("OsiGrbSolverInterface::clone(%d)\n", copyData);

  return( new OsiGrbSolverInterface( *this ) );
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiGrbSolverInterface::OsiGrbSolverInterface( const OsiGrbSolverInterface & source )
  : OsiSolverInterface(source),
    localenv_(NULL),
    lp_(NULL),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(source.hotStartMaxIteration_),
    obj_(NULL),
    collower_(NULL),
    colupper_(NULL),
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
    matrixByCol_(NULL),
    coltype_(NULL),
    coltypesize_(0),
    probtypemip_(false)
{
  debugMessage("OsiGrbSolverInterface::OsiGrbSolverInterface(%p)\n", (void*)&source);

  if (source.localenv_) {
    GUROBI_CALL( "OsiGrbSolverInterface", GRBloadenv( &localenv_, NULL ) );
		assert( localenv_ != NULL );
  } else
  	incrementInstanceCounter();
  
  gutsOfConstructor();
  gutsOfCopy( source );
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiGrbSolverInterface::~OsiGrbSolverInterface()
{
  debugMessage("OsiGrbSolverInterface::~OsiGrbSolverInterface()\n");

  gutsOfDestructor();
  if (localenv_) {
		GRBfreeenv( localenv_ );
  } else
  	decrementInstanceCounter();
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiGrbSolverInterface& OsiGrbSolverInterface::operator=( const OsiGrbSolverInterface& rhs )
{
  debugMessage("OsiGrbSolverInterface::operator=(%p)\n", (void*)&rhs);

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

void OsiGrbSolverInterface::applyColCut( const OsiColCut & cc )
{
	debugMessage("OsiGrbSolverInterface::applyColCut(%p)\n", (void*)&cc);

	const double * gurobiColLB = getColLower();
	const double * gurobiColUB = getColUpper();
	const CoinPackedVector & lbs = cc.lbs();
	const CoinPackedVector & ubs = cc.ubs();
	int i;

	for( i = 0; i < lbs.getNumElements(); ++i )
		if ( lbs.getElements()[i] > gurobiColLB[lbs.getIndices()[i]] )
			setColLower( lbs.getIndices()[i], lbs.getElements()[i] );
	for( i = 0; i < ubs.getNumElements(); ++i )
		if ( ubs.getElements()[i] < gurobiColUB[ubs.getIndices()[i]] )
			setColUpper( ubs.getIndices()[i], ubs.getElements()[i] );
}

//-----------------------------------------------------------------------------

void OsiGrbSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  debugMessage("OsiGrbSolverInterface::applyRowCut(%p)\n", (void*)&rowCut);

  double rhs = 0.0;
  char sns;
  double lb = rowCut.lb();
  double ub = rowCut.ub();
  if( lb <= -getInfinity() && ub >= getInfinity() )   // free constraint
  {
  	rhs = getInfinity();
  	sns = GRB_LESS_EQUAL;
  }
  else if( lb <= -getInfinity() )  // <= constraint
  {
  	rhs = ub;
  	sns = GRB_LESS_EQUAL;
  }
  else if( ub >= getInfinity() )  // >= constraint
  {
  	rhs = lb;
  	sns = GRB_GREATER_EQUAL;
  }
  else if( ub == lb )  // = constraint
  {
  	rhs = ub;
  	sns = GRB_EQUAL;
  }
  else  // range constraint
  {
		throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface", __FILE__, __LINE__);
  }
  
  GUROBI_CALL( "applyRowCut", GRBaddconstr( getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ),
  	rowCut.row().getNumElements(),
		const_cast<int*>( rowCut.row().getIndices() ), 
		const_cast<double*>( rowCut.row().getElements() ),
		sns, rhs, NULL) );
}

//#############################################################################
// Private methods (non-static and static) and static data
//#############################################################################

//------------------------------------------------------------------
// Static data
//------------------------------------------------------------------      
GRBenv* OsiGrbSolverInterface::globalenv_ = NULL;
bool OsiGrbSolverInterface::globalenv_is_ours = true;
unsigned int OsiGrbSolverInterface::numInstances_ = 0;
 
//-------------------------------------------------------------------
// Get pointer to GRBmodel*.
// const methods should use getMutableLpPtr().
// non-const methods should use getLpPtr().
//------------------------------------------------------------------- 
GRBmodel* OsiGrbSolverInterface::getMutableLpPtr() const
{
	if ( lp_ == NULL )
	{
		assert(getEnvironmentPtr() != NULL);
		
		std::string pn;
		getStrParam(OsiProbName, pn);
		GUROBI_CALL( "getMutableLpPtr", GRBnewmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()), 0, NULL, NULL, NULL, NULL, NULL) );
		assert( lp_ != NULL ); 
	}
  return lp_;
}

//-------------------------------------------------------------------

void OsiGrbSolverInterface::gutsOfCopy( const OsiGrbSolverInterface & source )
{
  // Set Rim and constraints
  const double* obj = source.getObjCoefficients();
  const double* rhs = source.getRightHandSide();
  const char* sense = source.getRowSense();
  const CoinPackedMatrix * cols = source.getMatrixByCol();
  const double* lb = source.getColLower();
  const double* ub = source.getColUpper();
  loadProblem(*cols,lb,ub,obj,sense,rhs,source.getRowRange());

  // Set Objective Sense
  setObjSense(source.getObjSense());

  // Set MIP information
  resizeColType(source.coltypesize_);
  CoinDisjointCopyN( source.coltype_, source.coltypesize_, coltype_ );
  
  // Set Solution
  setColSolution(source.getColSolution());
  setRowPrice(source.getRowPrice());
}

//-------------------------------------------------------------------
void OsiGrbSolverInterface::gutsOfConstructor()
{ }

//-------------------------------------------------------------------
void OsiGrbSolverInterface::gutsOfDestructor()
{  
  if ( lp_ != NULL )
  {
    GUROBI_CALL( "gutsOfDestructor", GRBfreemodel(lp_) );
  	lp_ = NULL;
  	freeAllMemory();
  }
  assert( lp_==NULL );
  assert( obj_==NULL );
  assert( collower_==NULL );
  assert( colupper_==NULL );
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
  assert( coltype_==NULL );
  assert( coltypesize_==0 );
}

//-------------------------------------------------------------------
/// free cached vectors

void OsiGrbSolverInterface::freeCachedColRim()
{
  freeCacheDouble( obj_ );  
  freeCacheDouble( collower_ ); 
  freeCacheDouble( colupper_ ); 
  assert( obj_==NULL );
  assert( collower_==NULL );
  assert( colupper_==NULL );
}

void OsiGrbSolverInterface::freeCachedRowRim()
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

void OsiGrbSolverInterface::freeCachedMatrix()
{
  freeCacheMatrix( matrixByRow_ );
  freeCacheMatrix( matrixByCol_ );
  assert( matrixByRow_==NULL ); 
  assert( matrixByCol_==NULL ); 
}

void OsiGrbSolverInterface::freeCachedResults()
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


void OsiGrbSolverInterface::freeCachedData( int keepCached )
{
  if( !(keepCached & OsiGrbSolverInterface::KEEPCACHED_COLUMN) )
    freeCachedColRim();
  if( !(keepCached & OsiGrbSolverInterface::KEEPCACHED_ROW) )
    freeCachedRowRim();
  if( !(keepCached & OsiGrbSolverInterface::KEEPCACHED_MATRIX) )
    freeCachedMatrix();
  if( !(keepCached & OsiGrbSolverInterface::KEEPCACHED_RESULTS) )
    freeCachedResults();
}

void OsiGrbSolverInterface::freeAllMemory()
{
  freeCachedData();
  if( hotStartCStat_ != NULL )
    delete[] hotStartCStat_;
  if( hotStartRStat_ != NULL )
    delete[] hotStartRStat_;
  hotStartCStat_     = NULL;
  hotStartCStatSize_ = 0;
  hotStartRStat_     = NULL;
  hotStartRStatSize_ = 0;
  freeColType();
}

//#############################################################################
// Resets as if default constructor
void 
OsiGrbSolverInterface::reset()
{
  setInitialData(); // clear base class
  gutsOfDestructor();
}

/**********************************************************************/
/* Returns 1 if can just do getBInv etc
   2 if has all OsiSimplex methods
   and 0 if it has none */
int OsiGrbSolverInterface::canDoSimplexInterface() const {
  return 0;
}

/**********************************************************************/
bool OsiGrbSolverInterface::basisIsAvailable() const {
  if (getNumCols() == 0)
  	return true;

  GUROBI_CALL( "basisIsAvailable", GRBupdatemodel(getMutableLpPtr()) );

  int status;
  GUROBI_CALL( "basisIsAvailable", GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &status) );

	if (status == GRB_LOADED || status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD || status == GRB_UNBOUNDED)
		return false;
	
  int dum;
  return GRBgetintattrelement(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, &dum) == 0;
}

/* Osi return codes:
0: free  
1: basic  
2: upper 
3: lower
*/
void OsiGrbSolverInterface::getBasisStatus(int* cstat, int* rstat) const {

	int numcols = getNumCols();
	int numrows = getNumRows();

	GUROBI_CALL( "getBasisStatus", GRBupdatemodel(getMutableLpPtr()) );

	GUROBI_CALL( "getBasisStatus", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, cstat) );
	
	for (int i = 0; i < numcols; ++i)
		switch (cstat[i])
		{
			case GRB_BASIC:
				cstat[i] = 1;
				break;
			case GRB_NONBASIC_LOWER:
				cstat[i] = 3;
				break;
			case GRB_NONBASIC_UPPER:
				cstat[i] = 2;
				break;
			case GRB_SUPERBASIC:
				cstat[i] = 0;
				break;
		}

	GUROBI_CALL( "getBasisStatus", GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, rstat) );

	char sense;
	for (int i = 0; i < numrows; ++i)
		switch (rstat[i])
		{
			case GRB_BASIC:
				rstat[i] = 1;
				break;
			case GRB_NONBASIC_LOWER:
			case GRB_NONBASIC_UPPER:
			  GUROBI_CALL( "getBasisStatus", GRBgetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, i, &sense) );
				rstat[i] = (sense == '>' ? 2 : 3);
				break;
		}
}
