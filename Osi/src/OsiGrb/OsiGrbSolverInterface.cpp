//-----------------------------------------------------------------------------
// name:     OSI Interface for Gurobi
// template: OSI Cplex Interface written by T. Achterberg
// author:   Stefan Vigerske
//           Humboldt University Berlin
// date:     09/02/2009
// comments: please scan this file for '???' and 'TODO' and read the comments
//-----------------------------------------------------------------------------
// Copyright (C) 2009 Humboldt University Berlin and others.
// All Rights Reserved.

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

inline void
checkGRBerror( int err, std::string grbfuncname, std::string osimethod )
{
  if( err != 0 )
    {
      char s[1001];
      sprintf( s, "%s returned error %d ", grbfuncname.c_str(), err );
      if (OsiGrbSolverInterface::globalenv_)
      	strncat(s, GRBgeterrormsg(OsiGrbSolverInterface::globalenv_), 1000);
      std::cout << "ERROR: " << s << " (" << osimethod << " in OsiGrbSolverInterface)" << std::endl;
      throw CoinError( s, osimethod, "OsiGrbSolverInterface" );
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
     
     int err = GRBsetcharattrarray(lp, GRB_CHAR_ATTR_VTYPE, 0, nc, contattr);
     
     delete[] contattr;

     checkGRBerror( err, "GRBsetcharattrarray", "switchToLP" );
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

     int err = GRBsetcharattrarray(lp, GRB_CHAR_ATTR_VTYPE, 0, nc, coltype_);
     checkGRBerror( err, "GRBsetcharattrarray", "switchToMIP" );

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
  int rc;
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

	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_LPMETHOD, algorithm);
	checkGRBerror( rc, "GRBsetintparam", "initialSolve" );

	/* set whether presolve or not */
  int presolve = GRB_PRESOLVE_AUTO;
  gotHint = (getHintParam(OsiDoPresolveInInitial,takeHint,strength));
  assert (gotHint);
  if (strength!=OsiHintIgnore)
  	presolve = takeHint ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF;

	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve);
	checkGRBerror( rc, "GRBsetintparam", "initialSolve" );

	/* set whether output or not */
 	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0));
	checkGRBerror( rc, "GRBsetintparam", "initialSolve" );

	/* optimize */
  rc = GRBoptimize(lp);
  checkGRBerror( rc, "GRBoptimize", "initialSolve" );
  
  /* reoptimize without presolve if status unclear */
  int stat;
  rc = GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "initialSolve" );

  if (stat == GRB_INF_OR_UNBD && presolve != GRB_PRESOLVE_OFF) {
  	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF);
  	checkGRBerror( rc, "GRBsetintparam", "initialSolve" );

    rc = GRBoptimize(lp);
    checkGRBerror( rc, "GRBoptimize", "initialSolve" );
  	
  	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve);
  	checkGRBerror( rc, "GRBsetintparam", "initialSolve" );
  }
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::resolve()
{
  debugMessage("OsiGrbSolverInterface::resolve()\n");
  int rc;
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

  rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_LPMETHOD, algorithm);
	checkGRBerror( rc, "GRBsetintparam", "resolve" );

	/* set whether presolve or not */
  int presolve = GRB_PRESOLVE_OFF;
  gotHint = getHintParam(OsiDoPresolveInResolve,takeHint,strength);
  assert (gotHint);
  if (strength != OsiHintIgnore)
  	presolve = takeHint ? GRB_PRESOLVE_AUTO : GRB_PRESOLVE_OFF;

	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve);
	checkGRBerror( rc, "GRBsetintparam", "resolve" );

	/* set whether output or not */
 	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0));
	checkGRBerror( rc, "GRBsetintparam", "resolve" );

	/* optimize */
  rc = GRBoptimize(lp);
  checkGRBerror( rc, "GRBoptimize", "resolve" );

  /* reoptimize if status unclear */
  int stat;
  rc = GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "resolve" );

  if (stat == GRB_INF_OR_UNBD && presolve != GRB_PRESOLVE_OFF) {
  	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, GRB_PRESOLVE_OFF);
  	checkGRBerror( rc, "GRBsetintparam", "resolve" );

    rc = GRBoptimize(lp);
    checkGRBerror( rc, "GRBoptimize", "resolve" );
  	
  	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_PRESOLVE, presolve);
  	checkGRBerror( rc, "GRBsetintparam", "resolve" );
  }
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::branchAndBound()
{
  debugMessage("OsiGrbSolverInterface::branchAndBound()\n");
  int rc;

  switchToMIP();

  GRBmodel* lp = getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS );

 	rc = GRBsetintparam(GRBgetenv(lp), GRB_INT_PAR_OUTPUTFLAG, (messageHandler()->logLevel() > 0));
	checkGRBerror( rc, "GRBsetintparam", "branchAndBound" );

  rc = GRBoptimize( lp );
  checkGRBerror( rc, "GRBoptimize", "resolve" );
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiGrbSolverInterface::setIntParam(OsiIntParam key, int value)
{
  debugMessage("OsiGrbSolverInterface::setIntParam(%d, %d)\n", key, value);

  int rc;
  bool retval = false;
  switch (key)
  {
    case OsiMaxNumIteration:
    	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, (double)value);
      checkGRBerror( rc, "GRBsetdblparam", "setIntParam" );
      retval = (rc == 0);
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

  bool retval = false;
  int rc;
  switch (key)
  {
//TODO
//  	case OsiDualObjectiveLimit:
//  		break;
  	case OsiPrimalObjectiveLimit:
    	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_CUTOFF, value);
      checkGRBerror( rc, "GRBsetdblparam", "setDblParam" );
  		break;
  	case OsiDualTolerance:
    	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_OPTIMALITYTOL, value);
      checkGRBerror( rc, "GRBsetdblparam", "setDblParam" );
      retval = (rc == 0);
  		break;
  	case OsiPrimalTolerance:
    	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_FEASIBILITYTOL, value);
      checkGRBerror( rc, "GRBsetdblparam", "setDblParam" );
      retval = (rc == 0);
  		break;
  	case OsiObjOffset:
  		retval = OsiSolverInterface::setDblParam(key,value);
  		break;
  	case OsiLastDblParam:
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
OsiGrbSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  debugMessage("OsiGrbSolverInterface::setStrParam(%d, %s)\n", key, value.c_str());

  bool retval=false;
  switch (key) {
  case OsiProbName:
    OsiSolverInterface::setStrParam(key,value);
    return retval = true;
  case OsiSolverName:
    return false;
  case OsiLastStrParam:
    return false;
  default:
    return false ;
  }
  return false;
}

//-----------------------------------------------------------------------------

bool
OsiGrbSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  debugMessage("OsiGrbSolverInterface::getIntParam(%d)\n", key);

  bool retval = false;
  int rc;
  switch (key)
    {
    case OsiMaxNumIteration:
    	double dblval;
      rc = GRBupdatemodel(getMutableLpPtr());
      checkGRBerror(rc, "GRBupdatemodel", "getIntParam");
    	rc = GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, &dblval);
      checkGRBerror( rc, "GRBgetdblparam", "getIntParam" );
    	value = (int) dblval;
      break;
    case OsiMaxNumIterationHotStart:
      value = hotStartMaxIteration_;
      retval = true;
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
OsiGrbSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  debugMessage("OsiGrbSolverInterface::getDblParam(%d)\n", key);

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getDblParam");

  bool retval = false;
  switch (key) 
  {
//TODO
//    	case OsiDualObjectiveLimit:
//    		break;
   	case OsiPrimalObjectiveLimit:
   		rc = GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_CUTOFF, &value);
   		checkGRBerror( rc, "GRBgetdblparam", "getDblParam" );
   		retval = (rc == 0);
  	case OsiDualTolerance:
  		rc = GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_OPTIMALITYTOL, &value);
  		checkGRBerror( rc, "GRBgetdblparam", "getDblParam" );
  		retval = (rc == 0);
  		break;
  	case OsiPrimalTolerance:
  		rc = GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_FEASIBILITYTOL, &value);
  		checkGRBerror( rc, "GRBgetdblparam", "getDblParam" );
  		retval = (rc == 0);
  		break;
    case OsiObjOffset:
      retval = OsiSolverInterface::getDblParam(key, value);
      break;
    case OsiLastDblParam:
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
OsiGrbSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  debugMessage("OsiGrbSolverInterface::getStrParam(%d)\n", key);

  switch (key) {
  	case OsiProbName:
  		OsiSolverInterface::getStrParam(key, value);
  		break;
  	case OsiSolverName:
  		value = "gurobi";
  		break;
  	case OsiLastStrParam:
  		return false;
  	default:
  		return false ;
  }

  return true;
}

//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiGrbSolverInterface::isAbandoned() const
{
  debugMessage("OsiGrbSolverInterface::isAbandoned()\n");

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isAbandoned");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isAbandoned" );
	
	return (
			stat == GRB_LOADED ||
			stat == GRB_NUMERIC || 
			stat == GRB_INTERRUPTED
		);
}

bool OsiGrbSolverInterface::isProvenOptimal() const
{
  debugMessage("OsiGrbSolverInterface::isProvenOptimal()\n");

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isProvenOptimal");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isProvenOptimal" );

	return (stat == GRB_OPTIMAL);
}

bool OsiGrbSolverInterface::isProvenPrimalInfeasible() const
{
  debugMessage("OsiGrbSolverInterface::isProvenPrimalInfeasible()\n");

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isProvenPrimalInfeasible");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isProvenPrimalInfeasible" );

	return (stat == GRB_INFEASIBLE);
}

bool OsiGrbSolverInterface::isProvenDualInfeasible() const
{
  debugMessage("OsiGrbSolverInterface::isProvenDualInfeasible()\n");

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isProvenDualInfeasible");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isProvenDualInfeasible" );

	return (stat == GRB_UNBOUNDED);
}

bool OsiGrbSolverInterface::isPrimalObjectiveLimitReached() const
{
  debugMessage("OsiGrbSolverInterface::isPrimalObjectiveLimitReached()\n");

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isPrimalObjectiveLimitReached");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isPrimalObjectiveLimitReached" );

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

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "isIterationLimitReached");

  int stat;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &stat);
	checkGRBerror( rc, "GRBgetintattr", "isIterationLimitReached" );

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
  int rc, i;

  assert(!probtypemip_);
  
  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getWarmStart");

  rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, cstat);
	checkGRBerror( rc, "GRBgetintattrarray", "getWarmStart" );
	if (rc)
	{
		delete[] cstat;
		delete[] rstat;
		return NULL;
	}

  rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, rstat);
	checkGRBerror( rc, "GRBgetintattrarray", "getWarmStart" );
	if (rc)
	{
		delete[] cstat;
		delete[] rstat;
		return NULL;
	}
  
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
				rc = GRBgetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, i, &sense);
				checkGRBerror( rc, "GRBgetcharattrelement", "getWarmStart" );
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
  int numcols, numrows, i, rc;
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
  
  rc = GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_CBASIS, 0, numrows, stat);
	checkGRBerror( rc, "GRBsetintattrarray", "setWarmStart" );
	if (rc)
		retval = false;

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
  
  rc = GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_VBASIS, 0, numcols, stat);
	checkGRBerror( rc, "GRBsetintattrarray", "setWarmStart" );
	if (rc)
		retval = false;

  delete[] stat;
  return retval;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiGrbSolverInterface::markHotStart()
{
  debugMessage("OsiGrbSolverInterface::markHotStart()\n");

  int rc;
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

  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "markHotStart");

  rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, hotStartCStat_);
	checkGRBerror( rc, "GRBgetintattrarray", "markHotStart" );
  rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, hotStartRStat_);
	checkGRBerror( rc, "GRBgetintattrarray", "markHotStart" );
}

void OsiGrbSolverInterface::solveFromHotStart()
{
  debugMessage("OsiGrbSolverInterface::solveFromHotStart()\n");

  int rc;
  double maxiter;

  switchToLP();

  assert( getNumCols() <= hotStartCStatSize_ );
  assert( getNumRows() <= hotStartRStatSize_ );

  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "solveFromHotStart");

  rc = GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_CBASIS, 0, getNumRows(), hotStartRStat_ );
	checkGRBerror( rc, "GRBsetintattrarray", "solveFromHotStart" );
  rc = GRBsetintattrarray(getLpPtr(OsiGrbSolverInterface::FREECACHED_RESULTS), GRB_INT_ATTR_VBASIS, 0, getNumCols(), hotStartCStat_ );
	checkGRBerror( rc, "GRBsetintattrarray", "solveFromHotStart" );

	rc = GRBgetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, &maxiter);
	checkGRBerror( rc, "GRBgetdblparam", "solveFromHotStart" );

	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, (double)hotStartMaxIteration_);
  checkGRBerror( rc, "GRBsetdblparam", "solveFromHotStart" );
  
  resolve();

	rc = GRBsetdblparam(GRBgetenv(getMutableLpPtr()), GRB_DBL_PAR_ITERATIONLIMIT, maxiter);
  checkGRBerror( rc, "GRBsetdblparam", "solveFromHotStart" );
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
  
  int numcols, rc;

  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getNumCols");

  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMVARS, &numcols);
  checkGRBerror( rc, "GRBgetintattr", "getNumCols" );

  return numcols;
}

int OsiGrbSolverInterface::getNumRows() const
{
  debugMessage("OsiGrbSolverInterface::getNumRows()\n");

  int numrows, rc;

  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getNumRows");

  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMCONSTRS, &numrows);
  checkGRBerror( rc, "GRBgetintattr", "getNumRows" );

  return numrows;
}

int OsiGrbSolverInterface::getNumElements() const
{
  debugMessage("OsiGrbSolverInterface::getNumElements()\n");

  int numnz, rc;

  rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getNumElements");

  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_NUMNZS, &numnz);
  checkGRBerror( rc, "GRBgetintattr", "getNumElements" );

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
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getColLower");
  		rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_LB, 0, ncols, collower_);
  	  checkGRBerror( rc, "GRBgetdblattrarray", "getColLower" );
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
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getColUpper");
  		rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_UB, 0, ncols, colupper_);
  	  checkGRBerror( rc, "GRBgetdblattrarray", "getColUpper" );
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
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getRowSense");
  		rc = GRBgetcharattrarray(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, 0, nrows, rowsense_);
  	  checkGRBerror( rc, "GRBgetcharattrarray", "getRowSense" );
  	  
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
  int rc;

  if ( rhs_==NULL )
  {
  	int nrows = getNumRows();
  	if( nrows > 0 )
  	{
  		rhs_ = new double[nrows];
  	  rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getRightHandSide");
  		rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_RHS, 0, nrows, rhs_);
  	  checkGRBerror( rc, "GRBgetdblattrarray", "getRightHandSide" );
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
    	const   char    *rowsense = getRowSense();
    	const   double  *rhs      = getRightHandSide();
    	
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
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getObjCoefficients");
  		rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_OBJ, 0, ncols, obj_);
  		checkGRBerror( rc, "GRBgetdblattrarray", "getObjCoefficients" );
  	}
  }
  
  return obj_;
}

//------------------------------------------------------------------
double OsiGrbSolverInterface::getObjSense() const
{
  debugMessage("OsiGrbSolverInterface::getObjSense()\n");

  int sense;
  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror( rc, "GRBupdatemodel", "getObjSense");
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &sense);
	checkGRBerror( rc, "GRBgetintattr", "getObjSense" );
 
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

  	int nelems, rc;
  	int *starts   = new int   [nrows + 1];
  	int *len      = new int   [nrows];

    rc = GRBupdatemodel(getMutableLpPtr());
    checkGRBerror(rc, "GRBupdatemodel", "getMatrixByRow");

  	rc = GRBgetconstrs(getMutableLpPtr(), &nelems, NULL, NULL, NULL, 0, nrows);
  	checkGRBerror( rc, "GRBgetconstrs", "getMatrixByRow" );

  	assert( nelems == getNumElements() );
  	int     *indices  = new int   [nelems];
  	double  *elements = new double[nelems]; 

  	rc = GRBgetconstrs(getMutableLpPtr(), &nelems, starts, indices, elements, 0, nrows); 
  	checkGRBerror( rc, "GRBgetconstrs", "getMatrixByRow" );

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
		int nelems, rc;
		int *starts = new int   [ncols + 1];
		int *len    = new int   [ncols];

	  rc = GRBupdatemodel(getMutableLpPtr());
	  checkGRBerror(rc, "GRBupdatemodel", "getMatrixByCol");

  	rc = GRBgetvars(getMutableLpPtr(), &nelems, NULL, NULL, NULL, 0, ncols); 
  	checkGRBerror( rc, "GRBgetvars", "getMatrixByCol" );

		int     *indices  = new int   [nelems];
		double  *elements = new double[nelems]; 

  	rc = GRBgetvars(getMutableLpPtr(), &nelems, starts, indices, elements, 0, ncols); 
  	checkGRBerror( rc, "GRBgetvars", "getMatrixByCol" );

		matrixByCol_ = new CoinPackedMatrix();

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

		  int rc = GRBupdatemodel(getMutableLpPtr());
		  checkGRBerror(rc, "GRBupdatemodel", "getColSolution");

			rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_X, 0, ncols, colsol_);
	  	checkGRBerror( rc, "GRBgetdblattrarray", "getColSolution" );
	  	
	  	if (rc)
	  		CoinFillN( colsol_, ncols, 0.0 );
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
  		
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getRowPrice");
  	  
			rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_PI, 0, nrows, rowsol_);
	  	checkGRBerror( rc, "GRBgetdblattrarray", "getRowPrice" );
  		
	  	if (rc)
	  		CoinFillN( rowsol_, nrows, 0.0 );
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
  		
  	  int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getReducedCost");

			rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_RC, 0, ncols, redcost_);
	  	checkGRBerror( rc, "GRBgetdblattrarray", "getReducedCost" );

	  	if (rc)
  			CoinFillN( redcost_, ncols, 0.0 );
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
  		
  		int rc = GRBupdatemodel(getMutableLpPtr());
  	  checkGRBerror(rc, "GRBupdatemodel", "getRowActivity");

  		rc = GRBgetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_SLACK, 0, nrows, rowact_);
	  	checkGRBerror( rc, "GRBgetdblattrarray", "getRowActivity" );
	  	
	  	if (rc)
				CoinFillN( rowact_, nrows, 0.0 );
	  	else
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
 
  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getObjValue");

	rc = GRBgetdblattr(getMutableLpPtr(), GRB_DBL_ATTR_OBJVAL, &objval);
	checkGRBerror( rc, "GRBgetdblattr", "getObjValue" );

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
  
  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getIterationCount");

	rc = GRBgetdblattr(getMutableLpPtr(), GRB_DBL_ATTR_ITERCOUNT, &itercnt);
	checkGRBerror( rc, "GRBgetdblattr", "getIterationCount" );

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
  
  int rc = GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_OBJ, elementIndex, elementValue);
	checkGRBerror( rc, "GRBsetdblattrelement", "setObjCoeff" );

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
  
  int rc = GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_OBJ, cnt, const_cast<int*>(indexFirst), const_cast<double*>(coeffList));
	checkGRBerror( rc, "GRBsetdblattrlist", "setObjCoeffSet" );

	if (obj_ != NULL)
		for (int i = 0; i < cnt; ++i)
			obj_[indexFirst[i]] = coeffList[i];
	
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColLower(int elementIndex, double elementValue)
{
  debugMessage("OsiGrbSolverInterface::setColLower(%d, %g)\n", elementIndex, elementValue);

  int rc = GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_LB, elementIndex, elementValue);
	checkGRBerror( rc, "GRBsetdblattrelement", "setColLower" );

  if(collower_ != NULL)
    collower_[elementIndex] = elementValue;
}

//-----------------------------------------------------------------------------
void OsiGrbSolverInterface::setColUpper(int elementIndex, double elementValue)
{  
  debugMessage("OsiGrbSolverInterface::setColUpper(%d, %g)\n", elementIndex, elementValue);

  int rc = GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_UB, elementIndex, elementValue);
	checkGRBerror( rc, "GRBsetdblattrelement", "setColUpper" );
	
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
  
  int rc = GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_LB, cnt, const_cast<int*>(indexFirst), lbList);
	checkGRBerror( rc, "GRBsetdblattrlist", "setColSetBounds" );

  rc = GRBsetdblattrlist(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_UB, cnt, const_cast<int*>(indexFirst), ubList);
	checkGRBerror( rc, "GRBsetdblattrlist", "setColSetBounds" );

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
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowLower", "OsiGrbSolverInterface");
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
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowUpper", "OsiGrbSolverInterface");
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
		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface");
  
  setRowType( elementIndex, sense, rhs, range );
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::setRowType(int i, char sense, double rightHandSide, double range)
{
  debugMessage("OsiGrbSolverInterface::setRowType(%d, %c, %g, %g)\n", i, sense, rightHandSide, range);
  
  if (range)
		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface");

  char grbsense;
  switch (sense)
  {
  	case 'R':
  		throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface");
  		
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
  
  int err = GRBsetcharattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_CHAR_ATTR_SENSE, i, grbsense );
  checkGRBerror(err, "GRBsetcharattrelement", "setRowType");

  if(rowsense_ != NULL)
    rowsense_[i] = sense;
  
  err = GRBsetdblattrelement(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_PROBLEM ), GRB_DBL_ATTR_RHS, i, rightHandSide);
  checkGRBerror(err, "GRBsetdblattrelement", "setRowType");

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
  		 throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface");
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
    		 throw CoinError("Ranged rows not supported by Gurobi", "setRowBounds", "OsiGrbSolverInterface");
    		 
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
   
   int err;
   
   err = GRBsetcharattrlist(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW), GRB_CHAR_ATTR_SENSE, cnt, const_cast<int*>(indexFirst), grbsense);
   checkGRBerror(err, "GRBsetcharattrlist", "setRowSetTypes");

   err = GRBsetdblattrlist(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW), GRB_DBL_ATTR_RHS, cnt, const_cast<int*>(indexFirst), rhs);
   checkGRBerror(err, "GRBsetdblattrlist", "setRowSetTypes");

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
  	int err = GRBsetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_VTYPE, index, GRB_CONTINUOUS);
  	checkGRBerror( err, "GRBsetcharattrelement", "setContinuous" );
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
  	int err = GRBsetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_VTYPE, index, coltype_[index]);
  	checkGRBerror( err, "GRBsetcharattrelement", "setInteger" );
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
  
  int rc = GRBsetintattr(getLpPtr( OsiGrbSolverInterface::FREECACHED_RESULTS ), GRB_INT_ATTR_MODELSENSE, (int)s);
  checkGRBerror(rc, "GRBsetintattr", "setObjSense");
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

  	std::cout << "OsiGrb::setColSolution: Gurobi does not allow setting the column solution. Command is ignored." << std::endl;
//  	int err = GRBsetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_X, 0, nc, const_cast<double*>(colsol_));
// 		checkGRBerror( err, "GRBsetdblattrarray", "setColSolution" );
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

   std::cout << "OsiGrb::setRowPrice: Gurobi does not allow setting the row price. Command is ignored." << std::endl;
//  	int err = GRBsetdblattrarray(getMutableLpPtr(), GRB_DBL_ATTR_PI, 0, nr, const_cast<double*>(rowsol_));
// 		checkGRBerror( err, "GRBsetdblattrarray", "setRowPrice" );
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
  
  int err = GRBaddvar(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_ROW ),
  		vec.getNumElements(),
  		const_cast<int*>(vec.getIndices()),
  		const_cast<double*>(vec.getElements()),
  		obj, collb, colub, coltype_[nc], NULL);
  checkGRBerror( err, "GRBaddvar", "addCol" );
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
  
  int err = GRBaddvars(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numcols, nz,
  		start, index, elem,
  		const_cast<double*>(obj),
  		const_cast<double*>(collb), const_cast<double*>(colub),
  		NULL, NULL);
  checkGRBerror( err, "GRBaddvars", "addCols" );

  delete[] start;
  delete[] elem;
  delete[] index;
}

//-----------------------------------------------------------------------------
void
OsiGrbSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  debugMessage("OsiGrbSolverInterface::deleteCols(%d, %p)\n", num, (void*)columnIndices);

  int rc = GRBdelvars(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_ROW ), num, const_cast<int*>(columnIndices));
  checkGRBerror(rc, "GRBdelvars", "deleteCols");
  
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
  	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface");
  
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
  	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface");

  int err;
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
    
  err = GRBaddconstr(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ),
  		vec.getNumElements(),
  		const_cast<int*>(vec.getIndices()),
  		const_cast<double*>(vec.getElements()),
  		grbsense, rhs, NULL);
  checkGRBerror( err, "GRBaddconstr", "addRow" );
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
    	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface");
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
  
  int err = GRBaddconstrs(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numrows, nz,
  		start, index, elem,
  		grbsense, rhs, NULL);
  checkGRBerror( err, "GRBaddconstrs", "addRows" );

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
    	throw CoinError("Ranged rows not supported by Gurobi", "addRow", "OsiGrbSolverInterface");
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
  
  int err = GRBaddconstrs(getLpPtr(OsiGrbSolverInterface::KEEPCACHED_ROW),
  		numrows, nz,
  		start, index, elem,
  		grbsense, rhs, NULL);
  checkGRBerror( err, "GRBaddconstrs", "addRows" );

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
  
  int err = GRBdelconstrs(getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ), num, const_cast<int*>(rowIndices));
  checkGRBerror( err, "GRBdelconstrs", "deleteRows" );
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
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");
		
		if (rowrhs)
			myrowrhs[i] = rowrhs[i];
		else
			myrowrhs[i] = 0.0;

		if (rowsen)
			switch (rowsen[i])
			{
				case 'R':
					throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");

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
  int err = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror( err, "GRBupdatemodel", "loadProblem");
  
	err = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &modelsense);
	checkGRBerror( err, "GRBgetintattr", "loadProblem" );

	gutsOfDestructor(); // kill old LP, if any

	std::string pn;
	getStrParam(OsiProbName, pn);
	err = GRBloadmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()),
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
			NULL, NULL, NULL);
	checkGRBerror( err, "GRBloadModel", "loadProblem" );

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
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");
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
			throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");
		
		if (rowrhs)
			myrowrhs[i] = rowrhs[i];
		else
			myrowrhs[i] = 0.0;

		if (rowsen)
			switch (rowsen[i])
			{
				case 'R':
					throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");

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

  int err = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror( err, "GRBupdatemodel", "loadProblem");

	int modelsense;
	err = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_MODELSENSE, &modelsense);
	checkGRBerror( err, "GRBgetintattr", "loadProblem" );

	gutsOfDestructor(); // kill old LP, if any

	std::string pn;
	getStrParam(OsiProbName, pn);
	err = GRBloadmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()),
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
			NULL, NULL, NULL);
	checkGRBerror( err, "GRBloadModel", "loadProblem" );

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
  int err = GRBwrite(getMutableLpPtr(), const_cast<char*>(fullname.c_str()));
  checkGRBerror( err, "GRBwrite", "writeMps" );
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
		int err = GRBloadenv( &globalenv_, NULL );
		checkGRBerror( err, "GRBloadenv", "incrementInstanceCounter" );
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
		throw CoinError("Cannot set global gurobi environment, since some OsiGrb instance is still using it.", "setEnvironment", "OsiGrbSolverInterface");
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
		int err = GRBloadenv( &localenv_, NULL );
		checkGRBerror( err, "GRBloadenv", "OsiGrbSolverInterface" );
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
		int err = GRBloadenv( &localenv_, NULL );
		checkGRBerror( err, "GRBloadenv", "OsiGrbSolverInterface" );
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
		int err = GRBloadenv( &localenv_, NULL );
		checkGRBerror( err, "GRBloadenv", "OsiGrbSolverInterface" );
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

  int err = 0;
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
		throw CoinError("Ranged rows not supported by Gurobi", "loadProblem", "OsiGrbSolverInterface");
  }
  
  err = GRBaddconstr( getLpPtr( OsiGrbSolverInterface::KEEPCACHED_COLUMN ), 
  	rowCut.row().getNumElements(),
		const_cast<int*>( rowCut.row().getIndices() ), 
		const_cast<double*>( rowCut.row().getElements() ),
		sns, rhs, NULL);
  checkGRBerror( err, "GRBaddconstr", "applyRowCut" );
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
		int err;
		assert(getEnvironmentPtr() != NULL);
		
		std::string pn;
		getStrParam(OsiProbName, pn);
		err = GRBnewmodel(getEnvironmentPtr(), &lp_, const_cast<char*>(pn.c_str()), 0, NULL, NULL, NULL, NULL, NULL);
		checkGRBerror( err, "GRBcreateprob", "getMutableLpPtr" );
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
  	int err = GRBfreemodel(lp_);
  	checkGRBerror( err, "GRBfreeprob", "gutsOfDestructor" );
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

  int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror( rc, "GRBupdatemodel", "basisIsAvailable" );

  int status;
  rc = GRBgetintattr(getMutableLpPtr(), GRB_INT_ATTR_STATUS, &status);
	checkGRBerror( rc, "GRBgetintattr", "basisIsAvailable" );

	if (status == GRB_LOADED || status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD || status == GRB_UNBOUNDED)
		return false;
	
  int dum;
  rc = GRBgetintattrelement(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, &dum);
  return (rc == 0);
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

	int rc = GRBupdatemodel(getMutableLpPtr());
  checkGRBerror(rc, "GRBupdatemodel", "getBasisStatus");

	rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_VBASIS, 0, numcols, cstat);
	checkGRBerror( rc, "GRBgetintattrarray", "getBasisStatus" );
	
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

	rc = GRBgetintattrarray(getMutableLpPtr(), GRB_INT_ATTR_CBASIS, 0, numrows, rstat);
	checkGRBerror( rc, "GRBgetintattrarray", "getBasisStatus" );

	char sense;
	for (int i = 0; i < numrows; ++i)
		switch (rstat[i])
		{
			case GRB_BASIC:
				rstat[i] = 1;
				break;
			case GRB_NONBASIC_LOWER:
			case GRB_NONBASIC_UPPER:
				rc = GRBgetcharattrelement(getMutableLpPtr(), GRB_CHAR_ATTR_SENSE, i, &sense);
				checkGRBerror( rc, "GRBgetcharattrelement", "getBasisStatus" );
				rstat[i] = (sense == '>' ? 2 : 3);
				break;
		}
}
