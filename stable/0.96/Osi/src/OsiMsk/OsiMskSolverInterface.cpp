/* Osi interface for Mosek ver. 4.0
   Lower versions are not supported
 -----------------------------------------------------------------------------
  name:     OSI Interface for MOSEK
  author:   Bo Jensen, Mads Jepsen
            
            email: support@MOSEK.com 
  date:     22 Dec 2005
 -----------------------------------------------------------------------------
*/
   
#if defined(_MSC_VER)
// Turn off compiler warning about long names 
#  pragma warning(disable:4786)
#endif
#include "mosekdl.h"
#include <iostream>
#include <cassert>
#include <string>
#include <numeric>
#include "CoinError.hpp"
#include "OsiMskSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

//#define DEBUG

#ifdef DEBUG
#define debugMessage printf
#else
#define debugMessage if( false ) printf
#endif
  
// Choose algorithm to be default in initial solve
// Only one of these flags should be set
//#define INITIAL_SOLVE MSK_OPTIMIZER_FREE_SIMPLEX
//#define INITIAL_SOLVE MSK_OPTIMIZER_DUAL_SIMPLEX
//#define INITIAL_SOLVE MSK_OPTIMIZER_PRIMAL_SIMPLEX
#define INITIAL_SOLVE MSK_OPTIMIZER_INTPNT

//Choose algorithm to be default in resolve
#define SOLVER MSK_OPTIMIZER_FREE_SIMPLEX
//#define SOLVER MSK_OPTIMIZER_INTPNT 
//#define SOLVER MSK_OPTIMIZER_DUAL_SIMPLEX
//#define INITIAL_SOLVE MSK_OPTIMIZER_PRIMAL_SIMPLEX
// Unset this flag to disable the warnings from interface.

//#define MSK_WARNING_ON

#undef getc

//#############################################################################
// A couple of helper functions
//#############################################################################

// Free memory pointet to by a double pointer

inline void freeCacheDouble( double*& ptr )
{
  if( ptr != NULL )
  {
    delete [] ptr;
    ptr = NULL;
  }
}

// Free memory pointet to by a char pointer

inline void freeCacheChar( char*& ptr )
{
  if( ptr != NULL )
  {
    delete [] ptr;
    ptr = NULL;
  }
}

// Free memory pointet to by a CoinPackedMatrix pointer 

inline void freeCacheMatrix( CoinPackedMatrix*& ptr )
{
  if( ptr != NULL )
  {
    delete ptr; 
    ptr = NULL;
  }
}

// Function used to connect MOSEK to stream
 
static void MSKAPI printlog(void *ptr,
                     char s[])
{
  printf("%s",s);
} 

// Prints a error message and throws a exception

static inline void
checkMSKerror( int err, std::string mskfuncname, std::string osimethod )
{
  if( err != MSK_RES_OK )
  {
    char s[100];
    sprintf( s, "%s returned error %d", mskfuncname.c_str(), err );
    std::cout << "ERROR: " << s << " (" << osimethod << 
	" in OsiMskSolverInterface)" << std::endl;
    throw CoinError( s, osimethod.c_str(), "OsiMskSolverInterface" );
  }
}

// Prints a warning message, can be shut off by undefining MSK_WARNING_ON

static inline void
OsiMSK_warning(std::string osimethod,  std::string warning)
{
   std::cout << "OsiMsk_warning: "<<warning<<" in "<<osimethod<< std::endl;
}

// Converts Range/Sense/Rhs to MSK bound structure

static inline void
MskConvertSenseToBound(const char rowsen, const double rowrng, 
					   const double rowrhs, double &rlb, double &rub, int &rtag)
{
    switch (rowsen) 
	{
	  case 'E':
	    rlb  = rub = rowrhs;
        rtag = MSK_BK_FX;
		break;
      case 'L':
        rlb  = -MSK_INFINITY;
        rub  = rowrhs;
        rtag = MSK_BK_UP;
        break;
      case 'G':
        rlb  = rowrhs;
        rub  = MSK_INFINITY;
        rtag = MSK_BK_LO;
	    break;
      case 'R':
		if( rowrng >= 0 )
		{
		  rlb = rowrhs - rowrng;
          rub = rowrhs;
		}
		else
		{
          rlb = rowrhs;
          rub = rowrhs + rowrng;
		}
        rtag = MSK_BK_RA;
        break;
      case 'N':
        rlb = -MSK_INFINITY;
        rub = MSK_INFINITY;
        rtag = MSK_BK_FR;
        break;
	  } 
}

// Converts a set of bounds to MSK boundkeys
 
static inline void
MskConvertColBoundToTag(const double collb, const double colub, double &clb, double &cub, int &ctag)
{
	if(collb > -MSK_INFINITY && colub < MSK_INFINITY)
	{
		ctag = MSK_BK_RA;
		clb = collb;
		cub = colub;
	}
	else if(collb <= - MSK_INFINITY && colub < MSK_INFINITY)
	{
		ctag = MSK_BK_UP;
		clb = -MSK_INFINITY;
		cub = colub;
	}
	else if(collb > - MSK_INFINITY && colub >= MSK_INFINITY)
	{
		ctag = MSK_BK_LO;
		clb = collb;
		cub = MSK_INFINITY;
	}
	else if(collb <= -MSK_INFINITY && colub >= MSK_INFINITY)
	{
		ctag = MSK_BK_FR;
		clb = -MSK_INFINITY;
		cub = MSK_INFINITY;
	} 
}

// Returns true if "solution" is defined in MOSEK, where solution can be basic, interior or 
// integer resp. (MSK_SOL_BAS), (MSK_SOL_ITR) or (MSK_SOL_ITG).

bool OsiMskSolverInterface::definedSolution(int solution) const
{
   int err, res;
   err = MSK_solutiondef(getMutableLpPtr(), solution, &res);
   checkMSKerror(err,"MSK_solutiondef","definedSolution");

   return ( res != MSK_RES_OK);
}

// Returns the flag for solver currently switched on in MOSEK resp. (MSK_OPTIMIZER_FREE),
// (MSK_OPTIMIZER_INTPNT), (MSK_OPTIMIZER_PRIMAL_SIMPLEX) or (MSK_OPTIMIZER_MIXED_INT).
// MOSEK also has Conic and nonconvex solvers, but these are for obvious reasons not 
// an option in the Osi interface.

int OsiMskSolverInterface::solverUsed() const
{
   int err, res;
   err = MSK_getintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, &res);
   checkMSKerror(err,"MSK_getintparam","definedSolution");
   return res;
}

// Switch the MOSEK solver to LP uses default solver specified by 'InitialSolver' 

void
OsiMskSolverInterface::switchToLP( void )
{
  debugMessage("OsiMskSolverInterface::switchToLP()\n");
  int err = MSK_putintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, InitialSolver);
  checkMSKerror(err,"MSK_putintparam","switchToLP");
  probtypemip_ = false;
}

// Switch the MOSEK solver to MIP. 

void
OsiMskSolverInterface::switchToMIP( void )
{
  debugMessage("OsiMskSolverInterface::switchToMIP()\n");
  int err = MSK_putintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_MIXED_INT);
  checkMSKerror(err,"MSK_putintparam","switchToMIP");
  probtypemip_ = true;
}

// Resize the coltype array. 

void
OsiMskSolverInterface::resizeColType( int minsize )
{
  debugMessage("OsiMskSolverInterface::resizeColType()\n");
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

// Free coltype array

void
OsiMskSolverInterface::freeColType()
{
  debugMessage("OsiMskSolverInterface::freeColType()\n");
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

//-----------------------------------------------------------------------------
// Free cached results and optimize the LP problem in task 

void OsiMskSolverInterface::initialSolve()
{
  debugMessage("OsiMskSolverInterface::initialSolve()\n");
  if( probtypemip_ == true )  
       switchToLP();
  else
  {
    if( definedSolution( MSK_SOL_BAS ) == true )
      {
        // Since the dual solver is not released yet, we choose primal simplex
        int err = MSK_putintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, SOLVER);
        checkMSKerror(err,"MSK_putintparam","initialSolve");
      }
    else
      {
        // No reoptimize possible use interior
        int err = MSK_putintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
        checkMSKerror(err,"MSK_putintparam","initialSolve");
      }
  }
  Mskerr = MSK_optimize(getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ));
}

//-----------------------------------------------------------------------------
// Resolves an LP problem. In MOSEK this is done automaticly in "optimize", so
// actually it makes no difference to call the above function or this
 
void OsiMskSolverInterface::resolve()
{
  debugMessage("OsiMskSolverInterface::resolve\n");
  
  if( probtypemip_ == true )  
    switchToLP();

  // We use free simplex since ??? 
  int err = MSK_putintparam(getMutableLpPtr(), MSK_IPAR_OPTIMIZER, SOLVER);
  checkMSKerror(err,"MSK_putintparam","initialSolve");

  Mskerr = MSK_optimize(getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ));
  MSK_solutionsummary(getMutableLpPtr(),0);
}

//-----------------------------------------------------------------------------
// Resolves an MIP problem with MOSEK MIP solver.

void OsiMskSolverInterface::branchAndBound()
{
  debugMessage("OsiMskSolverInterface::branchAndBound()\n");
  switchToMIP();
  Mskerr = MSK_optimize(getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ));
}

//#############################################################################
// Parameter related methods
//#############################################################################

//-----------------------------------------------------------------------------
// Sets a int parameter in MOSEK, had to make some logical assumptions in this function,
// see below.

bool
OsiMskSolverInterface::setIntParam(OsiIntParam key, int value)
{
  debugMessage("OsiMskSolverInterface::setIntParam(%d, %d)\n", key, value);
  bool retval = false;
  int solver;

  switch (key)
  {
      case OsiMaxNumIteration:
        // Check which solver is currently switched on change the max iteration
        // number allowed specific to that solver. 
        // Note that the iteration count is set to 0, every time you call "optimize".
        // So OsiMaxNumIteration=100, gives 100 iterations for every call to "optimize",
        // and not 100 in total since first time you invoked the optimizer.
        solver = solverUsed();
	    if( solver == MSK_OPTIMIZER_INTPNT )
          retval = (MSK_putintparam(
               getMutableLpPtr(), 
               MSK_IPAR_INTPNT_MAX_ITERATIONS, 
              value
          ) == MSK_RES_OK);
	   else if( solver == MSK_OPTIMIZER_PRIMAL_SIMPLEX 
             || solver == MSK_OPTIMIZER_DUAL_SIMPLEX 
	     || solver == MSK_OPTIMIZER_FREE_SIMPLEX )
          retval = (MSK_putintparam(
              getMutableLpPtr(), 
              MSK_IPAR_SIM_MAX_ITERATIONS, 
              value
          ) == MSK_RES_OK);
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
// Sets a double parameter in MOSEK.

bool
OsiMskSolverInterface::setDblParam(OsiDblParam key, double value)
{
  debugMessage("OsiMskSolverInterface::setDblParam(%d, %g)\n", key, value);
  bool retval = false;
  return true;
  switch (key) 
  {
    case OsiDualObjectiveLimit:
      if( getObjSense() == +1 )
	   retval = ( MSK_putdouparam( 
                       getMutableLpPtr(), 
                       MSK_DPAR_UPPER_OBJ_CUT, 
                       value 
       ) == MSK_RES_OK ); // min
      else
	   retval = ( MSK_putdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_LOWER_OBJ_CUT, 
                      value 
       ) == MSK_RES_OK ); // max
      break;
    case OsiPrimalObjectiveLimit:
      if( getObjSense() == +1 )
	   retval = ( MSK_putdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_LOWER_OBJ_CUT, 
                      value 
       ) == MSK_RES_OK ); // min
      else
	   retval = ( MSK_putdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_UPPER_OBJ_CUT, 
                      value 
       ) == MSK_RES_OK ); // max
      break;
    case OsiDualTolerance:
      retval = ( MSK_putdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_BASIS_TOL_S, 
                      value 
       ) == MSK_RES_OK );
      break;
	case OsiPrimalTolerance:
      retval = ( MSK_putdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_BASIS_TOL_X, 
                      value 
       ) == MSK_RES_OK );
      break;       
    case OsiObjOffset:
      retval = OsiSolverInterface::setDblParam(key, value);
      break;
    case OsiLastDblParam:
      retval = false;
      break;
   }

  return retval;
}


//-----------------------------------------------------------------------------
// Sets a string parameter in MOSEK.

bool
OsiMskSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  debugMessage("OsiMskSolverInterface::setStrParam(%d, %s)\n", key, value.c_str());

  bool retval=false;
  switch (key) 
  {
    case OsiProbName:
      OsiSolverInterface::setStrParam(key,value);
      return retval = true;
    case OsiSolverName:
      return false;
    case OsiLastStrParam:
      return false;
  }

  return false;
}

//-----------------------------------------------------------------------------
// Gets a int parameter in MOSEK. Same problems as in set int parameter, what 
// is max number iterations ? Simplex or Interior => choose the one switched on.

bool
OsiMskSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  debugMessage("OsiMskSolverInterface::getIntParam(%d)\n", key);

  bool retval = false;
  int solver;
  switch (key)
    {
    case OsiMaxNumIteration:
       solver = solverUsed();
	   if( solver == MSK_OPTIMIZER_INTPNT )
          retval = (MSK_getintparam(
                        getMutableLpPtr(), 
                        MSK_IPAR_INTPNT_MAX_ITERATIONS,  
                        &value
         ) == MSK_RES_OK);
	   else if( solver == MSK_OPTIMIZER_PRIMAL_SIMPLEX 
             || solver == MSK_OPTIMIZER_DUAL_SIMPLEX 
	     || solver == MSK_OPTIMIZER_DUAL_SIMPLEX  )
          retval = (MSK_getintparam(
                        getMutableLpPtr(), 
                        MSK_IPAR_SIM_MAX_ITERATIONS, 
                        &value
        ) == MSK_RES_OK);
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
// Get a double parameter in MOSEK.

bool
OsiMskSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  debugMessage("OsiMskSolverInterface::getDblParam(%d)\n", key);
  
  bool retval = false;

  switch (key) 
    {
    case OsiDualObjectiveLimit:
      if( getObjSense() == +1 )
	   retval = ( MSK_getdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_LOWER_OBJ_CUT, 
                      &value 
        ) == MSK_RES_OK ); // min
      else
	   retval = ( MSK_getdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_UPPER_OBJ_CUT, 
                      &value 
       ) == MSK_RES_OK ); // max
      break;
    case OsiPrimalObjectiveLimit:
      if( getObjSense() == +1 )
	   retval = ( MSK_getdouparam(  
                      getMutableLpPtr(), 
                      MSK_DPAR_UPPER_OBJ_CUT, 
                      &value 
       ) == MSK_RES_OK ); // min
      else
	   retval = ( MSK_getdouparam(  
                      getMutableLpPtr(), 
                      MSK_DPAR_LOWER_OBJ_CUT, 
                      &value 
       ) == MSK_RES_OK ); // max
      break;
    case OsiDualTolerance:
      retval = ( MSK_getdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_BASIS_TOL_S, 
                      &value 
       ) == MSK_RES_OK );
      break;
	case OsiPrimalTolerance:
      retval = ( MSK_getdouparam( 
                      getMutableLpPtr(), 
                      MSK_DPAR_BASIS_TOL_X, 
                      &value 
       ) == MSK_RES_OK );
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
// Gets a string parameter from MOSEK

bool
OsiMskSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  debugMessage("OsiMskSolverInterface::getStrParam(%d)\n", key);

  switch (key) 
  {
    case OsiProbName:
      OsiSolverInterface::getStrParam(key, value);
      break;
    case OsiSolverName:
      value = "MOSEK";
      break;
    case OsiLastStrParam:
      return false;
  }

  return true;
}

//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

//-----------------------------------------------------------------------------
// Returns true if solver abandoned in last call to solver.
// Mosek does not use this functionality

bool OsiMskSolverInterface::isAbandoned() const
{
  debugMessage("OsiMskSolverInterface::isAbandoned()\n");
  debugMessage("isAbandoned() %d\n",(Mskerr != MSK_RES_OK));
  return 0;
  //return (Mskerr != MSK_RES_OK);
}

//-----------------------------------------------------------------------------
// Returns true if "solution" available is proved to be optimal, where "solution" in LP
// could be both interior and basic, checks for both. 

bool OsiMskSolverInterface::isProvenOptimal() const
{
  debugMessage("OsiMskSolverInterface::isProvenOptimal()\n");
  int err, status, solution;
  if( probtypemip_ == false)
  {
	if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	else
	{
		if( definedSolution( MSK_SOL_ITR ) == true )
			solution = MSK_SOL_ITR;
		else
			return false;	
	}
  }
  else
  {
         if( definedSolution( MSK_SOL_ITG ) == true )
		solution = MSK_SOL_ITG;
	 else
		return false;     
  }
  err = MSK_getsolution(
	  getMutableLpPtr(),
	  solution,
	  NULL, 
	  &status, 
	  NULL, 
	  NULL, 
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL
	  );

  checkMSKerror(err,"MSK_getsolution","isProvenOptimal");
  debugMessage("Solution type , %d \n ",(status == MSK_SOL_STA_OPTIMAL)) ;
  return ( status == MSK_SOL_STA_OPTIMAL ); 	
}

//-----------------------------------------------------------------------------
// Returns true if a certificate of primal inf. exits

bool OsiMskSolverInterface::isProvenPrimalInfeasible() const
{
  debugMessage("OsiMskSolverInterface::isProvenPrimalInfeasible()\n");
  int err, status, solution;
  if( probtypemip_ == false)
  {
	if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	else
	{
		if( definedSolution( MSK_SOL_ITR ) == true )
			solution = MSK_SOL_ITR;
		else
			return false;	
	}
  }
  else
  {
	 if(  definedSolution( MSK_SOL_ITG ) == true)
		solution = MSK_SOL_ITG;
	 else
		return false;     
  }

  err = MSK_getsolution(
	  getMutableLpPtr(),
	  solution,
	  NULL, 
	  &status, 
	  NULL, 
	  NULL, 
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL
	  );

  checkMSKerror(err,"MSK_getsolution","isProvenPrimalInfeasible");
  debugMessage("isProvenPrimalInfeasible %d \n",(status == MSK_SOL_STA_PRIM_INFEAS_CER));
  return ( status == MSK_SOL_STA_PRIM_INFEAS_CER ); 	
}

//-----------------------------------------------------------------------------
// Should return true if a certificate of dual inf. exits
// But COIN does not support this feature thus we return false

bool OsiMskSolverInterface::isProvenDualInfeasible() const
{
  debugMessage("OsiMskSolverInterface::isProvenDualInfeasible()\n");

  int err, status, solution;
  if( probtypemip_ == false)
  {
	if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	else
	{
		if( definedSolution( MSK_SOL_ITR ) == true )
			solution = MSK_SOL_ITR;
		else
			return false;	
	}
  }
  else
  {
	 if(  definedSolution( MSK_SOL_ITG ) == true)
		solution = MSK_SOL_ITG;
	 else
		return false;     
  }

  err = MSK_getsolution(
	  getMutableLpPtr(),
	  solution,
	  NULL, 
	  &status, 
	  NULL, 
	  NULL, 
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL
	  );

  checkMSKerror(err,"MSK_getsolution","isProvenDualInfeasible");
  debugMessage("isProvenDualInfeasible %d \n",(status == MSK_SOL_STA_DUAL_INFEAS_CER));
  return ( status == MSK_SOL_STA_DUAL_INFEAS_CER ); 	
}

//-----------------------------------------------------------------------------
// Returns true if primal objective limit is reached. Checks the objective sense
// first. 

bool OsiMskSolverInterface::isPrimalObjectiveLimitReached() const
{
  debugMessage("OsiMskSolverInterface::isPrimalObjectiveLimitReached()\n");
  int err;
  double obj = getObjValue(),value;
  if( getObjSense() == +1 )
    {
	 err = MSK_getdouparam( 
				 getMutableLpPtr(), 
				 MSK_DPAR_UPPER_OBJ_CUT, 
				 &value);
    } 
  else
  {
	 err = MSK_getdouparam( 
				 getMutableLpPtr(), 
				 MSK_DPAR_LOWER_OBJ_CUT, 
				 &value 
				 );
	 obj = -obj;
	 value = -value;
  }
  
  checkMSKerror( err, "MSK_getdouparam", "isPrimalObjectiveLimitReached" );
  debugMessage("primal objective value  %f , lowerbound %f , reached %i \n",obj,value,(value <= obj));
  return ( value <= obj ); 	
}

//-----------------------------------------------------------------------------
// Returns true if dual objective limit is reached. Checks the objective sense
// first.  

bool OsiMskSolverInterface::isDualObjectiveLimitReached() const
{
  debugMessage("OsiMskSolverInterface::isDualObjectiveLimitReached()\n");
  int err;
  double obj = getObjValue(),value;
  if( getObjSense() == +1 )
    {
	 err = MSK_getdouparam( 
				 getMutableLpPtr(), 
				 MSK_DPAR_UPPER_OBJ_CUT, 
				 &value); 
    }
  else
  {
	 err = MSK_getdouparam( 
				 getMutableLpPtr(), 
				 MSK_DPAR_LOWER_OBJ_CUT, 
				 &value 
    );
    obj = -obj;
    value = -value;
  }
  checkMSKerror( err, "MSK_getdouparam", "isPrimalObjectiveLimitReached" );
  debugMessage("dual objective value  %f , lowerbound %f , reached %i \n",obj,value,(value <= obj));
  //debugMessage("dual obj reached %i \n",(obj <= value));
  return ( value <= obj ); 	
}

//-----------------------------------------------------------------------------
// Returns true if iteration number used in last call to optimize eq. max number for
// the solver used. 

bool OsiMskSolverInterface::isIterationLimitReached() const
{
  debugMessage("OsiMskSolverInterface::isIterationLimitReached()\n");
  int solver = solverUsed(), err = MSK_RES_OK, value = 0, iter;
  iter = getIterationCount();
  if( solver == MSK_OPTIMIZER_INTPNT )
     err = MSK_getintparam(
                getMutableLpPtr(), 
                MSK_IPAR_INTPNT_MAX_ITERATIONS, 
                &value
    );
  else if( solver == MSK_OPTIMIZER_PRIMAL_SIMPLEX ||
	   solver == MSK_OPTIMIZER_DUAL_SIMPLEX   ||
	   solver == MSK_OPTIMIZER_FREE_SIMPLEX 
	   )
     err = MSK_getintparam(
               getMutableLpPtr(), 
               MSK_IPAR_SIM_MAX_ITERATIONS, 
               &value
    );

  checkMSKerror( err, "MSK_getintparam", "isIterationLimitReached" );
  debugMessage("iteration limit reached %i \n",(value <= iter));
  return (value <= iter);
}

//#############################################################################
// WarmStart related methods
//#############################################################################

//-----------------------------------------------------------------------------
// Get warm start, returns NULL pointer if not availabel

CoinWarmStart* OsiMskSolverInterface::getWarmStart() const
{
  debugMessage("OsiMskSolverInterface::getWarmStart()\n");

  CoinWarmStartBasis* ws = NULL;
  int numcols = getNumCols();
  int numrows = getNumRows();
  int *cstat;
  int *rstat;
  int err, i;
  int *bk;
  bool skip = false;
  
  assert(!probtypemip_);
  if( definedSolution( MSK_SOL_BAS ) == true )
  {
      cstat = new int[numcols];
      rstat = new int[numrows];

      err  =   MSK_getsolution(getMutableLpPtr(),
	                             MSK_SOL_BAS,
                                 NULL,
                                 NULL,
	                             rstat, 
	                             cstat,
	                             NULL,
	                             NULL,
	                             NULL,
	                             NULL,
	                             NULL,
	                             NULL,
	                             NULL,
	                             NULL,
                                 NULL
      );
      
      checkMSKerror( err, "MSK_getsolution", "getWarmStart" );
  }
  else
  {
      /* No basic solution stored choose slack basis */
      /* Otherwise the unittest can not be passed    */
      ws = new CoinWarmStartBasis;
      ws->setSize( numcols, numrows );
      
      bk = new int[numcols];
      
      for( i = 0; i < numrows; ++i )
	      ws->setArtifStatus( i, CoinWarmStartBasis::basic );

      err = MSK_getboundslice(getMutableLpPtr(),
                                0,
                                0,
                                numcols,
                                bk,
                                NULL,
                                NULL);

      checkMSKerror( err, "MSK_getsolution", "getWarmStart" );
                          
      for( i = 0; i < numcols; ++i )
	  {
          switch(bk[i])
          {
                  case MSK_BK_RA:
                  case MSK_BK_LO:
                     ws->setStructStatus( i, CoinWarmStartBasis::atLowerBound );
                  break;
                  case MSK_BK_FX:
                  case MSK_BK_UP:
                     ws->setStructStatus( i, CoinWarmStartBasis::atUpperBound );
                  break;
                  case MSK_BK_FR:
                     ws->setStructStatus( i, CoinWarmStartBasis::isFree );
                  break;
                  default:
                  checkMSKerror( 1, "Wrong bound key", "getWarmStart" );
          }
	  }
      
      delete[] bk;
      
      return ws;
  }
  
  if( err == MSK_RES_OK )
  {
    ws = new CoinWarmStartBasis;
    ws->setSize( numcols, numrows );
    
    for( i = 0; i < numrows; ++i )
      {
        switch( rstat[i] )
          {
          case MSK_SK_BAS:
            ws->setArtifStatus( i, CoinWarmStartBasis::basic );
            break;
          case MSK_SK_LOW:
            ws->setArtifStatus( i, CoinWarmStartBasis::atLowerBound );
            break;
          case MSK_SK_FIX:
          case MSK_SK_UPR:
            ws->setArtifStatus( i, CoinWarmStartBasis::atUpperBound );
            break;
          case MSK_SK_SUPBAS:
            ws->setArtifStatus( i, CoinWarmStartBasis::isFree );
            break;
          default:  // unknown row status
            delete ws;
            ws   = NULL;
            skip = true;
            break;
          }
      }
    if( skip == false )
      {
        for( i = 0; i < numcols; ++i )
          {
            switch( cstat[i] )
              {
              case MSK_SK_BAS:
                ws->setStructStatus( i, CoinWarmStartBasis::basic );
                break;
              case MSK_SK_LOW:
                ws->setStructStatus( i, CoinWarmStartBasis::atLowerBound );
                break;
              case MSK_SK_FIX:
              case MSK_SK_UPR:
                ws->setStructStatus( i, CoinWarmStartBasis::atUpperBound );
                break;
              case MSK_SK_SUPBAS:
                ws->setStructStatus( i, CoinWarmStartBasis::isFree );
                break;
              default:  // unknown column status
                delete ws;
                ws = NULL;
                break;
              }
          }
      }
 }

 delete[] cstat;
 delete[] rstat;

 return ws;
}

//-----------------------------------------------------------------------------
// Set warm start

bool OsiMskSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  debugMessage("OsiMskSolverInterface::setWarmStart(%p)\n", warmstart);

  const CoinWarmStartBasis* ws = dynamic_cast<const CoinWarmStartBasis*>(warmstart);
  int numcols, numrows, i, restat;
  int *cstat, *rstat, *bkc, *bkx;
  bool retval = false, skip = false;

  if( !ws )
    return false;

  numcols = ws->getNumStructural();
  numrows = ws->getNumArtificial();
  
  if( numcols != getNumCols() || numrows != getNumRows() )
    return false;

  switchToLP();

  cstat = new int[numcols];
  rstat = new int[numrows];
  bkc   = new int[numrows];
      
  restat = MSK_getboundslice(getMutableLpPtr(),
                               1,
                               0,
                               numrows,
                               bkc,
                               NULL,
                               NULL);

  checkMSKerror( restat, "MSK_getboundslice", "setWarmStart" );

  for( i = 0; i < numrows; ++i )
  {
      switch( ws->getArtifStatus( i ) )
	  {
	    case CoinWarmStartBasis::basic:
	      rstat[i] = MSK_SK_BAS;
	    break;
	    case CoinWarmStartBasis::atLowerBound:
            switch(bkc[i])
            {
                    case MSK_BK_LO:
                    case MSK_BK_RA:
	                  rstat[i] = MSK_SK_LOW;
                    break;
                    case MSK_BK_UP:
                    case MSK_BK_FR:
	                  rstat[i] = MSK_SK_SUPBAS;
                    break;
                    case MSK_BK_FX:
	                  rstat[i] = MSK_SK_FIX;
                    break;
            }        
	    break;
	    case CoinWarmStartBasis::atUpperBound:
            switch(bkc[i])
            {
                    case MSK_BK_UP:
                    case MSK_BK_RA:
	                  rstat[i] = MSK_SK_UPR;
                    break;
                    case MSK_BK_LO:
                    case MSK_BK_FR:
	                  rstat[i] = MSK_SK_SUPBAS;
                    break;
                    case MSK_BK_FX:
	                  rstat[i] = MSK_SK_FIX;
                    break;
            }        
        break;
        case CoinWarmStartBasis::isFree:
	               rstat[i] = MSK_SK_SUPBAS;
	    break;
	    default:  // unknown row status
	     retval = false;
         skip   = true;
        break;
	  }
  }

  delete[] bkc;
  
  if( skip == false )
  {
          bkx = new int[numcols];
      
          restat = MSK_getboundslice(getMutableLpPtr(),
                                       0,
                                       0,
                                       numcols,
                                       bkx,
                                       NULL,
                                       NULL);

         checkMSKerror( restat, "MSK_getboundslice", "setWarmStart" );

         for( i = 0; i < numcols; ++i )
          {
            switch( ws->getStructStatus( i ) )
            {
              case CoinWarmStartBasis::basic:
               cstat[i] = MSK_SK_BAS;
               break;
              case CoinWarmStartBasis::atLowerBound:
                    switch(bkx[i])
                    {
                            case MSK_BK_LO:
                            case MSK_BK_RA:
                              cstat[i] = MSK_SK_LOW;
                            break;
                            case MSK_BK_UP:
                            case MSK_BK_FR:
                              cstat[i] = MSK_SK_SUPBAS;
                            break;
                            case MSK_BK_FX:
                              cstat[i] = MSK_SK_FIX;
                            break;
                    }        
               break;
             case CoinWarmStartBasis::atUpperBound:
                    switch(bkx[i])
                    {
                            case MSK_BK_UP:
                            case MSK_BK_RA:
                              cstat[i] = MSK_SK_UPR;
                            break;
                            case MSK_BK_LO:
                            case MSK_BK_FR:
                              cstat[i] = MSK_SK_SUPBAS;
                            break;
                            case MSK_BK_FX:
                              cstat[i] = MSK_SK_FIX;
                            break;
                    }        
               break;
             case CoinWarmStartBasis::isFree:
               cstat[i] = MSK_SK_SUPBAS;
               break;
             default:  // unknown col status
               retval   = false;
               skip     = true;
             break;
            }
          }
          
          delete[] bkx;
  }

  if( skip == false )
  {
          restat = MSK_putsolution( 
              getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ),
              MSK_SOL_BAS,
              rstat, 
              cstat,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL
          );
          
         delete[] cstat;
         delete[] rstat;
  }
  else
     return false;
  
  retval = (restat == MSK_RES_OK);
 
  return retval;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

//-----------------------------------------------------------------------------
// Mark hot start

void OsiMskSolverInterface::markHotStart()
{
  debugMessage("OsiMskSolverInterface::markHotStart()\n");

  int err;
  int numcols, numrows;

  assert(!probtypemip_);

  numcols = getNumCols();
  numrows = getNumRows();
  
  if( numcols > hotStartCStatSize_ )
  {
      delete[] hotStartCStat_;
      hotStartCStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numcols ) ); 
      // get some extra space for future hot starts
      hotStartCStat_ = new int[hotStartCStatSize_];
  }
    
  if( numrows > hotStartRStatSize_ )
  {
      delete[] hotStartRStat_;
      hotStartRStatSize_ = static_cast<int>( 1.2 * static_cast<double>( numrows ) ); 
      // get some extra space for future hot starts
      hotStartRStat_ = new int[hotStartRStatSize_];
  }

  err = MSK_getsolution( 	  
	  getMutableLpPtr(),
	  MSK_SOL_BAS,
	  NULL,
	  NULL,
	  hotStartRStat_, 
	  hotStartCStat_,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL
  );
  
  checkMSKerror( err, "MSK_getsolution", "markHotStart" );
}

//-----------------------------------------------------------------------------
// Solve from a hot start

void OsiMskSolverInterface::solveFromHotStart()
{
  debugMessage("OsiMskSolverInterface::solveFromHotStart()\n");

  int err;
  int maxiter;

  switchToLP();

  assert( getNumCols() <= hotStartCStatSize_ );
  assert( getNumRows() <= hotStartRStatSize_ );

  err = MSK_putsolution( 
	  getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ),
	  MSK_SOL_BAS,
	  hotStartRStat_, 
	  hotStartCStat_,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL,
	  NULL
  );
  
  checkMSKerror( err, "MSK_putsolution", "solveFromHotStart" );
  MSKtask_t task = getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS );

  err = MSK_getintparam( task, MSK_IPAR_SIM_MAX_ITERATIONS, &maxiter );
  checkMSKerror( err, "MSK_getintparam", "solveFromHotStart" );

  err = MSK_putintparam( task, MSK_IPAR_SIM_MAX_ITERATIONS, hotStartMaxIteration_ );
  checkMSKerror( err, "MSK_putintparam", "solveFromHotStart" );
  
  resolve();

  err = MSK_putintparam( task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter );
  checkMSKerror( err, "MSK_putintparam", "solveFromHotStart" );
}

//-----------------------------------------------------------------------------
// Unmark a hot start

void OsiMskSolverInterface::unmarkHotStart()
{
  debugMessage("OsiMskSolverInterface::unmarkHotStart()\n");
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################


//-----------------------------------------------------------------------------
// Returns number of columns in MOSEK task

int OsiMskSolverInterface::getNumCols() const
{
  debugMessage("OsiMskSolverInterface::getNumCols()\n");
  int numcol, err;
  err = MSK_getnumvar(getMutableLpPtr(),&numcol);
  checkMSKerror( err, "MSK_getnumvar", "getNumCols" );
  return numcol;
}

//-----------------------------------------------------------------------------
// Returns number of rows in MOSEK task

int OsiMskSolverInterface::getNumRows() const
{
  debugMessage("OsiMskSolverInterface::getNumRows()\n");
  int numrow, err;
  err = MSK_getnumcon(getMutableLpPtr(),&numrow);
  checkMSKerror( err, "MSK_getnumcon", "getNumRows" );
  return numrow;
}

//-----------------------------------------------------------------------------
// Returns number of non-zeroes (in matrix) in MOSEK task

int OsiMskSolverInterface::getNumElements() const
{
  debugMessage("OsiMskSolverInterface::getNumElements()\n");
  int numnon, err;
  err = MSK_getnumanz(getMutableLpPtr(),&numnon);
  checkMSKerror( err, "MSK_getnumanz", "getNumElements" );
  return numnon;
}


//-----------------------------------------------------------------------------
// Returns lower bounds on columns in MOSEK task

const double * OsiMskSolverInterface::getColLower() const
{
  debugMessage("OsiMskSolverInterface::getColLower()\n");

  if( collower_ == NULL )
  {
     assert(colupper_ == NULL);
	 int ncols = getNumCols();
     if( ncols > 0 )
	 {
	   colupper_       = new double[ncols];
	   collower_       = new double[ncols];
	   int *dummy_tags = new int[ncols];
	   int err = MSK_getboundslice( 
		 getMutableLpPtr(), 
		 0, 
		 0, 
		 ncols, 
		 dummy_tags, 
		 collower_, 
		 colupper_
       );
	   checkMSKerror(err, "MSK_getboundslice","getColUpper");
	   delete[] dummy_tags;
	 }
  }

  return collower_;
}

//-----------------------------------------------------------------------------
// Returns upper bounds on columns in MOSEK task

const double * OsiMskSolverInterface::getColUpper() const
{
  debugMessage("OsiMskSolverInterface::getColUpper()\n");

  if( colupper_ == NULL )
  {
	 assert( collower_ == NULL );
	 int ncols = getNumCols();
     if( ncols > 0 )
	 {
	   colupper_       = new double[ncols];
	   collower_       = new double[ncols];
	   int *dummy_tags = new int[ncols];
	   int err = MSK_getboundslice( 
		 getMutableLpPtr(), 
		 0, 
		 0, 
		 ncols, 
		 dummy_tags, 
		 collower_, 
		 colupper_
       );
	   checkMSKerror(err,"MSK_getboundslice","getColUpper");
	   delete[] dummy_tags;
	 }
  }

  return colupper_;
}


//-----------------------------------------------------------------------------
// Returns rowsense in MOSEK task, call getRightHandSide to produce triplets.

const char * OsiMskSolverInterface::getRowSense() const
{
  debugMessage("OsiMskSolverInterface::getRowSense()\n");

  if( rowsense_==NULL )
  {      
    getRightHandSide();
    assert( rowsense_!=NULL || getNumRows() == 0 );
  }

  return rowsense_;
}

//-----------------------------------------------------------------------------
// Returns the RHS in triplet form. MOSEK uses always boundkeys instead of the
// triplet, so we have to convert back to triplet. 

const double * OsiMskSolverInterface::getRightHandSide() const
{
  debugMessage("OsiMskSolverInterface::getRightHandSide()\n");

  if(rowsense_ == NULL) 
  {
    assert ((rhs_ == NULL) && (rowrange_ == NULL));      
    int nr = getNumRows();
    if ( nr != 0 ) 
	{
      rowsense_         = new char[nr];
      rhs_              = new double[nr];
      rowrange_         = new double[nr]; 
      const double * lb = getRowLower();
      const double * ub = getRowUpper();      
      int i;
      
      for ( i=0; i<nr; i++ )
        convertBoundToSense(lb[i], ub[i], rowsense_[i], rhs_[i], rowrange_[i]);
    }
  }

  return rhs_;
}

//-----------------------------------------------------------------------------
// Returns rowrange in MOSEK task, call getRightHandSide to produce triplets.

const double * OsiMskSolverInterface::getRowRange() const
{
  debugMessage("OsiMskSolverInterface::getRowRange()\n");

  if( rowrange_ == NULL ) 
  {
    getRightHandSide();
    assert( rowrange_!=NULL || getNumRows() == 0 );
  }

  return rowrange_;
}


//-----------------------------------------------------------------------------
// Returns lower bounds on rows in MOSEK task.

const double * OsiMskSolverInterface::getRowLower() const
{
  debugMessage("OsiMskSolverInterface::getRowLower()\n");

  if( rowlower_ == NULL )
  {
	  assert( rowupper_ == NULL );
      int nrows = getNumRows();

      if( nrows > 0 )
	  {
	    rowlower_       = new double[nrows];
		rowupper_       = new double[nrows];
	    int *dummy_tags = new int[nrows];
        
	    int err = MSK_getboundslice(
		  getMutableLpPtr(), 
		  1, 
		  0, 
		  nrows, 
		  dummy_tags, 
		  rowlower_, 
		  rowupper_
		);
        
	    checkMSKerror(err,"MSK_getboundslice","getRowLower");
	    delete[] dummy_tags;
	  }
  }

  return rowlower_;
}


//-----------------------------------------------------------------------------
// Returns upper bounds on rows in MOSEK task.

const double * OsiMskSolverInterface::getRowUpper() const
{  
  debugMessage("OsiMskSolverInterface::getRowUpper()\n");

  if( rowupper_ == NULL )
  {
      assert( rowlower_ == NULL );
      int nrows = getNumRows();
      
      if( nrows > 0 )
	  {
	    rowupper_       = new double[nrows];
		rowlower_       = new double[nrows];
	    int *dummy_tags = new int[nrows];
        
	    int err = MSK_getboundslice(
		  getMutableLpPtr(), 
		  1, 
		  0,
		  nrows, 
		  dummy_tags, 
		  rowlower_, 
		  rowupper_
		);
        
	    checkMSKerror(err,"MSK_getboundslice","getRowUpper");
	    delete[] dummy_tags;
	  }
  }

  return rowupper_;
}


//-----------------------------------------------------------------------------
// Returns objective coefficient in MOSEK task.


const double * OsiMskSolverInterface::getObjCoefficients() const
{
  debugMessage("OsiMskSolverInterface::getObjCoefficients()\n");

  if( obj_ == NULL )
  {
      int ncols = getNumCols();

      if( ncols > 0 )
	  {
	    obj_    = new double[ncols]; 
   	    int err = MSK_getc( getMutableLpPtr(), obj_ );
	    checkMSKerror( err, "MSK_getc", "getObjCoefficients" );
	  }
  }

  return obj_;
}


//-----------------------------------------------------------------------------
// Returns the direction of optimization

double OsiMskSolverInterface::getObjSense() const
{
  debugMessage("OsiMskSolverInterface::getObjSense()\n");
  int sense,err;
  
  err = MSK_getintparam(
              getMutableLpPtr(),
              MSK_IPAR_OBJECTIVE_SENSE,
              &sense
     );

  checkMSKerror(err,"MSK_getintparam","getObjSense");

  if( sense == MSK_OBJECTIVE_SENSE_MIN )
    return +1.0;
  else
    return -1.0;
}


//-----------------------------------------------------------------------------
// Returns true if variabel is set to continuous

bool OsiMskSolverInterface::isContinuous( int colNumber ) const
{
  debugMessage("OsiMskSolverInterface::isContinuous(%d)\n", colNumber);
  return getCtype()[colNumber] == 'C';
}


//-----------------------------------------------------------------------------
// Returns a Coin matrix by row

const CoinPackedMatrix * OsiMskSolverInterface::getMatrixByRow() const
{
  debugMessage("OsiMskSolverInterface::getMatrixByRow()\n");

  if ( matrixByRow_ == NULL ) 
  {
	int nc, nr, nz, *sub, *ptrb, *ptre, surp, *len;
	double *val;
	nc       = getNumCols();
	nr       = getNumRows();
	nz       = surp = getNumElements();
	ptrb     = new int[nr+1];
	ptre     = new int[nr];
    sub      = new int[nz];
	val      = new double[nz];
    len      = new int[nr];
    ptrb[nr] = nz;

	int err = MSK_getaslice(
		getMutableLpPtr(),
        1,
	    0,
	    nr,
	    nz,
	    &surp, 
	    ptrb, 
	    ptre, 
	    sub, 
	    val
	); 

	checkMSKerror(err, "MSK_getaslice", "getMatrixByRow");
	 
    for(int i=0; i<nr; i++)
		len[i] = ptre[i]-ptrb[i];

    matrixByRow_ = new CoinPackedMatrix();
    matrixByRow_->assignMatrix(false , nc, nr, nz, val, sub, ptrb, len);

    assert( matrixByRow_->getNumCols()==nc );
    assert( matrixByRow_->getNumRows()==nr );

	delete[] ptre;
  }

  debugMessage("End OsiMskSolverInterface::getMatrixByRow()\n");

  return matrixByRow_; 
} 


//-----------------------------------------------------------------------------
// Returns a Coin matrix by column

const CoinPackedMatrix * OsiMskSolverInterface::getMatrixByCol() const
{
  debugMessage("OsiMskSolverInterface::getMatrixByCol()\n");

  if ( matrixByCol_ == NULL ) 
  {
	int nc, nr, nz, *sub, *ptrb, *ptre, surp, *len;
	double *val;
	nc       = getNumCols();
	nr       = getNumRows();
	nz       = surp = getNumElements();
	ptrb     = new int[nc+1];
	ptre     = new int[nc];
    sub      = new int[nz];
	val      = new double[nz];
    len      = new int[nc];
    ptrb[nc] = nz;

	int err = MSK_getaslice(
		getMutableLpPtr(),
        0,
	    0,
	    nc,
	    nz,
	    &surp, 
	    ptrb, 
	    ptre, 
	    sub, 
	    val
	); 

	checkMSKerror(err, "MSK_getaslice", "getMatrixByCol");
	 
    for(int i=0; i<nc; i++)
		len[i] = ptre[i]-ptrb[i];

    matrixByCol_ = new CoinPackedMatrix();
    matrixByCol_->assignMatrix(true , nr, nc, nz, val, sub, ptrb, len);

    assert( matrixByCol_->getNumCols()==nc );
    assert( matrixByCol_->getNumRows()==nr );

	delete[] ptre;
  }

  return matrixByCol_; 
} 


//-----------------------------------------------------------------------------
// Returns the infinity level used in MOSEK.

double OsiMskSolverInterface::getInfinity() const
{
  debugMessage("OsiMskSolverInterface::getInfinity()\n");

  return MSK_INFINITY;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

//-----------------------------------------------------------------------------
// Returns the current col solution. 

const double * OsiMskSolverInterface::getColSolution() const
{
  debugMessage("OsiMskSolverInterface::getColSolution()\n");
  if( colsol_ != NULL )
    {
      debugMessage("colsol_ != NULL"); 
      delete[] colsol_;
    }
  int i;
  int nc = getNumCols();
  if( nc > 0 )
    {
      int solution = MSK_RES_ERR_UNDEF_SOLUTION;
      colsol_ = new double[nc];
      
      if( probtypemip_ == false)
	{
	  if( definedSolution( MSK_SOL_BAS ) == true )
	    solution = MSK_SOL_BAS;
	  else if( definedSolution( MSK_SOL_ITR) == true )
	    solution = MSK_SOL_ITR;
	}
      else if( definedSolution( MSK_SOL_ITG ) == true )
	solution = MSK_SOL_ITG;  
      
      if ( solution == MSK_RES_ERR_UNDEF_SOLUTION )
        {
	  for( i = 0; i < nc; ++i )
	    colsol_[i] = 0.0;
	  
	  return colsol_;
        }
      int err = MSK_getsolution( 	  
				  getMutableLpPtr(),
				  solution,
				  NULL,
				  NULL,
				  NULL, 
				  NULL,
				  NULL,
				  NULL,
				  colsol_,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL
				  );
      
      checkMSKerror(err,"MSK_getsolution","getColSolution");
    }
  return colsol_;
}

//-----------------------------------------------------------------------------
// Returns the row price / dual variabels in MOSEK task

const double * OsiMskSolverInterface::getRowPrice() const
{
  debugMessage("OsiMskSolverInterface::getRowPrice()\n");

  if( rowsol_ == NULL )
  {
      int i;
      int nr = getNumRows();
      if( nr > 0 )
	  {
		int solution = MSK_RES_ERR_UNDEF_SOLUTION;
	    rowsol_      = new double[nr];

		if( probtypemip_ == false)
		{
			if( definedSolution( MSK_SOL_BAS ) == true )
				solution = MSK_SOL_BAS;
			else if( definedSolution( MSK_SOL_ITR) == true )
				solution = MSK_SOL_ITR;
		 }
		 else if( definedSolution( MSK_SOL_ITG ) == true )
				solution = MSK_SOL_ITG;  
         
	     if ( solution == MSK_RES_ERR_UNDEF_SOLUTION )
         {
            for( i = 0; i < nr; ++i )
               rowsol_[i] = 0.0;

			return rowsol_;
         }

		int err = MSK_getsolution( 	  
				  getMutableLpPtr(),
				  solution,
				  NULL,
				  NULL,
				  NULL, 
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  rowsol_,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL
		);
	    checkMSKerror( err, "MSK_getsolution", "getRowPrice" );
	  }
   }

   return rowsol_;
}

//-----------------------------------------------------------------------------
// Returns the reduced cost in MOSEK task.

const double * OsiMskSolverInterface::getReducedCost() const
{
  debugMessage("OsiMskSolverInterface::getReducedCost()\n");

  if( redcost_ == NULL )
  {
      int ncols = getNumCols();

      if( ncols > 0 )
	  {
		int solution = MSK_RES_ERR_UNDEF_SOLUTION;
        if( probtypemip_ == false)
        {
            if( definedSolution( MSK_SOL_BAS ) == true )
				solution = MSK_SOL_BAS;
            else if( definedSolution( MSK_SOL_ITR) == true )
				solution = MSK_SOL_ITR;
		 }
         else if( definedSolution( MSK_SOL_ITG ) == true )
				solution = MSK_SOL_ITG;  
	     if ( solution == MSK_RES_ERR_UNDEF_SOLUTION )
			return NULL;
		 
	    redcost_    = new double[ncols];
		double *slx = new double[ncols];
		double *sux = new double[ncols];

        int err = MSK_getsolution( 	  
				  getMutableLpPtr(),
				  solution,
				  NULL,
				  NULL,
				  NULL, 
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  slx,
				  sux,
				  NULL
		);
        
        // Calculate reduced cost
		for(int i = 0; i < ncols; i++)
			redcost_[i] = slx[i]-sux[i];

		delete[] slx;
		delete[] sux;

	    checkMSKerror( err, "MSK_getsolution", "getReducedCost" );
	  }
  }

  return  redcost_;
}


//-----------------------------------------------------------------------------
// Returns the rowactivity in MOSEK task.

const double * OsiMskSolverInterface::getRowActivity() const
{
  debugMessage("OsiMskSolverInterface::getRowActivity()\n");
  if( rowact_ == NULL )
  {
      int i;
      int nrows = getNumRows();
      if( nrows > 0 )
	  {
		rowact_ = new double[nrows];
		int solution = MSK_RES_ERR_UNDEF_SOLUTION;
		if( probtypemip_ == false)
		{
			if( definedSolution( MSK_SOL_BAS ) == true )
				solution = MSK_SOL_BAS;
			else if( definedSolution( MSK_SOL_ITR) == true )
				solution = MSK_SOL_ITR;
		 }
		 else if( definedSolution( MSK_SOL_ITG ) == true )
				solution = MSK_SOL_ITG;
         
	     if ( solution == MSK_RES_ERR_UNDEF_SOLUTION )
         {
            for( i = 0; i < nrows; ++i )
               rowact_[i] = 0.0;

			return rowact_;
         }

		int err = MSK_getsolution( 	  
				  getMutableLpPtr(),
				  solution,
				  NULL,
				  NULL,
				  NULL, 
				  NULL,
				  NULL,
				  rowact_,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL
		);
        
	    checkMSKerror( err, "MSK_getsolution", "getRowActivity" );
	  }
 }

 return  rowact_;
}

//-----------------------------------------------------------------------------
// Returns the objective for defined solution in MOSEK task. 

double OsiMskSolverInterface::getObjValue() const
{
  debugMessage("OsiMskSolverInterface::getObjValue()\n");
  double value;
  int solution = MSK_RES_ERR_UNDEF_SOLUTION;
  if( probtypemip_ == false)
  {
	 if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	 else if( definedSolution( MSK_SOL_ITR) == true )
		solution = MSK_SOL_ITR;
  }
  else if( definedSolution( MSK_SOL_ITG ) == true )
	 solution = MSK_SOL_ITG;  

  if ( solution == MSK_RES_ERR_UNDEF_SOLUTION )
  {
     #if defined MSK_WARNING_ON
     OsiMSK_warning("OsiMskSolverInterface::getObjValue()",
     "Undefined solution but interface returns zero");
     #endif
     return 0.0;
  }

  int err = MSK_getprimalobj( getMutableLpPtr(), solution, &value );
  checkMSKerror(err,"MSK_getprimalobj","getObjValue");

  return value;
}

//-----------------------------------------------------------------------------
// Returns the iteration used in last call to optimize. Notice that the cross
// over phase in interior methods is not returned, when interior point is
// used, only the interior point iterations.

int OsiMskSolverInterface::getIterationCount() const
{
  debugMessage("OsiMskSolverInterface::getIterationCount()\n");

  int nr = 0, solver, err;
  int nrp=0;
  solver = solverUsed();
  if( solver == MSK_OPTIMIZER_PRIMAL_SIMPLEX )
    {
      err = MSK_getintinf(getMutableLpPtr(), MSK_IINF_SIM_PRIMAL_ITER, &nr);
      checkMSKerror(err,"MSK_getintinf","getIterationsCount");
    }
  
  if( solver == MSK_OPTIMIZER_DUAL_SIMPLEX  )
    {
      err = MSK_getintinf(getMutableLpPtr(), MSK_IINF_SIM_DUAL_ITER, &nr);
      checkMSKerror(err,"MSK_getintinf","getIterationsCount");
    }

  if( solver == MSK_OPTIMIZER_FREE_SIMPLEX  )
    {
      err = MSK_getintinf(getMutableLpPtr(), MSK_IINF_SIM_DUAL_ITER, &nr);
      checkMSKerror(err,"MSK_getintinf","getIterationsCount");
      err = MSK_getintinf(getMutableLpPtr(), MSK_IINF_SIM_PRIMAL_ITER, &nrp);
      checkMSKerror(err,"MSK_getintinf","getIterationsCount");
      nr = nr+nrp;
    }

  if( solver == MSK_OPTIMIZER_INTPNT )
    {
      err = MSK_getintinf(getMutableLpPtr(), MSK_IINF_INTPNT_ITER, &nr);
      checkMSKerror(err,"MSK_getintinf","getIterationsCount");
    }
  
  
  return nr;
}

//-----------------------------------------------------------------------------
// Returns one dual ray

std::vector<double*> OsiMskSolverInterface::getDualRays(int maxNumRays) const
{
  debugMessage("OsiMskSolverInterface::getDualRays(%d)\n", maxNumRays);

   OsiMskSolverInterface solver(*this);

   int numrows = getNumRows(), status, solution, r;
   if( probtypemip_ == false )
   {
	 if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	 else
	 {
		if( definedSolution( MSK_SOL_ITR ) == true )
			solution = MSK_SOL_ITR;
		else
			return std::vector<double*>();	
	 }
   }
   else
   {
	 if( definedSolution( MSK_SOL_ITG ) == true )
		solution = MSK_SOL_ITG;
	 else
        return std::vector<double*>();     
   }

   double *farkasray = new double[numrows];
   r = MSK_getsolution(
	      getMutableLpPtr(),
	      solution,
	      NULL, 
	      &status, 
	      NULL, 
	      NULL, 
	      NULL,
	      NULL,
	      NULL,
	      farkasray,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL
  );
  checkMSKerror( r, "MSK_getsolution", "getDualRays" );

  if( status != MSK_SOL_STA_DUAL_INFEAS_CER )
  { 
    delete[] farkasray;
    return std::vector<double*>();   
  }
  else
    return std::vector<double*>(1, farkasray);
}

//-----------------------------------------------------------------------------
// Returns one primal ray

std::vector<double*> OsiMskSolverInterface::getPrimalRays(int maxNumRays) const
{
  debugMessage("OsiMskSolverInterface::getPrimalRays(%d)\n", maxNumRays);

   OsiMskSolverInterface solver(*this);

   int numrows = getNumRows(), status, solution, r;
   if( probtypemip_ == false )
   {
	 if( definedSolution( MSK_SOL_BAS ) == true )
		solution = MSK_SOL_BAS;
	 else
	 {
		if( definedSolution( MSK_SOL_ITR ) == true )
			solution = MSK_SOL_ITR;
		else
			return std::vector<double*>();	
	 }
   }
   else
   {
	 if( definedSolution( MSK_SOL_ITG ) == true )
		solution = MSK_SOL_ITG;
	 else
        return std::vector<double*>();     
   }

   double *farkasray = new double[numrows];
   r = MSK_getsolution(
	      getMutableLpPtr(),
	      solution,
	      NULL, 
	      &status, 
	      NULL, 
	      NULL, 
	      NULL,
	      NULL,
	      farkasray,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL,
	      NULL
  );

  checkMSKerror( r, "MSK_getsolution", "getPrimalRays" );

  if( status != MSK_SOL_STA_PRIM_INFEAS_CER )
  { 
    delete[] farkasray;
    return std::vector<double*>();   
  }
  else
    return std::vector<double*>(1, farkasray);
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

//-----------------------------------------------------------------------------
// Sets a variabels objective coeff.

void OsiMskSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
  debugMessage("OsiMskSolverInterface::setObjCoeff(%d, %g)\n", elementIndex, elementValue);

  int err = MSK_putclist(
	  getLpPtr( OsiMskSolverInterface::FREECACHED_COLUMN ), 
	  1, 
	  &elementIndex, 
	  &elementValue
	);

    checkMSKerror(err, "MSK_putclist", "setObjCoeff");
}

//-----------------------------------------------------------------------------
// Sets a list of objective coeff.

void OsiMskSolverInterface::setObjCoeffSet(const int* indexFirst,
					   const int* indexLast,
					   const double* coeffList)
{
  debugMessage("OsiMskSolverInterface::setObjCoeffSet(%p, %p, %p)\n", indexFirst, indexLast, coeffList);

   const int cnt = indexLast - indexFirst;
   int err = MSK_putclist(
		getLpPtr(OsiMskSolverInterface::FREECACHED_COLUMN), cnt,
		const_cast<int*>(indexFirst),
		const_cast<double*>(coeffList)
	);

    checkMSKerror(err, "MSK_putclist", "setObjCoeffSet");
}

//-----------------------------------------------------------------------------
// Sets lower bound on one specific column

void OsiMskSolverInterface::setColLower(int elementIndex, double elementValue)
{
  debugMessage("OsiMskSolverInterface::setColLower(%d, %g)\n", elementIndex, elementValue);

  int finite = 1;
  if( elementValue <= -getInfinity() )
	  finite = 0;
  
  int err = MSK_chgbound(  
        getMutableLpPtr(), 
		0,
		elementIndex,
		1,
		finite,  
		elementValue 
		);

  checkMSKerror( err, "MSK_chgbound", "setColLower" );
   
  if( collower_ != NULL )
    collower_[elementIndex] = elementValue;
}

//-----------------------------------------------------------------------------
// Sets upper bound on one specific column. 

void OsiMskSolverInterface::setColUpper(int elementIndex, double elementValue)
{  
  debugMessage("OsiMskSolverInterface::setColUpper(%d, %g)\n", elementIndex, elementValue);

  int finite = 1;
  if( elementValue >= getInfinity() )
	  finite = 0;

  int err = MSK_chgbound( 
		getMutableLpPtr(), 
		0,
		elementIndex,
		0,
		finite,  
		elementValue 
    );

  checkMSKerror( err, "MSK_chgbound", "setColUpper" );
    
  if( colupper_ != NULL )
    colupper_[elementIndex] = elementValue;
} 

//-----------------------------------------------------------------------------
// Sets upper and lower bound on one specific column 

void OsiMskSolverInterface::setColBounds( int elementIndex, double lower, double upper )
{
  debugMessage("OsiMskSolverInterface::setColBounds(%d, %g, %g)\n", elementIndex, lower, upper);

  setColLower(elementIndex, lower);
  setColUpper(elementIndex, upper);
}

//-----------------------------------------------------------------------------
// Sets upper and lower bounds on a list of columns. Due to the strange storage of
// boundlist, it is not possible to change all the bounds in one call to MOSEK,
// so the standard method is used. 

void OsiMskSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
   debugMessage("OsiMskSolverInterface::setColSetBounds(%p, %p, %p)\n", indexFirst, indexLast, boundList);
   OsiSolverInterface::setColSetBounds( indexFirst, indexLast, boundList );
}

//-----------------------------------------------------------------------------
// Sets the lower bound on a row

void
OsiMskSolverInterface::setRowLower( int i, double elementValue )
{
  debugMessage("OsiMskSolverInterface::setRowLower(%d, %g)\n", i, elementValue);

  double rhs   = getRightHandSide()[i];
  double range = getRowRange()[i];
  char   sense = getRowSense()[i];
  double lower, upper;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( lower != elementValue ) 
  {
      convertBoundToSense( elementValue, upper, sense, rhs, range );
      setRowType( i, sense, rhs, range );
  }
}

//-----------------------------------------------------------------------------
// Sets the upper bound on a row 

void
OsiMskSolverInterface::setRowUpper( int i, double elementValue )
{
  debugMessage("OsiMskSolverInterface::setRowUpper(%d, %g)\n", i, elementValue);

  double rhs   = getRightHandSide()[i];
  double range = getRowRange()[i];
  char   sense = getRowSense()[i];
  double lower, upper;

  convertSenseToBound( sense, rhs, range, lower, upper );
  if( upper != elementValue ) 
  {
      convertBoundToSense( lower, elementValue, sense, rhs, range );
      setRowType( i, sense, rhs, range );
  }
}

//-----------------------------------------------------------------------------
// Sets the upper and lower bound on a row 

void
OsiMskSolverInterface::setRowBounds( int elementIndex, double lower, double upper )
{
  debugMessage("OsiMskSolverInterface::setRowBounds(%d, %g, %g)\n", elementIndex, lower, upper);

  double rhs, range;
  char sense;
  
  convertBoundToSense( lower, upper, sense, rhs, range );
  setRowType( elementIndex, sense, rhs, range );
}

//-----------------------------------------------------------------------------
// Sets the triplet on a row  

void
OsiMskSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  debugMessage("OsiMskSolverInterface::setRowType(%d, %c, %g, %g)\n", i, sense, rightHandSide, range);
  double rub,rlb;
  int rtag;

  MskConvertSenseToBound(sense, range, rightHandSide, rlb, rub, rtag); 

  int err = MSK_putbound(
                  getMutableLpPtr(),
                  1, 
                  i, 
                  rtag, 
                  rlb, 
                  rub);

  if( rowsense_ != NULL )
     rowsense_[i] = sense;

  if( rowrange_ != NULL )
     rowrange_[i] = range;

  if( rhs_      != NULL )
     rhs_[i]      = rightHandSide;
                  
  checkMSKerror( err, "MSK_putbound", "setRowType" );
}

//-----------------------------------------------------------------------------
// Set upper og lower bounds for a lisit of rows. Due to the strange storage of
// boundlist, it is not possible to change all the bounds in one call to MOSEK,
// so the standard method is used. 

void OsiMskSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  debugMessage("OsiMskSolverInterface::setRowSetBounds(%p, %p, %p)\n", indexFirst, indexLast, boundList);

   const int cnt = indexLast - indexFirst;
   if (cnt <= 0)
      return;

   for (int i = 0; i < cnt; ++i)
      setRowBounds(indexFirst[i], boundList[2*i], boundList[2*i+1]);
}


//-----------------------------------------------------------------------------
// Set triplets for a list of rows 

void
OsiMskSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  debugMessage("OsiMskSolverInterface::setRowSetTypes(%p, %p, %p, %p, %p)\n", 
     indexFirst, indexLast, senseList, rhsList, rangeList);
  
   const int cnt = indexLast - indexFirst;
   if (cnt <= 0)
      return;

   for (int i = 0; i < cnt; ++i)
      setRowType(indexFirst[i], senseList[i], rhsList[i], rangeList[i]);	
}


//-----------------------------------------------------------------------------
// Sets a variabel to continuous

void
OsiMskSolverInterface::setContinuous(int index)
{
  debugMessage("OsiMskSolverInterface::setContinuous(%d)\n", index);

  assert(coltype_ != NULL);
  assert(coltypesize_ >= getNumCols());

  coltype_[index] = 'C';
  
  int err = MSK_putvartype( getMutableLpPtr(), index, MSK_VAR_TYPE_CONT);
  checkMSKerror( err, "MSK_putvartype", "setContinuous" );    
}

//-----------------------------------------------------------------------------
// Sets a variabel to integer

void
OsiMskSolverInterface::setInteger(int index)
{
  debugMessage("OsiMskSolverInterface::setInteger(%d)\n", index);

  assert(coltype_ != NULL);
  assert(coltypesize_ >= getNumCols());
  
  coltype_[index] = 'I';
  
  int err = MSK_putvartype( getMutableLpPtr(), index, MSK_VAR_TYPE_INT);
  
  checkMSKerror( err, "MSK_putvartype", "setInteger" );
}

//-----------------------------------------------------------------------------
// Sets a list of variables to continuous

void
OsiMskSolverInterface::setContinuous(const int* indices, int len)
{
  debugMessage("OsiMskSolverInterface::setContinuous(%p, %d)\n", indices, len);

  for( int i = 0; i < len; ++i )
     setContinuous(indices[i]);
}

//-----------------------------------------------------------------------------
// Sets a list of variables to integer 

void
OsiMskSolverInterface::setInteger(const int* indices, int len)
{
  debugMessage("OsiMskSolverInterface::setInteger(%p, %d)\n", indices, len);

  for( int i = 0; i < len; ++i )
     setInteger(indices[i]);
}

//-----------------------------------------------------------------------------
// Sets the direction of optimization

void OsiMskSolverInterface::setObjSense(double s) 
{
	debugMessage("OsiMskSolverInterface::setObjSense(%g)\n", s);

	int err;
	if( s == +1.0 )
		err = MSK_putintparam(
                  getMutableLpPtr(),
                  MSK_IPAR_OBJECTIVE_SENSE,
                  MSK_OBJECTIVE_SENSE_MIN
      );
	else
		err = MSK_putintparam(
                  getMutableLpPtr(),
                  MSK_IPAR_OBJECTIVE_SENSE,
                  MSK_OBJECTIVE_SENSE_MAX
      );
      
	checkMSKerror(err,"MSK_putintparam","setObjSense");
}

//-----------------------------------------------------------------------------
// Sets the col solution. This is not used in MOSEK. Displays a warning.

void OsiMskSolverInterface::setColSolution(const double * cs) 
{
  debugMessage("OsiMskSolverInterface::setColSolution(%p)\n", cs);
 
  int nc = getNumCols();
  if( cs == NULL )
  {
    freeCachedResults();
  }
  else if( nc > 0 )
  {
      if( colsol_ != NULL )
          delete[] colsol_;
      
	  colsol_ = new double[nc];
      
      CoinDisjointCopyN( cs, nc, colsol_ );
      
      #if defined MSK_WARNING_ON
      OsiMSK_warning("OsiMskSolverInterface::setColSolution()",
      "Code has no effect in MOSEK");
      #endif
  }
}

//-----------------------------------------------------------------------------
// Sets the rowprices. This is not used in MOSEK. Displays a warning.

void OsiMskSolverInterface::setRowPrice(const double * rs) 
{
  debugMessage("OsiMskSolverInterface::setRowPrice(%p)\n", rs);

  int nr = getNumRows();

  if( rs == NULL )
    freeCachedResults();
  else if( nr > 0 )
  {
      if ( rowsol_ == NULL )
	    rowsol_ = new double[nr];
      
      CoinDisjointCopyN( rs, nr, rowsol_ ); 
      
     #if defined MSK_WARNING_ON
     OsiMSK_warning("OsiMskSolverInterface::setRowPrice()",
     "Code has no effect in MOSEK");
     #endif
   }
}


//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################

//-----------------------------------------------------------------------------
// Adds a column to the MOSEK task

void 
OsiMskSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  debugMessage("OsiMskSolverInterface::addCol(%p, %g, %g, %g)\n", &vec, collb, colub, obj);

  int nc = getNumCols();
  assert(coltypesize_ >= nc);

  resizeColType(nc + 1);
  coltype_[nc] = 'C';

  int tag,start = 0, ends = vec.getNumElements();

  double inf = getInfinity();
  if(collb > -inf && colub >= inf)
	tag = MSK_BK_LO;
  else if(collb <= -inf && colub < inf)
	tag = MSK_BK_UP;
  else if(collb > -inf && colub < inf)
	tag = MSK_BK_RA;
  else if(collb <= -inf && colub >= inf)
    tag = MSK_BK_FR;
  else
    throw CoinError("Bound error", "addCol", "OsiMSKSolverInterface");

  int err = MSK_appendvars(
	  getLpPtr(),
	  1,
	  const_cast<double*> (&obj),
	  &start,
	  &ends,
	  const_cast<int*>(vec.getIndices()),
	  const_cast<double*>(vec.getElements()),
	  &tag,
	  const_cast<double*> (&collb),
	  const_cast<double*> (&colub)
	);

  checkMSKerror( err, "MSK_appendvars", "addCol" );
}

//-----------------------------------------------------------------------------
// Adds a list of columns to the MOSEK task

void 
OsiMskSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  debugMessage("OsiMskSolverInterface::addCols(%d, %p, %p, %p, %p)\n", numcols, cols, collb, colub, obj);

  int i, nz = 0, err = MSK_RES_OK;
  
  // For efficiency we put hints on the total future size
  err = MSK_getmaxnumanz(
                     getLpPtr(),
                     &nz);
                     
  checkMSKerror( err, "MSK_getmaxanz", "addCols" );

  for( i = 0; i < numcols; ++i)
    nz += cols[i]->getNumElements();
  
  err = MSK_putmaxnumanz(
                     getLpPtr(),
                     nz);
                     
  checkMSKerror( err, "MSK_putmaxanz", "addCols" );
          
  err = MSK_putmaxnumvar(
                     getLpPtr(),
                     numcols+getNumCols());
                     
  checkMSKerror( err, "MSK_putmaxnumvar", "addCols" );

  for( i = 0; i < numcols; ++i )
   addCol( *(cols[i]), collb[i], colub[i], obj[i] );
}

//-----------------------------------------------------------------------------
// Deletes a list of columns from the MOSEK task 

void 
OsiMskSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  debugMessage("OsiMskSolverInterface::deleteCols(%d, %p)\n", num, columnIndices);
  int err;
  err = MSK_remove(
           getLpPtr( OsiMskSolverInterface::KEEPCACHED_ROW ),
           0,
           num,
           const_cast<int*>(columnIndices)
    );

  checkMSKerror( err, "MSK_remove", "deleteCols" );
}

//-----------------------------------------------------------------------------
// Adds a row in bound form to the MOSEK task

void 
OsiMskSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  debugMessage("OsiMskSolverInterface::addRow(%p, %g, %g)\n", &vec, rowlb, rowub);

  getNumRows();

  int tag,start = 0, ends = vec.getNumElements();
  double inf = getInfinity();
  if(rowlb > -inf && rowub >= inf)
	tag = MSK_BK_LO;
  else if(rowlb <= -inf && rowub < inf)
	tag = MSK_BK_UP;
  else if(rowlb > -inf && rowub < inf)
	tag = MSK_BK_RA;
  else if(rowlb <= -inf && rowub >= inf)
    tag = MSK_BK_FR;
  else
    throw CoinError("Bound error", "addRow", "OsiMSKSolverInterface");

  int err = MSK_appendcons(
	  getLpPtr( OsiMskSolverInterface::KEEPCACHED_COLUMN ),
	  1,
	  &start,
	  &ends,
	  const_cast<int*>(vec.getIndices()),
	  const_cast<double*>(vec.getElements()),
	  &tag,
	  const_cast<double*> (&rowlb),
	  const_cast<double*> (&rowub)
	);

  checkMSKerror( err, "MSK_appendvars", "addCol" );
}

//-----------------------------------------------------------------------------
// Adds a row in triplet form to the MOSEK task

void 
OsiMskSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  debugMessage("OsiMskSolverInterface::addRow(%p, %c, %g, %g)\n", &vec, rowsen, rowrhs, rowrng);

  double lb,ub;
  convertSenseToBound( rowsen, rowrhs, rowrng, lb, ub );
  addRow(vec, lb, ub);
}

//-----------------------------------------------------------------------------
// Adds a serie of rows in bound form to the MOSEK task

void 
OsiMskSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  debugMessage("OsiMskSolverInterface::addRows(%d, %p, %p, %p)\n", numrows, rows, rowlb, rowub);

  int i,nz = 0, err = MSK_RES_OK;
  
  // For efficiency we put hints on the total future size
  err = MSK_getmaxnumanz(
                     getLpPtr(),
                     &nz);
                     
  checkMSKerror( err, "MSK_getmaxanz", "addRows" );
  
  
  for( i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();
  
  err = MSK_putmaxnumanz(
                     getLpPtr(),
                     nz);
                     
  checkMSKerror( err, "MSK_putmaxanz", "addRows" );
          
  err = MSK_putmaxnumcon(
                     getLpPtr(),
                     numrows+getNumRows());
                    
  checkMSKerror( err, "MSK_putmaxnumcon", "addRows" );

  for( i = 0; i < numrows; ++i )
    addRow( *(rows[i]), rowlb[i], rowub[i] );
}

//-----------------------------------------------------------------------------
// Adds a list of rows in triplet form to the MOSEK task

void 
OsiMskSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  debugMessage("OsiMskSolverInterface::addRows(%d, %p, %p, %p, %p)\n", numrows, rows, rowsen, rowrhs, rowrng);

  int i, err = MSK_RES_OK, nz = 0;
  
    // For efficiency we put hints on the total future size
  for( i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();
  
  err = MSK_putmaxnumanz(
                     getLpPtr(),
                     nz);
                     
  checkMSKerror( err, "MSK_putmaxanz", "addRows" );
          
  err = MSK_putmaxnumcon(
                     getLpPtr(),
                     numrows);
                     
  checkMSKerror( err, "MSK_putmaxnumcon", "addRows" );

  for( i = 0; i < numrows; ++i )
    addRow( *(rows[i]), rowsen[i], rowrhs[i], rowrng[i] );
}

//-----------------------------------------------------------------------------
// Deletes a list of rows the MOSEK task 

void 
OsiMskSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  debugMessage("OsiMskSolverInterface::deleteRows(%d, %p)\n", num, rowIndices);
  int err;
  err = MSK_remove(
           getLpPtr( OsiMskSolverInterface::KEEPCACHED_COLUMN ),
           1,
           num,
           const_cast<int*>(rowIndices)
           );

  checkMSKerror( err, "MSK_remove", "deleteRows" );
}

//#############################################################################
// Methods to input a problem
//#############################################################################

//-----------------------------------------------------------------------------
// Loads a problem. Should have its "own" implementation so we don't have to convert
// to triplet, since this is convertet back in the load function called. But
// for simplicity, this is not done. 

void
OsiMskSolverInterface::loadProblem( 
                          const CoinPackedMatrix& matrix,
				          const double* collb, 
                          const double* colub,
				          const double* obj,
				          const double* rowlb, 
                          const double* rowub )
{
  debugMessage("OsiMskSolverInterface::loadProblem(%p, %p, %p, %p, %p, %p)\n", &matrix, collb, colub, obj, rowlb, rowub);

  const double inf = getInfinity();
  
  int nrows = matrix.getNumRows();
  char   * rowSense = new char  [nrows];
  double * rowRhs   = new double[nrows];
  double * rowRange = new double[nrows];
  
  int i;
  if( rowlb == NULL && rowub == NULL)
      for ( i = nrows - 1; i >= 0; --i )
          convertBoundToSense( -inf, inf, rowSense[i], rowRhs[i], rowRange[i] );
  else if( rowlb == NULL)
      for ( i = nrows - 1; i >= 0; --i )
          convertBoundToSense( -inf, rowub[i], rowSense[i], rowRhs[i], rowRange[i] );
  else if( rowub == NULL)
      for ( i = nrows - 1; i >= 0; --i )
          convertBoundToSense( rowlb[i], inf, rowSense[i], rowRhs[i], rowRange[i] );
  else
      for ( i = nrows - 1; i >= 0; --i )
          convertBoundToSense( rowlb[i], rowub[i], rowSense[i], rowRhs[i], rowRange[i] );

  loadProblem( matrix, collb, colub, obj, rowSense, rowRhs, rowRange );
  delete [] rowSense;
  delete [] rowRhs;
  delete [] rowRange;
}

//-----------------------------------------------------------------------------
// Loads a problem 

void
OsiMskSolverInterface::assignProblem( 
                      CoinPackedMatrix*& matrix,
				      double*& collb, 
                      double*& colub,
				      double*& obj,
				      double*& rowlb, 
                      double*& rowub )
{
  debugMessage("OsiMskSolverInterface::assignProblem()\n");

  loadProblem( *matrix, collb, colub, obj, rowlb, rowub );
  delete matrix;   matrix = 0;
  delete[] collb;  collb  = 0;
  delete[] colub;  colub  = 0;
  delete[] obj;    obj    = 0;
  delete[] rowlb;  rowlb  = 0;
  delete[] rowub;  rowub  = 0;
}

//-----------------------------------------------------------------------------
// Loads a problem 

void
OsiMskSolverInterface::loadProblem( 
                    const CoinPackedMatrix& matrix,
				    const double* collb, 
                    const double* colub,
				    const double* obj,
				    const char* rowsen, 
                    const double* rowrhs,
				    const double* rowrng )
{
  debugMessage("OsiMskSolverInterface::loadProblem(%p, %p, %p, %p, %p, %p, %p)\n",
     &matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
     
  int nc=matrix.getNumCols();
  int nr=matrix.getNumRows();
  
  if( nr == 0 || nc == 0 )   
    gutsOfDestructor();
  else
  {    
	assert( rowsen != NULL );
    assert( rowrhs != NULL );
    int i,j;
      
    double    * ob;
	int       * rtag  = NULL;
    double    * rlb   = NULL;
    double    * rub   = NULL;
	int       * ctag  = NULL;
	int       * cends = NULL;
	const int *len;
	const int *start;
    double    * clb   = NULL;
    double    * cub   = NULL;

    if( obj != NULL )
	  ob=const_cast<double*>(obj);
    else
	{
      ob = new double[nc];
	  CoinFillN(ob, nc, 0.0);
	}

  rtag = new int[nr];
  rlb  = new double[nr];
  rub  = new double[nr];
  
  for( i=0; i < nr; i++ )
	  MskConvertSenseToBound( rowsen[i], rowrng[i], rowrhs[i], rlb[i], rub[i], rtag[i]);

  bool freeMatrixRequired = false;
  CoinPackedMatrix * m = NULL;
  if( !matrix.isColOrdered() ) 
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
	
  double inf =getInfinity();
  ctag       = new int[nc];
  cends      = new int[nc];
  len        = (m->getVectorLengths());
  start      = (m->getVectorStarts());
  clb        = new double[nc];
  cub        = new double[nc];

  if( collb == NULL && colub == NULL )	   
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j] + len[j];
         MskConvertColBoundToTag(0, inf, clb[j], cub[j], ctag[j]);
      }
  else if( collb == NULL )
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j] + len[j];
         MskConvertColBoundToTag( 0, colub[j], clb[j], cub[j], ctag[j]);
      }
  else if( colub == NULL )
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j] + len[j];
         MskConvertColBoundToTag(collb[j], inf, clb[j], cub[j], ctag[j]);
      }
  else
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j] + len[j];
         MskConvertColBoundToTag( collb[j], colub[j], clb[j], cub[j], ctag[j]);
      }

  int err=MSK_inputdata(
	getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ),
	nr,
	nc,
	nr,
	nc,
	ob, 
	0.0,
	const_cast<int *>(m->getVectorStarts()),
	cends,
	const_cast<int *>(m->getIndices()),
	const_cast<double *>(m->getElements()), 
	rtag, 
	rlb, 
	rub, 
	ctag, 
	clb, 
	cub
  ); 

  checkMSKerror( err, "MSK_inputdata", "loadProblem" );
            
  if( obj   == NULL )
    delete[] ob;
  
  delete[] rtag;
  delete[] rlb;
  delete[] rub;
  delete[] ctag;
  delete[] clb;
  delete[] cub;
  delete[] cends;
  
  if ( freeMatrixRequired ) 
	delete m;

  resizeColType(nc);
  CoinFillN(coltype_, nc, 'C');
}

}
   
//-----------------------------------------------------------------------------
// Assigns a problem 

void
OsiMskSolverInterface::assignProblem( CoinPackedMatrix*& matrix,
				      double*& collb, double*& colub,
				      double*& obj,
				      char*& rowsen, double*& rowrhs,
				      double*& rowrng )
{
  debugMessage("OsiMskSolverInterface::assignProblem()\n");

   loadProblem( *matrix, collb, colub, obj, rowsen, rowrhs, rowrng );
   delete matrix;   matrix = 0;
   delete[] collb;  collb = 0;
   delete[] colub;  colub = 0;
   delete[] obj;    obj = 0;
   delete[] rowsen; rowsen = 0;
   delete[] rowrhs; rowrhs = 0;
   delete[] rowrng; rowrng = 0;
}

//-----------------------------------------------------------------------------
// Loads a problem 

void
OsiMskSolverInterface::loadProblem(
                   const int numcols, 
                   const int numrows,
				   const int* start, 
                   const int* index,
				   const double* value,
				   const double* collb, 
                   const double* colub,   
				   const double* obj,
				   const double* rowlb, 
                   const double* rowub )
{
  debugMessage("OsiMskSolverInterface::loadProblem()\n");

  const double inf = getInfinity();
  
  char   * rowSense = new char  [numrows];
  double * rowRhs   = new double[numrows];
  double * rowRange = new double[numrows];
  
  for ( int i = numrows - 1; i >= 0; --i ) 
  {
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
// Loads a problem 

void
OsiMskSolverInterface::loadProblem(
                   const int numcols, 
                   const int numrows,
				   const int* start, 
                   const int* index,
				   const double* value,
				   const double* collb, 
                   const double* colub,   
				   const double* obj,
				   const char* rowsen, 
                   const double* rowrhs,
				   const double* rowrng )
{
  debugMessage("OsiMskSolverInterface::loadProblem(%d, %d, %p, %p, %p, %p, %p, %p, %p, %p, %p)\n",
     numcols, numrows, start, index, value, collb, colub, obj, rowsen, rowrhs, rowrng);

  const int nc = numcols;
  const int nr = numrows;

  if( nr == 0 || nc == 0 )  
    gutsOfDestructor();     
  else
  {
    assert( rowsen != NULL );
    assert( rowrhs != NULL );

    int i,j;
      
    double   * ob;
	int      * rtag =NULL;
    double   * rlb = NULL;
    double   * rub = NULL;
	int      * ctag =NULL;
	int      * cends = NULL;
    double   * clb = NULL;
    double   * cub = NULL;

  if( obj != NULL )
	ob=const_cast<double*>(obj);
  else
  {
     ob = new double[nc];
	 CoinFillN(ob, nc, 0.0);
  }

  rtag = new int[nr];
  rlb  = new double[nr];
  rub  = new double[nr];

  for( i=0; i < nr; i++ )
	  MskConvertSenseToBound( rowsen[i], rowrng[i], rowrhs[i], rlb[i], rub[i], rtag[i]);
	
  double inf = getInfinity();
  ctag       = new int[nc];
  cends      = new int[nc];
  clb        = new double[nc];
  cub        = new double[nc];

  if( collb == NULL && colub == NULL )	   
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j+1];
         MskConvertColBoundToTag(0, inf, clb[j], cub[j], ctag[j]);
      }
  else if( collb == NULL )
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j+1];
         MskConvertColBoundToTag( 0, colub[j], clb[j], cub[j], ctag[j]);
      }
  else if( colub == NULL )
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j+1];
         MskConvertColBoundToTag(collb[j], inf, clb[j], cub[j], ctag[j]);
      }
  else
      for(j=0; j < nc; j++)
      {
         cends[j] = start[j+1];
         MskConvertColBoundToTag( collb[j], colub[j], clb[j], cub[j], ctag[j]);
      }

  int err=MSK_inputdata(
    getLpPtr( OsiMskSolverInterface::FREECACHED_RESULTS ),
	nr,
	nc,
	nr,
	nc,
	ob, 
	0.0,
	const_cast<int *>(start),
	cends,
	const_cast<int *>(index),
	const_cast<double *>(value), 
	rtag, 
	rlb, 
	rub, 
	ctag, 
	clb, 
	cub
  );

  checkMSKerror( err, "MSK_inputdata", "loadProblem3" );
            
  if( obj   == NULL )
    delete[] ob;
  
  delete[] rtag;
  delete[] rlb;
  delete[] rub;
  delete[] ctag;
  delete[] clb;
  delete[] cub;

  resizeColType(nc);
  CoinFillN(coltype_, nc, 'C');
  }
}
 
//-----------------------------------------------------------------------------
// Reads a MPS file with Coin native MPS reader. If marked code is switch on
// then MOSEK file reader is used, and .gz files can be read aswell.

int OsiMskSolverInterface::readMps( const char * filename,
				     const char * extension )
{
  debugMessage("OsiMskSolverInterface::readMps(%s, %s)\n", filename, extension);
  #if 0
  std::string f(filename);
  std::string e(extension);
  std::string fullname = f + "." + e;
  int err=MSK_readdata(getMutableLpPtr(),const_cast<char*>(fullname.c_str()));
  checkMSKerror(err, "MSK_readdatafile", "readMps" );
  return err;
  #else
  return OsiSolverInterface::readMps(filename,extension);
  #endif
}

//-----------------------------------------------------------------------------
// Writes the problem in MPS format, uses MOSEK writer. Notice that generic names is
// switched on, this could be avoided by using MOSEK reader instead of Coin native reader.

void OsiMskSolverInterface::writeMps( const char * filename,
				      const char * extension,
				      double objSense ) const
{
  debugMessage("OsiMskSolverInterface::writeMps(%s, %s, %g)\n", filename, extension, objSense);
  std::string f(filename);
  std::string e(extension);  
  std::string fullname = f + "." + e;
  int err = MSK_writedata( getMutableLpPtr(), const_cast<char*>( fullname.c_str() ));
  checkMSKerror( err, "MSK_writedatafile", "writeMps" );
}

//#############################################################################
// MSK specific public interfaces
//#############################################################################

//-----------------------------------------------------------------------------
// Returns MOSEK task in the interface object

MSKenv_t OsiMskSolverInterface::getEnvironmentPtr()
{
  assert( env_ != NULL );
  return env_;
}

//-----------------------------------------------------------------------------
// Returns MOSEK task in the interface object

MSKtask_t OsiMskSolverInterface::getLpPtr( int keepCached )
{
  freeCachedData( keepCached );
  return getMutableLpPtr();
}

//-----------------------------------------------------------------------------
// Returns the coltype_ array

const char * OsiMskSolverInterface::getCtype() const
{
  debugMessage("OsiMskSolverInterface::getCtype()\n");
  return coltype_;
}

//#############################################################################
// Static instance counter methods
//#############################################################################

//-----------------------------------------------------------------------------
// Increment the instance count, so we know when to close and open MOSEK.

void OsiMskSolverInterface::incrementInstanceCounter()
{
    debugMessage("OsiMskSolverInterface::incrementInstanceCounter()\n");

    if ( numInstances_ == 0 )
    {
      int err=0;
      char file[] = "mosek.alloc";
      checkMSKerror( err, "MSK_openmosek", "incrementInstanceCounter" );
      err = MSK_makeenv(&env_,NULL, NULL,NULL,file);
      checkMSKerror( err, "MSK_makeenv", "incrementInstanceCounter" );
      err = MSK_linkfunctoenvstream(env_, MSK_STREAM_LOG, NULL, printlog); 
      checkMSKerror( err, "MSK_linkfunctoenvstream", "incrementInstanceCounter" );
      err = MSK_initenv(env_);
      checkMSKerror( err, "MSK_initenv", "incrementInstanceCounter" );
    }
    numInstances_++;
}

//-----------------------------------------------------------------------------
// Decrement the instance count, so we know when to close and open MOSEK.

void OsiMskSolverInterface::decrementInstanceCounter()
{
  debugMessage("OsiMskSolverInterface::decrementInstanceCounter()\n");

  assert( numInstances_ != 0 );
  numInstances_--;
  if ( numInstances_ == 0 )
  {
     int err = MSK_deleteenv(&env_);
     checkMSKerror( err, "MSK_deleteenv", "decrementInstanceCounter" );
     env_ = NULL;
     MSK_ = NULL;
  }
}

//-----------------------------------------------------------------------------
// Returns the number of OsiMskSolverInterface objects in play

unsigned int OsiMskSolverInterface::getNumInstances()
{
    debugMessage("OsiMskSolverInterface::getNumInstances()\n");
    return numInstances_;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-----------------------------------------------------------------------------
// Constructor

OsiMskSolverInterface::OsiMskSolverInterface()
  : OsiSolverInterface(),
    Mskerr(MSK_RES_OK),
    InitialSolver(INITIAL_SOLVE),
    task_(NULL),
    hotStartCStat_(NULL),
    hotStartCStatSize_(0),
    hotStartRStat_(NULL),
    hotStartRStatSize_(0),
    hotStartMaxIteration_(1000000), 
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
	debugMessage("OsiMskSolverInterface::OsiMskSolverinterface()\n");
	incrementInstanceCounter();
	gutsOfConstructor();
}



//-----------------------------------------------------------------------------
// Clone from one to another object

OsiSolverInterface * OsiMskSolverInterface::clone(bool copyData) const
{
  debugMessage("OsiMskSolverInterface::clone(%d)\n", copyData);
  return( new OsiMskSolverInterface( *this ) );
}


//-----------------------------------------------------------------------------
// Copy constructor

OsiMskSolverInterface::OsiMskSolverInterface( const OsiMskSolverInterface & source )
  : OsiSolverInterface(source),
    Mskerr(MSK_RES_OK),
    InitialSolver(INITIAL_SOLVE),
    task_(NULL),
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
  debugMessage("OsiMskSolverInterface::OsiMskSolverInterface(%p)\n", &source);

  incrementInstanceCounter();  
  gutsOfConstructor();
  gutsOfCopy( source );
}

//-----------------------------------------------------------------------------
// Destructor

OsiMskSolverInterface::~OsiMskSolverInterface()
{
  debugMessage("OsiMskSolverInterface::~OsiMskSolverInterface()\n");

  gutsOfDestructor();
  decrementInstanceCounter();
}

//-----------------------------------------------------------------------------
// Assign operator

OsiMskSolverInterface& OsiMskSolverInterface::operator=( const OsiMskSolverInterface& rhs )
{
  debugMessage("OsiMskSolverInterface::operator=(%p)\n", &rhs);

  if (this != &rhs)
  {    
      OsiSolverInterface::operator=( rhs );
      gutsOfDestructor();
      gutsOfConstructor();
      if ( rhs.task_ !=NULL )
	    gutsOfCopy( rhs );
  }

  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

//-----------------------------------------------------------------------------
// Apply col cut

void OsiMskSolverInterface::applyColCut( const OsiColCut & cc )
{
  debugMessage("OsiMskSolverInterface::applyColCut(%p)\n", &cc);

  double * MskColLB = new double[getNumCols()];
  double * MskColUB = new double[getNumCols()];
  const CoinPackedVector & lbs = cc.lbs();
  const CoinPackedVector & ubs = cc.ubs();
  int i;
  
  for( i = 0; i < getNumCols(); ++i )
  {
          MskColLB[i] = getColLower()[i];
          MskColUB[i] = getColUpper()[i];
  }
    
  for( i = 0; i < lbs.getNumElements(); ++i ) 
    if ( lbs.getElements()[i] > MskColLB[lbs.getIndices()[i]] )
      setColLower( lbs.getIndices()[i], lbs.getElements()[i] );
  for( i = 0; i < ubs.getNumElements(); ++i )
    if ( ubs.getElements()[i] < MskColUB[ubs.getIndices()[i]] )
      setColUpper( ubs.getIndices()[i], ubs.getElements()[i] );
  
  delete[] MskColLB;
  delete[] MskColUB;
}

//-----------------------------------------------------------------------------
// Apply row cut

void OsiMskSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  debugMessage("OsiMskSolverInterface::applyRowCut(%p)\n", &rowCut);
  const CoinPackedVector & row=rowCut.row();
  addRow(row ,  rowCut.lb(),rowCut.ub());
}

//#############################################################################
// Private methods (non-static and static) and static data
//#############################################################################   

unsigned int OsiMskSolverInterface::numInstances_ = 0;
MSKhand_t OsiMskSolverInterface::MSK_=NULL;
MSKenv_t OsiMskSolverInterface::env_=NULL;
 
//-----------------------------------------------------------------------------
// Returns MOSEK task in object

MSKtask_t OsiMskSolverInterface::getMutableLpPtr() const
{
	//std::cout << "Task " << task_ << std::endl;
	if ( task_ == NULL )
	  {
	    debugMessage("OsiMskSolverInterface::getMutableLpPtr()\n");
	    char file[] = "MOSEK.log";          
	    assert(env_ != NULL);

	    int err = MSK_makeemptytask(env_,&task_);
	    checkMSKerror(err, "MSK_makeemptytask","getMutableLpPtr");

	    err = MSK_linkfunctotaskstream(task_, MSK_STREAM_LOG, NULL, printlog); 
	    checkMSKerror( err, "MSK_linkfunctotaskstream", "getMutableLpPtr" );

	    err = MSK_linkfiletotaskstream(task_, MSK_STREAM_LOG, file, 0);
	    checkMSKerror( err, "MSK_linkfiletotaskstream", "getMutableLpPtr" );

	    err = MSK_putdouparam( task_, MSK_DPAR_BASIS_REL_TOL_S, 0 );
	    checkMSKerror(err,"MSK_putdouparam","getMutableLpPtr()");

	    err = MSK_putintparam(task_, MSK_IPAR_OPTIMIZER, SOLVER);
	    checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");

	    // 	  err = MSK_putintparam(task_, MSK_IPAR_WRITE_GENERIC_NAMES, MSK_ON);
	    //           checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");  

	    err = MSK_putintparam(task_, MSK_IPAR_PRESOLVE_USE, MSK_ON);
	    checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");  
	    
	    err = MSK_putintparam(task_, MSK_IPAR_WRITE_DATA_FORMAT, MSK_DATA_FORMAT_MPS);
	    checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");  
	    
	    err = MSK_putintparam(task_, MSK_IPAR_WRITE_GENERIC_NAMES, MSK_ON);
	    checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");  
	    
 	    err = MSK_putintparam(task_, MSK_IPAR_LOG, MSK_OFF);
 	    checkMSKerror(err,"MSK_putintparam","getMutableLpPtr()");  
	    
	    std::string pn;
	    getStrParam(OsiProbName,pn);
	    MSK_puttaskname( task_, const_cast<char*>(pn.c_str()) );
	    checkMSKerror(err,"MSK_puttaskname","getMutableLpPtr()");   
	  }
	//int err= MSK_writedata(task_, "tmp.mbt");
	//checkMSKerror(err,"MSK_puttaskname","getMutableLpPtr()");   
	
	return task_;
}

//-----------------------------------------------------------------------------
// Makes a copy

void OsiMskSolverInterface::gutsOfCopy( const OsiMskSolverInterface & source )
{
    debugMessage("OsiMskSolverInterface::gutsOfCopy()\n");

	// Set Objective Sense
    setObjSense(source.getObjSense());

    InitialSolver = source.InitialSolver;

    // Set Rim and constraints
    const double* obj             = source.getObjCoefficients();
    const double* rhs             = source.getRightHandSide();
    const char* sense             = source.getRowSense();
    const CoinPackedMatrix * cols = source.getMatrixByCol();
    const double* lb              = source.getColLower();
    const double* ub              = source.getColUpper();

    loadProblem(*cols,lb,ub,obj,sense,rhs,source.getRowRange());

    // Set MIP information
    resizeColType(source.coltypesize_);
    CoinDisjointCopyN( source.coltype_, source.coltypesize_, coltype_ );

    // Set Solution
    setRowPrice(source.getRowPrice());
    setColSolution(source.getColSolution());
}


//-----------------------------------------------------------------------------
// Empty function

void OsiMskSolverInterface::gutsOfConstructor()
{ 
  debugMessage("OsiMskSolverInterface::gutsOfConstructor()\n");
}


//-----------------------------------------------------------------------------
// Function called from destructor

void OsiMskSolverInterface::gutsOfDestructor()
{  
  if ( task_ != NULL )
  {
      int err = MSK_deletetask(&task_);
      checkMSKerror( err, "MSK_deletetask", "gutsOfDestructor" );
      task_ = NULL;
      freeAllMemory();
  }
    
  assert( task_==NULL );
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


//-----------------------------------------------------------------------------
// Free function

void OsiMskSolverInterface::freeCachedColRim()
{
  freeCacheDouble( obj_ );  
  freeCacheDouble( collower_ ); 
  freeCacheDouble( colupper_ ); 
  assert( obj_==NULL );
  assert( collower_==NULL );
  assert( colupper_==NULL );
}

//-----------------------------------------------------------------------------
// Free function 

void OsiMskSolverInterface::freeCachedRowRim()
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

//-----------------------------------------------------------------------------
// Free function

void OsiMskSolverInterface::freeCachedMatrix()
{
  freeCacheMatrix( matrixByRow_ );
  freeCacheMatrix( matrixByCol_ );
  assert( matrixByRow_==NULL ); 
  assert( matrixByCol_==NULL ); 
}

//-----------------------------------------------------------------------------
//  Free function

void OsiMskSolverInterface::freeCachedResults()
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

//-----------------------------------------------------------------------------
//  Free function

void OsiMskSolverInterface::freeCachedData( int keepCached )
{
  if( !(keepCached & OsiMskSolverInterface::KEEPCACHED_COLUMN) )
    freeCachedColRim();
  if( !(keepCached & OsiMskSolverInterface::KEEPCACHED_ROW) )
    freeCachedRowRim();
  if( !(keepCached & OsiMskSolverInterface::KEEPCACHED_MATRIX) )
    freeCachedMatrix();
  if( !(keepCached & OsiMskSolverInterface::KEEPCACHED_RESULTS) )
    freeCachedResults();
}

//-----------------------------------------------------------------------------
//  Free function

void OsiMskSolverInterface::freeAllMemory()
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
