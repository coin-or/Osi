// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiSimplexInterface_H
#define OsiSimplexInterface_H

#include <string>
#include <vector>

#include "CoinMessageHandler.hpp"
#include "CoinPackedVectorBase.hpp"

#include "OsiSolverParameters.hpp"

class CoinPackedMatrix;
class CoinWarmStart;

//#############################################################################

/** Simplex Interface Abstract Base Class

Abstract Base Class for describing an advanced interface to a simplex solver.
When switched on allows great control of simplex iterations.  Also allows
access to tableau.
*/

class OsiSimplexInterface  {
   friend void OsiSimplexInterfaceCommonUnitTest(
      const OsiSolverInterface* emptySi,
      const std::string & mpsDir);

public:
  ///@name OsiSimplexInterface methods 
  //@{
  /**Enables normal operation of subsequent functions.
     This method is supposed to ensure that all typical things (like
     reduced costs, etc.) are updated when individual pivots are executed
     and can be queried by other methods.  says whether will be
     doing primal or dual
  */
  virtual void enableSimplexInterface(bool doingPrimal) {};

  ///Undo whatever setting changes the above method had to make
  virtual void disableSimplexInterface() {};

  ///Returns true if a basis is available
  virtual bool basisIsAvailable() = 0;

  /** The following two methods may be replaced by the
     methods of OsiSolverInterface using OsiWarmStartBasis if:
     1. OsiWarmStartBasis resize operation is implemented
     more efficiently and
     2. It is ensured that effects on the solver are the same

     Returns a basis status of the structural/artificial variables 
     At present as warm start i.e 0 free, 1 basic, 2 upper, 3 lower

     NOTE  artificials are treated as +1 elements so for <= rhs
     artificial will be at lower bound if constraint is tight
  */
  virtual void getBasisStatus(int* cstat, int* rstat) = 0;

  /** Set the status of structural/artificial variables and
      factorize, update solution etc 

     NOTE  artificials are treated as +1 elements so for <= rhs
     artificial will be at lower bound if constraint is tight
  */
  virtual int setBasisStatus(const int* cstat, const int* rstat) = 0;

  /** Perform a pivot by substituting a colIn for colOut in the basis. 
     The status of the leaving variable is given in statOut. Where
     1 is to upper bound, -1 to lower bound
  */
  virtual int pivot(int colIn, int colOut, int outStatus) = 0;

  /** Obtain a result of the primal pivot 
      Outputs: colOut -- leaving column, outStatus -- its status,
      t -- step size, and, if dx!=NULL, *dx -- primal ray direction.
      Inputs: colIn -- entering column, sign -- direction of its change (+/-1).
      Both for colIn and colOut, artificial variables are index by
      the negative of the row index minus 1.
      Return code (for now): 0 -- leaving variable found, 
      -1 -- everything else?
      Clearly, more informative set of return values is required 
  */
  virtual int primalPivotResult(int colIn, int sign, 
				int& colOut, int& outStatus, 
				double& t, CoinPackedVector* dx) = 0;

  /** Obtain a result of the dual pivot (similar to the previous method)
      Differences: entering variable and a sign of its change are now
      the outputs, the leaving variable and its statuts -- the inputs
      If dx!=NULL, then *dx contains dual ray
      Return code: same
  */
  virtual int dualPivotResult(int& colIn, int& sign, 
			      int colOut, int outStatus, 
			      double& t, CoinPackedVector* dx) = 0;

  ///Get the reduced gradient for the cost vector c 
  virtual void getReducedGradient(double* columnReducedCosts, 
				  double * duals,
				  const double * c) = 0;

  /** Set a new objective and apply the old basis so that the
      reduced costs are properly updated */
  virtual void setObjectiveAndRefresh(double* c) = 0;

  ///Get a row of the tableau (slack part in slack if not NULL)
  virtual void getBInvARow(int row, double* z, double * slack=NULL) = 0;

  ///Get a row of the basis inverse
  virtual void getBInvRow(int row, double* z) = 0;

  ///Get a column of the tableau
  virtual void getBInvACol(int col, double* vec) = 0;

  ///Get a column of the basis inverse
  virtual void getBInvCol(int col, double* vec) = 0;

  /** Get basic indices (order of indices corresponds to the
      order of elements in a vector retured by getBInvACol() and
      getBInvCol()).
  */
  virtual void getBasics(int* index) = 0;
  //@}


  ///@name Constructors and destructors
  //@{
    /// Default Constructor
    OsiSimplexInterface(); 
    
    /// Copy constructor 
    OsiSimplexInterface(const OsiSimplexInterface & rhs);
  
    /// Assignment operator 
    OsiSimplexInterface & operator=(const OsiSimplexInterface& rhs);
  
    /// Destructor 
    virtual ~OsiSimplexInterface();
  //@}

  //---------------------------------------------------------------------------

private:
  ///@name Private member data 
  //@{
 //@}
};

//#############################################################################
/** A function that tests the methods in the OsiSimplexInterface class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
OsiSimplexInterfaceCommonUnitTest(
   const OsiSolverInterface* emptySi,
   const std::string & mpsDir);

#endif
