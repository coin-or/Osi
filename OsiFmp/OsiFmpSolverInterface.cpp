// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "CoinPackedMatrix.hpp"
#include "OsiFmpSolverInterface.hpp"
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
	 m = new CoinPackedMatrix();
	 m->setExtraGap(0.0);
	 m->setExtraMajor(0.0);
	 m->operator=(matrix);
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
void OsiFmpSolverInterface::initialSolve()
{
	int tctn=0;
  char loglev[5];
	char scmd[80];


  // Switch off printing if asked to
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = getHintParam(OsiDoReducePrint,takeHint,strength);
  assert(gotHint);
  
	// Set the appropriate FortMP log level
	if (strength!=OsiHintIgnore&&takeHint) 
	{
		itoa(messageHandler()->logLevel(),loglev,10);
		strcpy(scmd,"LOG DISPLAY LEVEL ");
		strcat(scmd,loglev);
		SPECMDC(scmd,&tctn);
	}

	
	char tmpname[128];
	strcpy(tmpname,modelname_.c_str());

	// Set direction of the optimisation...
	if(isMax_)
	{
		strcpy(scmd,"MAX");
		SPECMDC(scmd,&tctn);
	}
  else
  {
		strcpy(scmd,"MIN");
		SPECMDC(scmd,&tctn);
	}

	strcpy(scmd,"MIP OFF");
	SPECMDC(scmd,&tctn);

	//PAT removed it temporarily
	// strcpy(scmd,"PRESOLVE ON");
	// SPECMDC(scmd,&tctn);

	// if there is no basis, allocates and use crash start
	if(basis_==NULL)
	{
		basis_=new int[nc_+nr_+1];
		strcpy(scmd,"SIMPLEX START CRASH");
		SPECMDC(scmd,&tctn);
	}


	// Call the solver
	SUBLP2C(&nr_, &nc_, &nz_, &ns_,
					tmpname,"NOSPECS ",
          aij_, rowin_, colin_,
					upb_, lob_, rhs_, lhs_, 
					cost_, mitype_, sref_, 
					sfun_, sbeg_, send_,
          &objconstant_, &obj_, 
					sol_, dsl_, basis_, &stsl_, &tctn);


	// If presolve fails, try without
/*	if (tctn==550)
	{
		tctn=0;
		strcpy(scmd,"PRESOLVE OFF");
	  SPECMDC(scmd,&tctn);

			// Call the solver
	  SUBLP2C(&nr_, &nc_, &nz_, &ns_,
					tmpname,"NOSPECS ",
          aij_, rowin_, colin_,
					upb_, lob_, rhs_, lhs_, 
					cost_, mitype_, sref_, 
					sfun_, sbeg_, send_,
          &objconstant_, &obj_, 
					sol_, dsl_, basis_, &stsl_, &tctn);

		
	}
*/
	// If an error occoured, set stsl to -1, used as abandoned state
  if (tctn!=0)
		stsl_=-1;

}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::resolve()
{
	int tctn=0;
  char loglev[5];
	char scmd[80];

  // Switch off printing if asked to
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = getHintParam(OsiDoReducePrint,takeHint,strength);
  assert(gotHint);
  
	if (strength!=OsiHintIgnore&&takeHint) 
	{
		itoa(messageHandler()->logLevel(),loglev,10);
		strcpy(scmd,"LOG DISPLAY LEVEL ");
		strcat(scmd,loglev);
		SPECMDC(scmd,&tctn);
	}

	char tmpname[128];
	strcpy(tmpname,modelname_.c_str());

	// Set direction of the optimisation...
	if(isMax_)
	{
		strcpy(scmd,"MAX");
		SPECMDC(scmd,&tctn);
	}
  else
  {
		strcpy(scmd,"MIN");
		SPECMDC(scmd,&tctn);
	}

	strcpy(scmd,"MIP OFF");
	SPECMDC(scmd,&tctn);


	if(basis_==NULL)
		basis_=new int[nc_+nr_+1];
	else
	{
		strcpy(scmd,"SIMPLEX START INPUT BASIS");
		SPECMDC(scmd,&tctn);
	}

	SUBLP2C(&nr_, &nc_, &nz_, &ns_,
					tmpname,"NOSPECS ",
          aij_, rowin_, colin_,
					upb_, lob_, rhs_, lhs_, 
					cost_, mitype_, sref_, 
					sfun_, sbeg_, send_,
          &objconstant_, cost_, 
					sol_, dsl_ , basis_, &stsl_, &tctn);

	// If an error occoured, set stsl to -1, used as abandoned state
  if (tctn!=0)
		stsl_=-1;

}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::branchAndBound()
{

	int tctn=0;
  char loglev[5];
	char scmd[80];

  // Switch off printing if asked to
  bool takeHint;
  OsiHintStrength strength;
  bool gotHint = getHintParam(OsiDoReducePrint,takeHint,strength);
  assert(gotHint);
  
	if (strength!=OsiHintIgnore&&takeHint) 
	{
		itoa(messageHandler()->logLevel(),loglev,10);
		strcpy(scmd,"LOG DISPLAY LEVEL ");
		strcat(scmd,loglev);
		SPECMDC(scmd,&tctn);
	}

	char tmpname[128];
	strcpy(tmpname,modelname_.c_str());

	// Set direction of the optimisation...
	if(isMax_)
	{
		strcpy(scmd,"MAX");
		SPECMDC(scmd,&tctn);
	}
  else
  {
		strcpy(scmd,"MIN");
		SPECMDC(scmd,&tctn);
	}

	// Enable branch and bound
	strcpy(scmd,"MIP ON");
	SPECMDC(scmd,&tctn);


	if(basis_==NULL)
		basis_=new int[nc_+nr_+1];
	else
	{
		strcpy(scmd,"SIMPLEX START INPUT BASIS");
		SPECMDC(scmd,&tctn);
	}

	SUBLP2C(&nr_, &nc_, &nz_, &ns_,
					tmpname,"NOSPECS ",
          aij_, rowin_, colin_,
					upb_, lob_, rhs_, lhs_, 
					cost_, mitype_, sref_, 
					sfun_, sbeg_, send_,
          &objconstant_, cost_, 
					sol_, dsl_ , basis_, &stsl_, &tctn);

	// If an error occoured, set stsl to -1, used as abandoned state
  if (tctn!=0)
		stsl_=-1;

}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiFmpSolverInterface::setIntParam(OsiIntParam key, int value)
{

  int tctn;
	char scmd[80];
	char val[20];

	// transforms the argument into a string
	itoa(value,val,10);

  switch (key) {
  case OsiMaxNumIteration:
    if (value < 0)
      return false;
		// send SPEC command to FortMP
		strcpy(scmd,"MAXIMUM SIMPLEX ITERATIONS = ");
		strcat(scmd,val);
		SPECMDC(scmd,&tctn);
		// also store the value for retrival
		OsiSolverInterface::setIntParam(key, value);
    return tctn == 0;
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
OsiFmpSolverInterface::setDblParam(OsiDblParam key, double value)
{
  int tctn=0;
	char scmd[80];
	char val[20];

	// transform the argument into a string
	gcvt(value,15,val);
  
  switch (key) {
  case OsiDualObjectiveLimit:
    return OsiSolverInterface::setDblParam(key, value);

  case OsiPrimalObjectiveLimit:
    return OsiSolverInterface::setDblParam(key, value);

  case OsiDualTolerance: 
		// send spec command to FortMP
		strcpy(scmd,"DUAL PIVOT ZERO TOLERANCE = ");
		strcat(scmd,val);
		SPECMDC(scmd,&tctn);
		// only store the value if successful
		if(tctn==0)
			OsiSolverInterface::setDblParam(key, value);

    return tctn == 0;

  case OsiPrimalTolerance:

		// Primal tolerance fixed in FortMP
		if(value==1e-35)
			return true;
		else
			return false;

  case OsiObjOffset: 
    objconstant_=value;
    return true;

  case OsiLastDblParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiFmpSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{

  switch (key) {
  case OsiProbName:

		modelname_=value.c_str();
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
OsiFmpSolverInterface::getIntParam(OsiIntParam key, int& value) const 
{
  switch (key) {
  case OsiMaxNumIteration:
    ;//PAT value = ekk_getImaxiter(getMutableModelPtr());
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
OsiFmpSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  // *THINK* : how to set opt limits for OSL?
  switch (key) {
  case OsiDualObjectiveLimit:
    OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiPrimalObjectiveLimit:
    OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiDualTolerance:
    OsiSolverInterface::getDblParam(key, value);
    ;//PAT value = ekk_getRtoldinf(getMutableModelPtr());
    break;
  case OsiPrimalTolerance:
		value=1e-35;
//    OsiSolverInterface::getDblParam(key, value);
    ;//PAT value = ekk_getRtolpinf(getMutableModelPtr());
    break;
  case OsiObjOffset:
    value=objconstant_;
    break;
  case OsiLastDblParam:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiFmpSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  switch (key) {
  case OsiProbName:
    ;//PAT value = ekk_getCname(getMutableModelPtr());
    break;
  case OsiSolverName:
    value = "FortMP";
    break;
  case OsiLastStrParam:
    return false;
  }
  return true;
}


//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiFmpSolverInterface::isAbandoned() const
{
  ;//PAT EKKModel* model = getMutableModelPtr();
  ;//PAT return (ekk_getIprobstat(model)==6);
	return true;
}

bool OsiFmpSolverInterface::isProvenOptimal() const
{
  return (stsl_ == 3 || stsl_==5);
}

bool OsiFmpSolverInterface::isProvenPrimalInfeasible() const
{
  return (stsl_ == 1);
}

bool OsiFmpSolverInterface::isProvenDualInfeasible() const
{
	return (stsl_ == 2);
}

bool OsiFmpSolverInterface::isPrimalObjectiveLimitReached() const
{
  // *TEST*
  double limit = 0.0;
  getDblParam(OsiPrimalObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  ;//PAT EKKModel* model = getMutableModelPtr();
  ;//PAT if (ekk_getIprobstat(model)==6)
  ;//PAT   return false;
  ;//PAT const int lastalgo = ekk_lastAlgorithm(model);
  ;//PAT const double obj = ekk_getRobjvalue(model);
  ;//PAT const double maxmin = ekk_getRmaxmin(model);

  ;//PAT switch (lastalgo) {
  ;//PAT  case 0: // no simplex was needed
  ;//PAT    return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
  ;//PAT  case 2: // dual simplex
  ;//PAT    if (ekk_getIprobstat(model) == 0) // optimal
	;//PAT return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
  ;//PAT    return false;
  ;//PAT  case 1: // primal simplex
  ;//PAT    return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
  ;//PAT }
  return false; // fake return
}

bool OsiFmpSolverInterface::isDualObjectiveLimitReached() const
{
  // *TEST*
  double limit = 0.0;
  getDblParam(OsiDualObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  ;//PAT EKKModel* model = getMutableModelPtr();
  ;//PAT if (ekk_getIprobstat(model)==6)
  ;//PAT   return false;
  ;//PAT const int lastalgo = ekk_lastAlgorithm(model);
  ;//PAT const double obj = ekk_getRobjvalue(model);
  ;//PAT const double maxmin = ekk_getRmaxmin(model);
  
	;//PAT switch (lastalgo) {
  ;//PAT  case 0: // no simplex was needed
  ;//PAT    return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
  ;//PAT  case 1: // primal simplex
  ;//PAT    if (ekk_getIprobstat(model) == 0) // optimal
	;//PAT return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
  ;//PAT    return false;
  ;//PAT  case 2: // dual simplex
  ;//PAT    if (ekk_getIprobstat(model) != 0 && ekk_getIprobstat2(model) == 11)
	// over dual limit
	;//PAT return true;
  ;//PAT    return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
  ;//PAT }
  return false; // fake return
}

bool OsiFmpSolverInterface::isIterationLimitReached() const
{
  // *TEST*
  ;//PAT EKKModel* model = getMutableModelPtr();
  ;//PAT if (ekk_getIprobstat(model)==6)
  ;//PAT   return false;
  ;//PAT const int stat = ekk_getIprobstat(model);
  ;//PAT return (stat == 3);
	return true;
}

//#############################################################################
// WarmStart related methods
//#############################################################################

CoinWarmStart* OsiFmpSolverInterface::getWarmStart() const
{
  // *TEST*
  ;//PAT EKKModel* model = getMutableModelPtr();
  const int numcols = getNumCols();
  const int numrows = getNumRows();

  CoinWarmStartBasis* ws = new CoinWarmStartBasis;

  if(basis_==NULL)
  {
	  return ws;
  }

  ws->setSize(numcols, numrows);

  int i;

	for (i = 0; i < numrows; ++i) 
	{
		switch(basis_[i+1])
		{ 
		case  0:
			ws->setArtifStatus(i, CoinWarmStartBasis::basic);
			break;
		case  -1:
			ws->setArtifStatus(i, CoinWarmStartBasis::atLowerBound);
			break;
		case  1:
			ws->setArtifStatus(i, CoinWarmStartBasis::atUpperBound);
			break;
		}
  }

  for (i = 0; i < numcols; ++i) 
	{
    
		switch(basis_[numrows+i+1])
		{
    case  0:
      ws->setStructStatus(i, CoinWarmStartBasis::basic);
      break;
    case  -1: 
      ws->setStructStatus(i, CoinWarmStartBasis::atLowerBound);
      break;
    case  1:
      ws->setStructStatus(i, CoinWarmStartBasis::atUpperBound);
      break;
    }
  }
  
  return ws;
}

//-----------------------------------------------------------------------------

bool OsiFmpSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  // *TEST*
  const CoinWarmStartBasis* ws =
    dynamic_cast<const CoinWarmStartBasis*>(warmstart);

  if (! ws)
    return false;

  const int numcols = ws->getNumStructural();
  const int numrows = ws->getNumArtificial();

// The number of structural and artificial 
// must be consistent with the current model
	if(numcols!=nc_ || numrows!=nr_)
		return false;

	if(nc_ ==0 && nr_==0)
		return true;

	// Allocates a new basis vector if it isn't there
	if(basis_==NULL)
		basis_=new int[1+nc_+nr_];

	int i;

	basis_[0]=0; // The first row represents the objective

  for (i = 0; i < numrows; ++i) {
    switch (ws->getArtifStatus(i)) {
    case CoinWarmStartBasis::atUpperBound:
      basis_[i+1]=1;
      break;
    case CoinWarmStartBasis::atLowerBound:
      basis_[i+1]=-1;
      break;
    case CoinWarmStartBasis::basic:
      basis_[i+1]=0;
      break;
    case CoinWarmStartBasis::isFree: // FortMP cannot distinguish them
      basis_[i+1]=-1;
      break;
    }
  }

  for (i = 0; i < numcols; ++i) {
    switch (ws->getStructStatus(i)) {
    case CoinWarmStartBasis::atUpperBound:
      basis_[nr_+i+1]=1;
      break;
    case CoinWarmStartBasis::atLowerBound:
      basis_[nr_+i+1]=-1;
      break;
    case CoinWarmStartBasis::basic:
      basis_[nr_+i+1]=0;
      break;
    case CoinWarmStartBasis::isFree: // superbasic
			basis_[nr_+i+1]=-1;
      break;
    }
  }

	return true;
}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiFmpSolverInterface::markHotStart()
{
  // PAT todo
  
  delete ws_;
  ws_ = dynamic_cast<CoinWarmStartBasis*>(getWarmStart());
  //  ekk_startFastDualSimplex(model, itlim);

}

void OsiFmpSolverInterface::solveFromHotStart()
{
  ;//PAT EKKModel* model = getMutableModelPtr();
  // *TEST*
  ;//PAT ekk_setIiternum(model, 0);
  setWarmStart(ws_);
  //  ekk_fastDualSimplex(model);
  ;//PAT itlimOrig_ = ekk_getImaxiter(model);
  int itlim;
  OsiSolverInterface::getIntParam(OsiMaxNumIterationHotStart, itlim);
  ;//PAT ekk_setImaxiter(model, itlim);
  resolve();
  ;//PAT ekk_setImaxiter(model, itlimOrig_);
}

void OsiFmpSolverInterface::unmarkHotStart()
{
  // *TEST*
  ;//PAT EKKModel* model = getMutableModelPtr();
  //  ekk_stopFastDualSimplex(model);
  ;//PAT ekk_setImaxiter(model, itlimOrig_);
  delete ws_;
  ws_ = NULL;
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
// Get number of rows, columns and elements
//------------------------------------------------------------------
int OsiFmpSolverInterface::getNumCols() const
{
		return nc_;
}
int OsiFmpSolverInterface::getNumRows() const
{
		return nr_;
}
int OsiFmpSolverInterface::getNumElements() const
{
		return nz_;
}

//------------------------------------------------------------------
// Get pointer to rim vectors
//------------------------------------------------------------------  
const double * OsiFmpSolverInterface::getColLower() const
{
	return lob_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getColUpper() const
{
	return upb_;
}
//------------------------------------------------------------------
const char * OsiFmpSolverInterface::getRowSense() const
{
  extractSenseRhsRange();
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getRightHandSide() const
{
	extractSenseRhsRange();
  return rowrhs_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getRowRange() const
{
  extractSenseRhsRange();
  return rowrange_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getRowLower() const
{
	return lhs_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getRowUpper() const
{
	return rhs_;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getObjCoefficients() const
{
	return cost_;
}
//------------------------------------------------------------------
double OsiFmpSolverInterface::getObjSense() const
{
	if(isMax_)
		return -1;
	else
		return 1;
}

//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------
bool OsiFmpSolverInterface::isContinuous(int colNumber) const
{
		
  if(mitype_[colNumber]!=0) return false;
	return true;
}
//------------------------------------------------------------------
#if 0
bool OsiFmpSolverInterface::isInteger(
      int columnNumber ) const
{
  return !(isContinuous(columnNumber));
}
//------------------------------------------------------------------
bool OsiFmpSolverInterface::isBinary(
  int columnNumber ) const
{
		//PAT check...
		if(mitype_[columnNumber]==1)
				return true;
		return false;
}
//------------------------------------------------------------------
bool OsiFmpSolverInterface::isIntegerNonBinary(
  int columnNumber ) const
{
  if(mitype_[columnNumber]==2)
				return true;
		return false;
}
//------------------------------------------------------------------

bool OsiFmpSolverInterface::isFreeBinary(
  int columnNumber ) const
{
  if ( isContinuous(columnNumber)) 
		return false;
  
	if ((upb_[columnNumber]== 1) && (lob_[columnNumber]== 0))
		return true;
  else 
		return false;
}
#endif

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------
const CoinPackedMatrix * OsiFmpSolverInterface::getMatrixByRow() const
{
 if ( matrixByRow_ == NULL ) {

	 int i,p;
	 
	 CoinBigIndex *starts=new CoinBigIndex[nr_+1];
	 int *counts=new int[nr_];
	 int *indices=new int[nz_];
	 double *values=new double[nz_];

	 CoinFillN(indices,nz_,0);
	 CoinFillN(counts,nr_,0);
	 // counts the non zeroes for each row
	 for(i=0;i<nz_;i++)
		 counts[rowin_[i]-1]++;

	 // determines the start of each row
	 starts[0]=0;
	 for(i=1;i<nr_+1;i++)
		 starts[i]=starts[i-1]+counts[i-1];

	 // Reuse the counts as offsets
	 CoinFillN(counts,nr_,0);
	 for(i=0;i<nz_;i++)
	 {
		   p=starts[rowin_[i]-1]+counts[rowin_[i]-1];
			 indices[p]=colin_[i]-1;
			 values[p]=aij_[i];
			 counts[rowin_[i]-1]++;
	 }

    matrixByRow_ = new CoinPackedMatrix();
    matrixByRow_->copyOf(false,
												 getNumCols(), getNumRows(), getNumElements(),
												 values,indices,starts, NULL);
 
  delete[] starts;
	delete[] indices;
	delete[] counts;
	delete[] values;
 }
  return matrixByRow_;
}
//------------------------------------------------------------------
const CoinPackedMatrix * OsiFmpSolverInterface::getMatrixByCol() const
{
  if ( matrixByColumn_ == NULL ) 
	{
	 int i,p;
	 
	 CoinBigIndex *starts=new CoinBigIndex[nc_+1];
	 int *counts=new int[nc_];
	 int *indices=new int[nz_];
	 double *values=new double[nz_];

	 CoinFillN(indices,nz_,0);
	 CoinFillN(counts,nc_,0);
	 // counts the non zeroes for each column
	 for(i=0;i<nz_;i++)
		 counts[colin_[i]-1]++;

	 // determines the start of each column
	 starts[0]=0;
	 for(i=1;i<nc_+1;i++)
		 starts[i]=starts[i-1]+counts[i-1];

	 // Reuse the counts as offsets
	 CoinFillN(counts,nc_,0);
	 for(i=0;i<nz_;i++)
	 {
		   p=starts[colin_[i]-1]+counts[colin_[i]-1];
			 indices[p]=rowin_[i]-1;
			 values[p]=aij_[i];
			 counts[colin_[i]-1]++;

	 }

    matrixByColumn_= new CoinPackedMatrix();
    matrixByColumn_->copyOf(true,
												 getNumRows(),getNumCols(), getNumElements(),
												 values,indices,starts, NULL);
  
  delete[] starts;
	delete[] indices;
	delete[] counts;
	delete[] values;
	}
  return matrixByColumn_;
}
//------------------------------------------------------------------
// Get solver's value for infinity
//------------------------------------------------------------------
double OsiFmpSolverInterface::getInfinity() const
{
		;//PAT
  return 1e30;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

const double * OsiFmpSolverInterface::getColSolution() const
{
	if(sol_)
		return &sol_[nr_+1];
	else return NULL;
}
//------------------------------------------------------------------

const double * OsiFmpSolverInterface::getRowPrice() const
{
	if(dsl_)
		return &dsl_[1];
	else return NULL;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getReducedCost() const
{
	if(dsl_)
		return &dsl_[nr_+1];
	else return NULL;
}
//------------------------------------------------------------------
const double * OsiFmpSolverInterface::getRowActivity() const
{
	if(sol_)
		return &sol_[1];
	else return NULL;
}
//------------------------------------------------------------------
double OsiFmpSolverInterface::getObjValue() const
{

	if(stsl_!=0)
	{
		return obj_;
	}
	else
		return OsiSolverInterface::getObjValue();

}
//------------------------------------------------------------------
int OsiFmpSolverInterface::getIterationCount() const
{
	int ssxiter, ipmiter, mipnode;
	double ssxtim, ipmtim, xovtim, miptim;

	GTVALS(&ssxiter, &ipmiter, &mipnode,
         &ssxtim, &ipmtim, &xovtim, &miptim);

	return ssxiter;
}
//------------------------------------------------------------------
std::vector<double*> OsiFmpSolverInterface::getDualRays(int maxNumRays) const
{
  const int m = getNumRows();
  double* ray = new double[m];
;//PAT   const double* negray = ekk_rowaux(getMutableModelPtr());
;//PAT   for (int i = 0; i < m; ++i) {
;//PAT      ray[i] = -negray[i];
;//PAT   }
  return std::vector<double*>(1, ray);
}
//------------------------------------------------------------------
std::vector<double*> OsiFmpSolverInterface::getPrimalRays(int maxNumRays) const
{
  const int n = getNumCols();
  double* ray = new double[n];
;//PAT   CoinDisjointCopyN(ekk_colaux(getMutableModelPtr()), n, ray);
  return std::vector<double*>(1, ray);
}
//------------------------------------------------------------------
#if 0
OsiVectorInt
OsiFmpSolverInterface::getFractionalIndices(const double etol) const
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
OsiFmpSolverInterface::setObjCoeff( int elementIndex, double elementValue )
{
		assert(cost_ != NULL);
		cost_[elementIndex]=elementValue;

}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setColLower( int elementIndex, double elementValue )
{

		assert(lob_ != NULL);
		lob_[elementIndex]=forceIntoRange(elementValue, -FMP_INFINITY, FMP_INFINITY);
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setColUpper( int elementIndex, double elementValue )
{
		assert(upb_ != NULL);
		upb_[elementIndex]=forceIntoRange(elementValue, -FMP_INFINITY, FMP_INFINITY);
}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
		while (indexFirst != indexLast) 
		{
			const int iCol=*indexFirst++;
      lob_[iCol]= forceIntoRange(*boundList++, -FMP_INFINITY, FMP_INFINITY);
      upb_[iCol]= forceIntoRange(*boundList++, -FMP_INFINITY, FMP_INFINITY);
    }
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setRowLower( int i, double elementValue )
{
		assert(lhs_ != NULL);
		lhs_[i]=forceIntoRange(elementValue, -FMP_INFINITY, FMP_INFINITY);
	
		// Update cached data
		if (rowsense_ != NULL) 
		{
			assert ((rowrhs_ != NULL) && (rowrange_ != NULL));
			convertBoundToSense(lhs_[i], rhs_[i],
			rowsense_[i], rowrhs_[i], rowrange_[i]);
		}
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setRowUpper( int i, double elementValue )
{
		assert(rhs_ != NULL);
		rhs_[i]=forceIntoRange(elementValue, -FMP_INFINITY, FMP_INFINITY);

		// Update cached data
		if (rowsense_ != NULL) 
		{
			assert ((rowrhs_ != NULL) && (rowrange_ != NULL));
			convertBoundToSense(lhs_[i], rhs_[i],
			rowsense_[i], rowrhs_[i], rowrange_[i]);
		}
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setRowBounds( int i,
				     double lower, double upper )
{

		assert(rhs_ != NULL);
		assert(lhs_ != NULL);
		
		rhs_[i]=forceIntoRange(upper, -FMP_INFINITY, FMP_INFINITY);
		lhs_[i]=forceIntoRange(lower, -FMP_INFINITY, FMP_INFINITY);

		// Update cached data
		if (rowsense_ != NULL) 
		{
			assert ((rowrhs_ != NULL) && (rowrange_ != NULL));
			convertBoundToSense(lhs_[i], rhs_[i],
			rowsense_[i], rowrhs_[i], rowrange_[i]);
		}
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
	double lower, upper;
	convertSenseToBound(sense, rightHandSide, range, lower, upper);
	setRowBounds(i, lower, upper);
}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{

	assert(rhs_ != NULL);
	assert(lhs_ != NULL);

  const int len = indexLast - indexFirst;
 	while (indexFirst != indexLast) 
	{
    const int iRow=*indexFirst++;
    lhs_[iRow]= forceIntoRange(*boundList++, -FMP_INFINITY, FMP_INFINITY);
    rhs_[iRow]= forceIntoRange(*boundList++, -FMP_INFINITY, FMP_INFINITY);
  }
  
	if (rowsense_ != NULL) 
	{
    assert ((rowrhs_ != NULL) && (rowrange_ != NULL));
    indexFirst -= len;
    while (indexFirst != indexLast) 
		{
			const int iRow=*indexFirst++;
			convertBoundToSense(lhs_[iRow], rhs_[iRow],rowsense_[iRow], rowrhs_[iRow], rowrange_[iRow]);
    }
  }
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
 
/*
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
*/
}
//#############################################################################
void
OsiFmpSolverInterface::setContinuous(int index)
{
	assert(index>=0 && index<nc_);
	mitype_[index]=0;
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setInteger(int index)
{
	// FortMP distingushes between binary and general integer
	if(upb_[index]==1)
		mitype_[index]=1;
	else
		mitype_[index]=2;
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setContinuous(const int* indices, int len)
{
	int i;
	for(i=0;i<len;i++)
		mitype_[indices[i]]=0;
}
//-----------------------------------------------------------------------------
void
OsiFmpSolverInterface::setInteger(const int* indices, int len)
{
	int i;
	for(i=0;i<len;i++)
	{
		if(upb_[indices[i]]==1)
			mitype_[indices[i]]=1;
		else
			mitype_[indices[i]]=2;
	}
}
//#############################################################################
void OsiFmpSolverInterface::setObjSense(double s) 
{
	isMax_=!(s>0);
}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::setColSolution(const double * cs) 
{
	int i;
	if(sol_==NULL)
	{
		sol_=new double[1+nr_+nc_];
		CoinFillN(sol_,1+nr_+nc_,0.0);
	}
	for(i=0;i<nc_;i++)
		sol_[nr_+1+i]=cs[i];
	
}
//-----------------------------------------------------------------------------
void OsiFmpSolverInterface::setRowPrice(const double * rs) 
{
	int i;
	if(dsl_==NULL)
	{
		dsl_=new double[1+nr_+nc_];
		CoinFillN(dsl_,1+nr_+nc_,0.0);
	}
	for(i=0;i<nc_;i++)
		dsl_[1+i]=rs[i];
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiFmpSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  
	const CoinPackedVectorBase* pvec= &vec;
	addCols(1,&pvec,&collb,&colub,&obj);
	
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{

	freeCachedResults();

	double *nobj,*nlob,*nupb,*naij;
	int *nrowin,*ncolin,*nmitype;

  int i,j;
  int nz = 0;
  for (i = 0; i < numcols; ++i)
    nz += cols[i]->getNumElements();

	// allocates new vectors
	nrowin= new int[nz_+nz];
	ncolin= new int[nz_+nz];
	naij= new double[nz_+nz];

	nlob=new double[nc_+numcols];
	nupb=new double[nc_+numcols];
	nobj=new double[nc_+numcols];
	nmitype=new int[nc_+numcols];

	// copy old problem into new vectors
	memcpy(nrowin,rowin_,sizeof(int)*nz_);
	memcpy(ncolin,colin_,sizeof(int)*nz_);
	memcpy(naij,aij_,sizeof(double)*nz_);
	memcpy(nlob,lob_,sizeof(double)*nc_);
	memcpy(nupb,upb_,sizeof(double)*nc_);
	memcpy(nobj,cost_,sizeof(double)*nc_);
	memcpy(nmitype,mitype_,sizeof(int)*nc_);

	// release old copy
	delete[] rowin_;
	delete[] colin_;
	delete[] aij_;
	delete[] lob_;
	delete[] upb_;
	delete[] mitype_;
	delete[] cost_;

	// Set pointers to new vectors
	rowin_=nrowin;
	colin_=ncolin;
	aij_=naij;
	lob_=nlob;
	upb_=nupb;
	mitype_=nmitype;
	cost_=nobj;

	// Finally, copy contents
	// NOTE: The internal matrix is no longer column ordered...

	CoinDisjointCopyN(collb,numcols,lob_ + nc_);
	CoinDisjointCopyN(colub,numcols,upb_ + nc_);
	CoinDisjointCopyN(obj,numcols,cost_ + nc_);
	CoinFillN(mitype_ + nc_,numcols,0);

	nz = 0;
  for (i = 0; i < numcols; ++i) 
	{
    const CoinPackedVectorBase* col = cols[i];
    const int len = col->getNumElements();
		// copy the column indices
    CoinDisjointCopyN(col->getIndices(), len, rowin_ + nz_ +nz);
    // copy the non zero values
		for(j=0;j<len;j++)
			rowin_[nz_+j]++;
		CoinDisjointCopyN(col->getElements(), len, aij_ + nz_ +nz);
		CoinFillN(colin_ + nz_ + nz, len, nc_ + i +1);
    nz += len;
 }
	 // Update dimensions..
	 nz_+=nz;
	 nc_+=numcols;

	 // Allocate solution space...
	 AllocSolutionSpace();
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  // *TEST*
//  ekk_deleteColumns(getModelPtr(), num, columnIndices);
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
	const CoinPackedVectorBase *pvec = &vec;
	addRows(1,&pvec,&rowlb,&rowub);
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{

		const CoinPackedVectorBase *pvec = &vec;
	addRows(1,&pvec,&rowsen,&rowrhs,&rowrng);
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{

	freeCachedResults();

	double *naij,*nlhs,*nrhs;
	int *nrowin,*ncolin;

  int i,j;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

	// allocates new vectors
	nrowin= new int[nz_+nz];
	ncolin= new int[nz_+nz];
	naij= new double[nz_+nz];
	nlhs=new double[nr_+numrows];
	nrhs=new double[nr_+numrows];

	// copy old problem into new vectors
	memcpy(nrowin,rowin_,sizeof(int)*nz_);
	memcpy(ncolin,colin_,sizeof(int)*nz_);
	memcpy(naij,aij_,sizeof(double)*nz_);
	memcpy(nlhs,lhs_,sizeof(double)*nr_);
	memcpy(nrhs,rhs_,sizeof(double)*nr_);

	// release old copy
	delete[] rowin_;
	delete[] colin_;
	delete[] aij_;
	delete[] lhs_;
	delete[] rhs_;

	// Set pointers to new vectors
	rowin_=nrowin;
	colin_=ncolin;
	aij_=naij;
	lhs_=nlhs;
	rhs_=nrhs;

	// Finally, copy contents
	// NOTE: The matrix is no longer column ordered...

	nz = 0;
   for (i = 0; i < numrows; ++i) 
	 {
		  rhs_[nr_ + i]=rowub[i]; 
			lhs_[nr_ + i]=rowlb[i]; 
			const CoinPackedVectorBase* row = rows[i];
			const int len = row->getNumElements();
		
			// copy the column indices
			CoinDisjointCopyN(row->getIndices(), len, colin_ + nz_ +nz);
			// Add offset of 1 used by FortMP 
			for(j=0;j<len;j++)
				colin_[nz_+j]++;
			// copy the non zero values
			CoinDisjointCopyN(row->getElements(), len, aij_ + nz_ +nz);
			CoinFillN(rowin_ + nz_ + nz, len, nr_ + i + 1);
			nz += len;
		}
	 // Update dimensions..
	 nz_+=nz;
	 nr_+=numrows;
	 AllocSolutionSpace();
}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
 
	freeCachedResults();

	double *naij,*nlhs,*nrhs;
	int *nrowin,*ncolin;

  int i,j;
  int nz = 0;
  for (i = 0; i < numrows; ++i)
    nz += rows[i]->getNumElements();

	// allocates new vectors
	nrowin= new int[nz_+nz];
	ncolin= new int[nz_+nz];
	naij= new double[nz_+nz];
	nlhs=new double[nr_+numrows];
	nrhs=new double[nr_+numrows];

	// copy old problem into new vectors
	memcpy(nrowin,rowin_,sizeof(int)*nz_);
	memcpy(ncolin,colin_,sizeof(int)*nz_);
	memcpy(naij,aij_,sizeof(double)*nz_);
	memcpy(nlhs,lhs_,sizeof(double)*nr_);
	memcpy(nrhs,rhs_,sizeof(double)*nr_);

	// release old copy
	delete[] rowin_;
	delete[] colin_;
	delete[] aij_;
	delete[] lhs_;
	delete[] rhs_;

	// Set pointers to new vectors
	rowin_=nrowin;
	colin_=ncolin;
	aij_=naij;
	lhs_=nlhs;
	rhs_=nrhs;

	// Finally, copy contents
	// NOTE: The cached matrix is no longer column ordered...

	nz = 0;
   for (i = 0; i < numrows; ++i) {
    convertSenseToBound(rowsen[i], rowrhs[i], rowrng[i], rhs_[nr_ + i], lhs_[nr_ + i]);
    const CoinPackedVectorBase* row = rows[i];
    const int len = row->getNumElements();
		// copy the column indices
    CoinDisjointCopyN(row->getIndices(), len, colin_ + nz_ +nz);
		// Add offset of 1 used by FortMP 
			for(j=0;j<len;j++)
				colin_[nz_+j]++;

    // copy the non zero values
		CoinDisjointCopyN(row->getElements(), len, aij_ + nz_ +nz);
		CoinFillN(rowin_ + nz_ + nz, len, nr_ + i +1);
    nz += len;
 }
	 // Update dimensions..
	 nz_+=nz;
	 nr_+=numrows;

 	 AllocSolutionSpace();

}
//-----------------------------------------------------------------------------
void 
OsiFmpSolverInterface::deleteRows(const int num, const int * rowIndices)
{
/*
	This is going to be expensive!
*/


}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiFmpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
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
OsiFmpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
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
OsiFmpSolverInterface::loadProblem(const CoinPackedMatrix& matrix,
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
OsiFmpSolverInterface::assignProblem(CoinPackedMatrix*& matrix,
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
OsiFmpSolverInterface::loadProblem(const int numcols, const int numrows,
				   const CoinBigIndex * start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const double* rowlb, const double* rowub)
{

	int i;

	//freeCachedResults();

	nc_=numcols;
	nr_=numrows;
	nz_= start[numcols]; //TEST
	
	// allocates memory	
	if(nr_ && nc_ && nz_)
	{	
		lhs_=	new double[nr_];
		rhs_=	new double[nr_];
		lob_=	new double[nc_];
		upb_=	new double[nc_];
		cost_=	new double[nc_];
		aij_=	new double[nz_];
		rowin_=	new int[nz_];
		colin_=	new int[nz_];
		mitype_= new int[nc_];
		
		//all variables are continuous unless otherwise specified 
		CoinFillN(mitype_,nc_,0);
	}
	else
		return;

	// copy column information
	if(collb!=NULL)
		CoinDisjointCopyN(collb,nc_,lob_);
	else
		CoinFillN(lob_,nc_,0.0);
	
	if(colub!=NULL)
		CoinDisjointCopyN(colub,nc_,upb_);
	else
		CoinFillN(upb_,nc_,getInfinity());

	if(obj!=NULL)
		CoinDisjointCopyN(obj,nc_,cost_);
	else
		CoinFillN(cost_,nc_,0.0);

	// copy row information
	if(rowlb!=NULL)
		CoinDisjointCopyN(rowlb,nr_,lhs_);
	else
		CoinFillN(lhs_,nr_,-getInfinity());

	if(rowub!=NULL)
		CoinDisjointCopyN(rowub,nr_,rhs_);
	else
		CoinFillN(rhs_,nr_,getInfinity());

	int j=0;
	for(i=0;i<nz_;i++)
	{
		while(i>=start[j+1])
			j++;
		colin_[i]=j+1;
		rowin_[i]=index[i]+1;
		aij_[i]=value[i];
	}

	AllocSolutionSpace();
}
//-----------------------------------------------------------------------------

void
OsiFmpSolverInterface::loadProblem(const int numcols, const int numrows,
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
   for (int i = numrows-1; i >= 0; --i) 
	 {   
      convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
   }
  
	loadProblem(numcols, numrows, start, index,
							value,	collb, colub, obj, rowlb, rowub);

   delete[] rowlb;
   delete[] rowub;
}

//-----------------------------------------------------------------------------
// Read mps files
//-----------------------------------------------------------------------------

int OsiFmpSolverInterface::readMps(const char *filename,
				   const char *extension)
{
   return OsiSolverInterface::readMps(filename, extension);
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------

void OsiFmpSolverInterface::writeMps(const char * filename,
				     const char * extension,
				     double objSense) const
{
  std::string f(filename);
  std::string e(extension);
  std::string fullname;
  if (e!="") 
	{
    fullname = f + "." + e;
  }
	else 
	{ 
    fullname = f;
  }
  
	OsiFmpSolverInterface::writeMpsNative(fullname.c_str(),NULL,NULL);
}

//#############################################################################
// OSL specific public interfaces
//#############################################################################
#if 0
EKKModel * OsiFmpSolverInterface::getModelPtr()
{
  freeCachedResults();
  return getMutableModelPtr();
}

//------------------------------------------------------------------- 

void OsiFmpSolverInterface::incrementInstanceCounter()
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

void OsiFmpSolverInterface::decrementInstanceCounter()
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

unsigned int OsiFmpSolverInterface::getNumInstances()
{
  return numInstances_;
}
#endif

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiFmpSolverInterface::OsiFmpSolverInterface ()
  :
OsiSolverInterface()
{
	stsl_=0;
	nr_=0;						
	nc_=0;						
	nz_=0;						
	ns_=0;
	nq_=0;
	isMax_=0;
	objconstant_=0;
	obj_=0;
	// Reset pointers
	send_=NULL;
	sref_=NULL;
	sfun_=NULL;
	sbeg_=NULL;
	rowin_=NULL;
	colin_=NULL;
	aij_=NULL;
	mitype_=NULL;
	upb_=NULL;
  lob_=NULL;
	cost_=NULL;
	rhs_=NULL;
	lhs_=NULL;
	basis_=NULL;
	sol_=NULL;
	dsl_=NULL;
	modelname_="";
	rowsense_=NULL;
	rowrhs_=NULL;
	rowrange_=NULL;
  matrixByRow_=NULL;
  matrixByColumn_=NULL;
  ws_ = NULL;

	int tctn;
	BLDFMPC(&tctn);
	
	// PAT: The following doesn't work... need to find out why 
	//BFMPCBC(NULL,NULL,NULL,(Xph2cb) CBIterCountP2,NULL,NULL);

}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * OsiFmpSolverInterface::clone(bool CopyData) const
{
   if (CopyData) {
      return new OsiFmpSolverInterface(*this);
   } else {
      return new OsiFmpSolverInterface();
   }
}

#if 0
//-------------------------------------------------------------------
// Alternate Constructor 
//-------------------------------------------------------------------
OsiFmpSolverInterface::OsiFmpSolverInterface (EKKModel * m)
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
OsiFmpSolverInterface::OsiFmpSolverInterface (
                  const OsiFmpSolverInterface & source)
:
OsiSolverInterface(source)
{
	stsl_=0;
	nr_=0;						
	nc_=0;						
	nz_=0;						
	ns_=0;
	nq_=0;
	isMax_=0;
	objconstant_=0;
	obj_=0;

	// Reset pointers
	send_=NULL;
	sref_=NULL;
	sfun_=NULL;
	sbeg_=NULL;
	rowin_=NULL;
	colin_=NULL;
	aij_=NULL;
	mitype_=NULL;
	upb_=NULL;
  lob_=NULL;
	cost_=NULL;
	rhs_=NULL;
	lhs_=NULL;
	basis_=NULL;
	sol_=NULL;
	dsl_=NULL;
	modelname_="";
	rowsense_=NULL;
	rowrhs_=NULL;
	rowrange_=NULL;
  matrixByRow_=NULL;
  matrixByColumn_=NULL;
  ws_ = NULL;

	nc_=source.nc_;
	nr_=source.nr_;						
	nz_=source.nz_;		
	ns_=source.ns_;
	nq_=source.nq_;
	isMax_=source.isMax_;
	objconstant_=source.objconstant_;

	if(nr_)
	{
		lhs_=	new double[nr_];
		rhs_=	new double[nr_];
		// copy row information
		CoinDisjointCopyN(source.lhs_,nr_,lhs_);
		CoinDisjointCopyN(source.rhs_,nr_,rhs_);
	}

	if(nc_)
	{
		lob_=	new double[nc_];
		upb_=	new double[nc_];
		cost_=new double[nc_];
		mitype_= new int[nc_];
		// copy column information
		CoinDisjointCopyN(source.lob_,nc_,lob_);
		CoinDisjointCopyN(source.upb_,nc_,upb_);
		CoinDisjointCopyN(source.cost_,nc_,cost_);
		CoinDisjointCopyN(source.mitype_,nc_,mitype_);
	}
	if(ns_)
	{
		send_=	new int[ns_];
		sbeg_=	new int[ns_];
		sfun_=new int[ns_];
		sref_= new int[ns_];
		// copy column information
		CoinDisjointCopyN(source.send_,ns_,send_);
		CoinDisjointCopyN(source.sref_,ns_,sref_);
		CoinDisjointCopyN(source.sbeg_,ns_,sbeg_);
		CoinDisjointCopyN(source.sfun_,ns_,sfun_);
	}

	if(nz_)
	{
		aij_=	new double[nz_];
		rowin_=	new int[nz_];
		colin_=	new int[nz_];
		// copy matrix
		CoinDisjointCopyN(source.colin_,nz_,colin_);
		CoinDisjointCopyN(source.rowin_,nz_,rowin_);
		CoinDisjointCopyN(source.aij_,nz_,aij_);
	}	
	modelname_=source.modelname_;

	AllocSolutionSpace();
/*
	if(source.basis_!=NULL)
	{
			basis_=new int[1+nr_+nc_];
			CoinDisjointCopyN(source.basis_,1+nr_+nc_,basis_);
	}
	if(source.sol_!=NULL)
	{
			sol_=new double[1+nr_+nc_];
			CoinDisjointCopyN(source.sol_,1+nr_+nc_,sol_);
	}
	if(source.dsl_!=NULL)
	{
			dsl_=new double[1+nr_+nc_];
			CoinDisjointCopyN(source.dsl_,1+nr_+nc_,dsl_);
	}

	
	if(source.rowsense_!=NULL)
	{
	
			rowsense_ = new char[nr_];
			rowrhs_ = new double[nr_];
			rowrange_ = new double[nr_];

			CoinDisjointCopyN(source.rowsense_,nr_,rowsense_);
			CoinDisjointCopyN(source.rowrhs_,nr_,rowrhs_);
			CoinDisjointCopyN(source.rowrange_,nr_,rowrange_);
	}
	

	if(ws_!=NULL)
		setWarmStart(source.getWarmStart());
*/

	// PAT todo...

}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiFmpSolverInterface::~OsiFmpSolverInterface ()
{
	int tctn;
  gutsOfDestructor();
;//PAT   decrementInstanceCounter();
	UNBFMPC(&tctn);
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiFmpSolverInterface &
OsiFmpSolverInterface::operator=(const OsiFmpSolverInterface& rhs)
{
  if (this != &rhs) {    
    OsiSolverInterface::operator=(rhs);
    gutsOfDestructor();

	nc_=rhs.nc_;
	nr_=rhs.nr_;						
	nz_=rhs.nz_;		
	ns_=rhs.ns_;
	nq_=rhs.nq_;
	isMax_=rhs.isMax_;
	objconstant_=rhs.objconstant_;

	if(nr_)
	{
		lhs_=	new double[nr_];
		rhs_=	new double[nr_];
		// copy row information
		CoinDisjointCopyN(rhs.lhs_,nr_,lhs_);
		CoinDisjointCopyN(rhs.rhs_,nr_,rhs_);
	}

	if(nc_)
	{
		lob_=	new double[nc_];
		upb_=	new double[nc_];
		cost_=new double[nc_];
		mitype_= new int[nc_];
		// copy column information
		CoinDisjointCopyN(rhs.lob_,nc_,lob_);
		CoinDisjointCopyN(rhs.upb_,nc_,upb_);
		CoinDisjointCopyN(rhs.cost_,nc_,cost_);
		CoinDisjointCopyN(rhs.mitype_,nc_,mitype_);
	}
	if(ns_)
	{
		send_=	new int[ns_];
		sbeg_=	new int[ns_];
		sfun_=new int[ns_];
		sref_= new int[ns_];
		// copy column information
		CoinDisjointCopyN(rhs.send_,ns_,send_);
		CoinDisjointCopyN(rhs.sref_,ns_,sref_);
		CoinDisjointCopyN(rhs.sbeg_,ns_,sbeg_);
		CoinDisjointCopyN(rhs.sfun_,ns_,sfun_);
	}

	if(nz_)
	{
		aij_=	new double[nz_];
		rowin_=	new int[nz_];
		colin_=	new int[nz_];
		// copy matrix
		CoinDisjointCopyN(rhs.colin_,nz_,colin_);
		CoinDisjointCopyN(rhs.rowin_,nz_,rowin_);
		CoinDisjointCopyN(rhs.aij_,nz_,aij_);
	}	
	
	modelname_=rhs.modelname_;

	AllocSolutionSpace();
	/* should not copy cached information ... why?


	if(rhs.basis_!=NULL)
	{
			basis_=new int[1+nr_+nc_];
			CoinDisjointCopyN(rhs.basis_,1+nr_+nc_,basis_);
	}
	if(rhs.sol_!=NULL)
	{
			sol_=new double[1+nr_+nc_];
			CoinDisjointCopyN(rhs.sol_,1+nr_+nc_,sol_);
	}
	if(rhs.dsl_!=NULL)
	{
			dsl_=new double[1+nr_+nc_];
			CoinDisjointCopyN(rhs.dsl_,1+nr_+nc_,dsl_);
	}

	if(rhs.rowsense_!=NULL)
	{
	
			rowsense_ = new char[nr_];
			rowrhs_ = new double[nr_];
			rowrange_ = new double[nr_];

			CoinDisjointCopyN(rhs.rowsense_,nr_,rowsense_);
			CoinDisjointCopyN(rhs.rowrhs_,nr_,rowrhs_);
			CoinDisjointCopyN(rhs.rowrange_,nr_,rowrange_);
	}
	

	if(ws_!=NULL)
		setWarmStart(rhs.getWarmStart());
*/
	// PAT todo...



  }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiFmpSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
	freeCachedResults();
  const CoinPackedVector & row=rowCut.row();
	addRow(row,rowCut.lb(),rowCut.ub());
}

//-----------------------------------------------------------------------------

void OsiFmpSolverInterface::applyColCut( const OsiColCut & cc )
{

	const CoinPackedVector lbs = cc.lbs();
  const CoinPackedVector ubs = cc.ubs();
  int i;

  for ( i=0; i<lbs.getNumElements(); i++ ) 
	{
    if ( lbs.getElements()[i] > lob_[lbs.getIndices()[i]])
			lob_[lbs.getIndices()[i]]=lbs.getElements()[i];
  }
  for ( i=0; i<ubs.getNumElements(); i++ ) 
	{
    if ( ubs.getElements()[i] < upb_[ubs.getIndices()[i]])
			upb_[ubs.getIndices()[i]]=ubs.getElements()[i];
  }

}

//#############################################################################
// Static methods and data
//#############################################################################
#if 0
EKKContext * OsiFmpSolverInterface::getContextPtr()
{
  assert( contextPtr_ != NULL );
  return contextPtr_;
}

//------------------------------------------------------------------- 

EKKContext * OsiFmpSolverInterface::contextPtr_  = NULL;

//-------------------------------------------------------------------

unsigned int OsiFmpSolverInterface::numInstances_ = 0;
#endif
 
//#############################################################################
// Private methods
//#############################################################################

//-------------------------------------------------------------------
// Get pointer to EKKModel.
// const methods should use getMutableModelPtr().
// non-const methods should use getModelPtr().
//------------------------------------------------------------------- 
#if 0
EKKModel * OsiFmpSolverInterface::getMutableModelPtr() const
{
  if ( modelPtr_ == NULL ) {
    modelPtr_ = ekk_newModel(getContextPtr(),NULL);
  }
  return modelPtr_;
}

//-------------------------------------------------------------------


void OsiFmpSolverInterface::gutsOfCopy( const OsiFmpSolverInterface & source )
{
  modelPtr_ = source.modelPtr_;
}
#endif

//-------------------------------------------------------------------

void OsiFmpSolverInterface::gutsOfDestructor()
{  
	freeCachedResults();
	// Reset scalars
	
	nr_=0;						
	nc_=0;						
	nz_=0;						
	ns_=0;
	nq_=0;
	isMax_=0;
	objconstant_=0;
	modelname_="";
	// Deallocate memory
	if(send_!=NULL)
		delete[] send_;
	if(sref_!=NULL)
		delete[] sref_;
	if(sfun_!=NULL)
		delete[] sfun_;
	if(sbeg_!=NULL)
		delete[] sbeg_;
	if(rowin_!=NULL)
		delete[] rowin_;
	if(colin_!=NULL)
		delete[] colin_;
	if(aij_!=NULL)
		delete[] aij_;
	if(mitype_!=NULL)
		delete[] mitype_;
	if(upb_!=NULL)
		delete[] upb_;
  if(lob_!=NULL)
		delete[] lob_;
	if(cost_!=NULL)
		delete[] cost_;
	if(rhs_!=NULL)
		delete[] rhs_;
	if(lhs_!=NULL)
		delete[] lhs_;

	// Reset pointers
	send_=NULL;
	sref_=NULL;
	sfun_=NULL;
	sbeg_=NULL;
	rowin_=NULL;
	colin_=NULL;
	aij_=NULL;
	mitype_=NULL;
	upb_=NULL;
  lob_=NULL;
	cost_=NULL;
	rhs_=NULL;
	lhs_=NULL;
}

//------------------------------------------------------------------- 

void OsiFmpSolverInterface::freeCachedResults()
{  

	obj_=0;
	stsl_=0;
	if(basis_!=NULL)
		delete[] basis_;
	if(sol_!=NULL)
		delete[] sol_;
	if(dsl_!=NULL)
		delete[] dsl_;
	
	if(rowsense_!=NULL)
		delete[] rowsense_;
	if(rowrhs_!=NULL)
		delete[] rowrhs_;
	if(rowrange_!=NULL)
		delete[] rowrange_;
 
	if(matrixByRow_!=NULL)
		delete matrixByRow_;
	if(matrixByColumn_!=NULL)
		delete matrixByColumn_;
	if(ws_!=NULL)
		delete ws_;

	basis_=NULL;
	sol_=NULL;
	dsl_=NULL;
	rowsense_=NULL;
	rowrhs_=NULL;
	rowrange_=NULL;
  matrixByRow_=NULL;
  matrixByColumn_=NULL;
  ws_ = NULL;
}

//------------------------------------------------------------------
void OsiFmpSolverInterface::extractSenseRhsRange() const
{
		if (rowsense_ == NULL) 
		{
    // all three must be NULL
			assert ((rowrhs_ == NULL) && (rowrange_ == NULL));
   
			if ( nr_ != 0 ) 
			{
				rowsense_ = new char[nr_];
				rowrhs_ = new double[nr_];
				rowrange_ = new double[nr_];
      
				int i;
				for ( i=0; i<nr_; i++ ) 
				{
					convertBoundToSense(lhs_[i], rhs_[i], rowsense_[i], rowrhs_[i], rowrange_[i]);
				}
			}
	  }
}

//#############################################################################
// Resets as if default constructor
void 
OsiFmpSolverInterface::reset()
{
  setInitialData(); // clear base class
  gutsOfDestructor();
  itlimOrig_=9999999;
}

void OsiFmpSolverInterface::AllocSolutionSpace()
{
	// Only allocate space if the model is not empty
	if(nc_>0 || nr_>0)
	{
		sol_=new double[nc_+nr_+1];
		dsl_=new double[nc_+nr_+1];
		CoinFillN(dsl_,nc_+nr_+1,0.0);
		//sets variables to their lower bound 
		CoinCopyN(lob_,nc_,&sol_[nr_+1]);

	}

}

int OsiFmpSolverInterface::CBIterCountP2(int iter,double obj)
{
	return 1;
}


