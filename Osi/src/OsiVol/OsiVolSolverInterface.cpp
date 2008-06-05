// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cstdlib>
#include <numeric>
#include <cassert>
#include <cmath>

#include "CoinHelperFunctions.hpp"
#include "CoinWarmStartDual.hpp"

#include "OsiVolSolverInterface.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

//#######################################################################
//  Private helper methods
//#######################################################################

void
OsiVolSolverInterface::updateRowMatrix_() const 
{
   if (! rowMatrixCurrent_) {
      rowMatrix_.reverseOrderedCopyOf(colMatrix_);
      rowMatrixCurrent_ = true;
   }
}

void
OsiVolSolverInterface::updateColMatrix_() const
{
   if (! colMatrixCurrent_) {
      colMatrix_.reverseOrderedCopyOf(rowMatrix_);
      colMatrixCurrent_ = true;
   }
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::checkData_() const
{
   int i;
   for (i = getNumRows() - 1; i >= 0; --i) {
      if (rowlower_[i] > -1.0e20 &&
	  rowupper_[i] < 1.0e20 &&
	  rowlower_[i] != rowupper_[i])
	 throw CoinError("Volume algorithm is unable to handle ranged rows",
			"checkData_", "OsiVolSolverInterface");
   }

   for (i = getNumCols() - 1; i >= 0; --i)  {
      if (collower_[i] < -1.0e20 || colupper_[i] > 1.0e20)
	 throw CoinError("Volume algorithm is unable to handle infinite bounds",
			"checkData_", "OsiVolSolverInterface");
   }
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::compute_rc_(const double* u, double* rc) const 
{
  if (isZeroOneMinusOne_) {
    rowMatrixOneMinusOne_->timesMajor(u, rc);
  } else {
    rowMatrix_.transposeTimes(u, rc);
  }

  const int psize = getNumCols();
  std::transform(rc, rc+psize, objcoeffs_, rc, std::minus<double>());
  std::transform(rc, rc+psize, rc, std::negate<double>());
}

//#############################################################################

bool
OsiVolSolverInterface::test_zero_one_minusone_(const CoinPackedMatrix& m) const
{
  const int vecnum = m.getMajorDim();
  const double* elem = m.getElements();
  const int* start = m.getVectorStarts();
  const int* length = m.getVectorLengths();
  int i, j;
  for (i = 0; i < vecnum; ++i) {
    for (j = start[i] + length[i] - 1; j >= start[i]; --j) {
      const double val = elem[j];
      if (val != 1.0 && val != 0.0 && val != -1.0) {
	return false;
      }
    }
  }
  return true;
}

//-----------------------------------------------------------------------------

OsiVolSolverInterface::OsiVolMatrixOneMinusOne_::
OsiVolMatrixOneMinusOne_(const CoinPackedMatrix& m) {
  const int major = m.getMajorDim();
  const double* elem = m.getElements();
  const int* ind = m.getIndices();
  const int* start = m.getVectorStarts();
  const int* length = m.getVectorLengths();

  majorDim_ = major;
  minorDim_ = m.getMinorDim();

  plusSize_ = 0;
  minusSize_ = 0;
  int i, j;
  for (i = 0; i < major; ++i) {
    for (j = start[i] + length[i] - 1; j >= start[i]; --j) {
      const double val = elem[j];
      if (val == 1.0) {
	++plusSize_;
      } else if (val == -1.0) {
	++minusSize_;
      }
    }
  }
  if (plusSize_ > 0) {
    plusInd_ = new int[plusSize_];
  }
  if (minusSize_ > 0) {
    minusInd_ = new int[minusSize_];
  }
  plusStart_ = new int[major];
  plusLength_ = new int[major];
  minusStart_ = new int[major];
  minusLength_ = new int[major];

  plusSize_ = 0;
  minusSize_ = 0;
  for (i = 0; i < major; ++i) {
    plusStart_[i] = plusSize_;
    minusStart_[i] = minusSize_;
    const int last = start[i] + length[i];
    for (j = start[i]; j < last; ++j) {
      const double val = elem[j];
      if (val == 1.0) {
	plusInd_[plusSize_++] = ind[j];
      } else if (val == -1.0) {
	minusInd_[minusSize_++] = ind[j];
      }
    }
    plusLength_[i] = plusSize_ - plusStart_[i];
    minusLength_[i] = minusSize_ - minusStart_[i];
  }
  if (plusSize_ == 0) {
    delete[] plusStart_;    plusStart_  = NULL;
    delete[] plusLength_;   plusLength_ = NULL;
  }
  if (minusSize_ == 0) {
    delete[] minusStart_;    minusStart_  = NULL;
    delete[] minusLength_;   minusLength_ = NULL;
  }
}

//-----------------------------------------------------------------------------

OsiVolSolverInterface::OsiVolMatrixOneMinusOne_::
~OsiVolMatrixOneMinusOne_() {
  if (plusSize_ > 0) {
    delete[] plusInd_;      plusInd_    = NULL;
    delete[] plusStart_;    plusStart_  = NULL;
    delete[] plusLength_;   plusLength_ = NULL;
  }
  if (minusSize_ > 0) {
    delete[] minusInd_;      minusInd_    = NULL;
    delete[] minusStart_;    minusStart_  = NULL;
    delete[] minusLength_;   minusLength_ = NULL;
  }
}

//-----------------------------------------------------------------------------

void OsiVolSolverInterface::OsiVolMatrixOneMinusOne_::
timesMajor(const double* x, double* y) const
{
  memset(y, 0, minorDim_ * sizeof(double));
  int i;

  if (plusSize_ > 0 && minusSize_ > 0) {
    for (i = majorDim_ - 1; i >= 0; --i) {
      const double x_i = x[i];
      if (x_i != 0.0) {
	const int* vecInd = plusInd_ + plusStart_[i];
        int j;
	for ( j = plusLength_[i] - 1; j >= 0; --j)
	  y[vecInd[j]] += x_i;
	vecInd = minusInd_ + minusStart_[i];
	for ( j = minusLength_[i] - 1; j >= 0; --j)
	  y[vecInd[j]] -= x_i;
      }
    }
    return;
  }
  if (plusSize_ > 0) {
    for (i = majorDim_ - 1; i >= 0; --i) {
      const double x_i = x[i];
      if (x_i != 0.0) {
	const int* vecInd = plusInd_ + plusStart_[i];
	for (int j = plusLength_[i] - 1; j >= 0; --j)
	  y[vecInd[j]] += x_i;
      }
    }
    return;
  }
  if (minusSize_ > 0) {
    for (i = majorDim_ - 1; i >= 0; --i) {
      const double x_i = x[i];
      if (x_i != 0.0) {
	const int* vecInd = minusInd_ + minusStart_[i];
	for (int j = minusLength_[i] - 1; j >= 0; --j)
	  y[vecInd[j]] -= x_i;
      }
    }
    return;
  }
}

//#############################################################################

void
OsiVolSolverInterface::gutsOfDestructor_()
{
  rowMatrix_.clear();
  colMatrix_.clear();
  rowMatrixCurrent_ = true;
  colMatrixCurrent_ = true;

  delete[] colupper_;   colupper_ = 0;  
  delete[] collower_;	collower_ = 0;  
  delete[] continuous_; continuous_ = 0;
  delete[] rowupper_;	rowupper_ = 0; 
  delete[] rowlower_;	rowlower_ = 0; 
  delete[] rowsense_;	rowsense_ = 0; 
  delete[] rhs_;	rhs_ = 0;      
  delete[] rowrange_;	rowrange_ = 0; 
  delete[] objcoeffs_;	objcoeffs_ = 0;

  delete[] colsol_;	        colsol_ = 0;   
  delete[] rowprice_;	        rowprice_ = 0;
  delete[] rowpriceHotStart_;	rowpriceHotStart_ = 0;
  delete[] rc_;       	        rc_ = 0;
  delete[] lhs_;       	        lhs_ = 0;

  lagrangeanCost_ = 0.0;

  maxNumrows_ = 0;
  maxNumcols_ = 0;
}

//#############################################################################

void
OsiVolSolverInterface::rowRimAllocator_()
{
   rowupper_ = new double[maxNumrows_];
   rowlower_ = new double[maxNumrows_];
   rowsense_ = new char[maxNumrows_];
   rhs_      = new double[maxNumrows_];
   rowrange_ = new double[maxNumrows_];
   rowprice_ = new double[maxNumrows_];
   lhs_      = new double[maxNumrows_];
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::colRimAllocator_()
{
   colupper_  = new double[maxNumcols_];
   collower_  = new double[maxNumcols_];
   continuous_ = new bool[maxNumcols_];
   objcoeffs_ = new double[maxNumcols_];
   colsol_    = new double[maxNumcols_];
   rc_        = new double[maxNumcols_];
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::rowRimResize_(const int newSize)
{
   if (newSize > maxNumrows_) {
      double* rub   = rowupper_;
      double* rlb   = rowlower_;
      char*   sense = rowsense_;
      double* right = rhs_;
      double* range = rowrange_;
      double* dual  = rowprice_;
      double* left  = lhs_;
      maxNumrows_ = CoinMax(1000, (newSize * 5) / 4);
      rowRimAllocator_();
      const int rownum = getNumRows();
      CoinDisjointCopyN(rub  , rownum, rowupper_);
      CoinDisjointCopyN(rlb  , rownum, rowlower_);
      CoinDisjointCopyN(sense, rownum, rowsense_);
      CoinDisjointCopyN(right, rownum, rhs_);
      CoinDisjointCopyN(range, rownum, rowrange_);
      CoinDisjointCopyN(dual , rownum, rowprice_);
      CoinDisjointCopyN(left , rownum, lhs_);
      delete[] rub;
      delete[] rlb;
      delete[] sense;
      delete[] right;
      delete[] range;
      delete[] dual;
      delete[] left;
   }
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::colRimResize_(const int newSize)
{
   if (newSize > maxNumcols_) {
      double* cub = colupper_;
      double* clb = collower_;
      bool* cont  = continuous_;
      double* obj = objcoeffs_;
      double* sol = colsol_;
      double* rc  = rc_;
      maxNumcols_ = CoinMax(1000, (maxNumcols_ * 5) / 4);
      colRimAllocator_();
      const int colnum = getNumCols();
      CoinDisjointCopyN(cub , colnum, colupper_);
      CoinDisjointCopyN(clb , colnum, collower_);
      CoinDisjointCopyN(cont, colnum, continuous_);
      CoinDisjointCopyN(obj , colnum, objcoeffs_);
      CoinDisjointCopyN(sol , colnum, colsol_);
      CoinDisjointCopyN(rc  , colnum, rc_);
      delete[] cub;
      delete[] clb;
      delete[] cont;
      delete[] obj;
      delete[] sol;
      delete[] rc;
   }
}

//#############################################################################

void
OsiVolSolverInterface::convertBoundsToSenses_()
{
   for (int i = getNumRows() - 1; i >= 0; --i ) {
      convertBoundToSense(rowlower_[i], rowupper_[i],
			  rowsense_[i], rhs_[i], rowrange_[i]);
   }
}

//-----------------------------------------------------------------------------

void
OsiVolSolverInterface::convertSensesToBounds_()
{
   for (int i = getNumRows() - 1; i >= 0; --i) {
      convertSenseToBound(rowsense_[i], rhs_[i], rowrange_[i],
			  rowlower_[i], rowupper_[i]);
   }
}

//#############################################################################
//  Here are the routines implementing the virtual methods inherited from
//  VOL_user_hooks.
//#############################################################################

int
OsiVolSolverInterface::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
   compute_rc_(u.v, rc.v);
   return 0;
}

//-----------------------------------------------------------------------

int
OsiVolSolverInterface::solve_subproblem(const VOL_dvector& dual,
					const VOL_dvector& rc,
					double& lcost,
					VOL_dvector& x, VOL_dvector& v,
					double& pcost)
{
  int i;

  const int psize = x.size();
  for (i = 0; i < psize; ++i) {
    x[i] = (rc[i] >= 0.0) ? collower_[i] : colupper_[i];
  }
   
  const int dsize = v.size();
  lcost = (std::inner_product(rhs_, rhs_ + dsize, dual.v, 0.0) +
	   std::inner_product(x.v, x.v + psize, rc.v, 0.0) );

  if (isZeroOneMinusOne_) {
    colMatrixOneMinusOne_->timesMajor(x.v, v.v);
  } else {
    colMatrix_.times(x.v, v.v);
  }

  std::transform(v.v, v.v+dsize, rhs_, v.v, std::minus<double>());
  std::transform(v.v, v.v+dsize, v.v, std::negate<double>());

  pcost = std::inner_product(x.v, x.v + psize, objcoeffs_, 0.0);

  return 0;
}

//#############################################################################
// Solve methods
//#############################################################################

void 
OsiVolSolverInterface::initialSolve()
{
   // set every entry to 0.0 in the dual solution
   CoinFillN(rowprice_, getNumRows(), 0.0);
   resolve();
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::resolve()
{
  int i;
  
  checkData_();

  // Only one of these can do any work
  updateRowMatrix_();
  updateColMatrix_();

  const int dsize = getNumRows();
  const int psize = getNumCols();

  // Negate the objective coefficients if necessary
  if (objsense_ < 0) {
    std::transform(objcoeffs_, objcoeffs_+psize, objcoeffs_,
		   std::negate<double>());
  }

  // Set the lb/ub on the duals
  volprob_.dual_lb.allocate(dsize);
  volprob_.dual_ub.allocate(dsize);
  double * dlb = volprob_.dual_lb.v;
  double * dub = volprob_.dual_ub.v;
  for (i = 0; i < dsize; ++i) {
    dlb[i] = rowupper_[i] <  OsiVolInfinity ? -1.0e31 : 0.0;
    dub[i] = rowlower_[i] > -OsiVolInfinity ?  1.0e31 : 0.0;
  }
  volprob_.dsize = dsize;
  volprob_.psize = psize;
  
  // Set the dual starting point
  VOL_dvector& dsol = volprob_.dsol;
  dsol.allocate(dsize);
  std::transform(rowprice_, rowprice_+dsize, dsol.v,
		 std::bind2nd(std::multiplies<double>(), objsense_));

  // adjust the dual vector (if necessary) to be sure it's feasible
  double * dv = dsol.v;
  for (i = 0; i < dsize; ++i) {
    if (dv[i] < dlb[i]) {
      dv[i] = dlb[i];
    } else if (dv[i] > dub[i]) {
      dv[i] = dub[i];
    }
  }

  // If requested, check whether the matrix contains anything but 0/1/-1
#if 0
  isZeroOneMinusOne_ = false;
#else
  isZeroOneMinusOne_ = test_zero_one_minusone_(colMatrix_);
  if (isZeroOneMinusOne_) {
    colMatrixOneMinusOne_ = new OsiVolMatrixOneMinusOne_(colMatrix_);
    rowMatrixOneMinusOne_ = new OsiVolMatrixOneMinusOne_(rowMatrix_);
  }
#endif

  volprob_.solve(*this, true);

  // extract the solution

  // the lower bound on the objective value
  lagrangeanCost_ = objsense_ * volprob_.value;
  // the primal solution
  CoinDisjointCopyN(volprob_.psol.v, psize, colsol_);

  // Reset the objective coefficients if necessary
  if (objsense_ < 0) {
    std::transform(objcoeffs_, objcoeffs_ + psize, objcoeffs_,
		   std::negate<double>());
    // also, multiply the dual solution by -1
    std::transform(volprob_.dsol.v, volprob_.dsol.v+dsize, rowprice_,
		   std::negate<double>());
  } else {
    // now we just have to copy the dual
    CoinDisjointCopyN(volprob_.dsol.v, dsize, rowprice_);
  }

  // Compute the reduced costs
  compute_rc_(rowprice_, rc_);

  // Compute the left hand side (row activity levels)
  if (isZeroOneMinusOne_) {
    colMatrixOneMinusOne_->timesMajor(colsol_, lhs_);
  } else {
    colMatrix_.times(colsol_, lhs_);
  }

  if (isZeroOneMinusOne_) {
    delete colMatrixOneMinusOne_;
    colMatrixOneMinusOne_ = NULL;
    delete rowMatrixOneMinusOne_;
    rowMatrixOneMinusOne_ = NULL;
  }
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiVolSolverInterface::setIntParam(OsiIntParam key, int value)
{
  switch (key) {
  case OsiMaxNumIteration:
    if (value < 0)
      return false;
    volprob_.parm.maxsgriters = value;
    break;
  case OsiMaxNumIterationHotStart:
    if (value < 0)
      return false;
    OsiSolverInterface::setIntParam(key, value);
    break;
  case OsiLastIntParam:
    return false;
  default:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::setDblParam(OsiDblParam key, double value)
{
  switch (key) {
  case OsiDualObjectiveLimit:
    volprob_.parm.ubinit = value;
    break;
  case OsiPrimalObjectiveLimit: // not applicable
    return false;
  case OsiDualTolerance: // only ~0 is applicable, so accept only 1e-50 ...
    return (value == 1e-50);
  case OsiPrimalTolerance:
    if (value < 1e-04 || value > 1e-1)
      return false;
    volprob_.parm.primal_abs_precision = value;
    break;
  case OsiObjOffset: 
    return OsiSolverInterface::setDblParam(key, value);
  case OsiLastDblParam:
    return false;
  default:
      return false;
  }
  return true;
}


//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  bool retval=false;
  switch (key) {
  case OsiSolverName:
    return false;

  case OsiProbName:
    OsiSolverInterface::setStrParam(key,value);
    return retval = true;

  case OsiLastStrParam:
    return false;

  default:
      return false;
  }
  return false;
}

//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::getIntParam(OsiIntParam key, int& value) const
{
  switch (key) {
  case OsiMaxNumIteration:
    value = volprob_.parm.maxsgriters;
    break;
  case OsiMaxNumIterationHotStart:
    OsiSolverInterface::getIntParam(key, value);
    break;
  case OsiLastIntParam:
    return false;
  default:
    return false;
  }
  return true;
}

//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  switch (key) {
  case OsiDualObjectiveLimit:
    value = volprob_.parm.ubinit;
    break;
  case OsiPrimalObjectiveLimit: // not applicable
    return false;
  case OsiDualTolerance: // not applicable, but must return almost 0
    value = 1e-50;
    break;
  case OsiPrimalTolerance:
    value = volprob_.parm.primal_abs_precision; 
    break;
  case OsiObjOffset:
    OsiSolverInterface::getDblParam(key, value);
    break;
  case OsiLastDblParam:
    return false;
  default:
    return false;
  }
  return true;
}


//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  switch (key) {
  case OsiProbName:
    OsiSolverInterface::getStrParam(key, value);
    return true;
  case OsiSolverName:
    value = "vol";
    return true;
  case OsiLastStrParam:
    return false;
  default:
    return false;
  }
  return false;
}
//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool
OsiVolSolverInterface::isAbandoned() const
{
  // *THINK*: see the *THINK* in isProvenOptimal()
  return false;
}

bool
OsiVolSolverInterface::isProvenOptimal() const
{
  // if exited before reaching the iteration limit it is declared optimal
  // *THINK*: Granted, it can exit because the dual value is not improving.
  // *THINK*: Should that be "abandoned"? But then it'll be abandoned way too
  // *THINK*: frequently... 
  return (! isDualObjectiveLimitReached() &&
	  volprob_.iter() < volprob_.parm.maxsgriters);
}

bool
OsiVolSolverInterface::isProvenPrimalInfeasible() const
{
  // LL: *FIXME* : at the moment the volume can't detect primal infeasibility.
  // LL: *FIXME* : The dual will go to infinity.
  return false;
}

bool
OsiVolSolverInterface::isProvenDualInfeasible() const
{
  // LL: *FIXME* : at the moment the volume assumes dual feasibility...
  return false;
}

bool
OsiVolSolverInterface::isPrimalObjectiveLimitReached() const
{
  // The volume algorithm doesn't know anything about the primal; only the
  // dual is monotone
  return false;
}

bool
OsiVolSolverInterface::isDualObjectiveLimitReached() const
{
  return volprob_.parm.ubinit - volprob_.value < volprob_.parm.granularity;
}

bool
OsiVolSolverInterface::isIterationLimitReached() const
{
  return volprob_.iter() >= volprob_.parm.maxsgriters;
}

//#############################################################################
// WarmStart related methods
//#############################################################################

CoinWarmStart*
OsiVolSolverInterface::getEmptyWarmStart () const
{ return (dynamic_cast<CoinWarmStart *>(new CoinWarmStartDual())) ; }

CoinWarmStart*
OsiVolSolverInterface::getWarmStart() const
{
  return new CoinWarmStartDual(getNumRows(), rowprice_);
}

//-----------------------------------------------------------------------------

bool
OsiVolSolverInterface::setWarmStart(const CoinWarmStart* warmstart)
{
  const CoinWarmStartDual* ws =
    dynamic_cast<const CoinWarmStartDual*>(warmstart);

  if (! ws)
    return false;

  const int ws_size = ws->size();
  if (ws_size != getNumRows() && ws_size != 0) {
    throw CoinError("wrong dual warmstart size", "setWarmStart",
		   "OsiVolSolverInterface");
  }

  CoinDisjointCopyN(ws->dual(), ws_size, rowprice_);
  return true;
}

//#############################################################################
// HotStart related methods
//#############################################################################

void
OsiVolSolverInterface::markHotStart()
{
  delete[] rowpriceHotStart_;
  rowpriceHotStart_ = new double[getNumRows()];
  CoinDisjointCopyN(rowprice_, getNumRows(), rowpriceHotStart_);
}

void
OsiVolSolverInterface::solveFromHotStart()
{
   int itlimOrig = volprob_.parm.maxsgriters;
   getIntParam(OsiMaxNumIterationHotStart, volprob_.parm.maxsgriters);
   CoinDisjointCopyN(rowpriceHotStart_, getNumRows(), rowprice_);
   resolve();
   volprob_.parm.maxsgriters = itlimOrig;
}

void
OsiVolSolverInterface::unmarkHotStart()
{
  delete[] rowpriceHotStart_;
  rowpriceHotStart_ = NULL;
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

bool 
OsiVolSolverInterface::isContinuous(int colNumber) const
{
  assert( continuous_!=NULL );
  if ( continuous_[colNumber] ) return true;
  return false;
}

//-----------------------------------------------------------------------------

const CoinPackedMatrix *
OsiVolSolverInterface::getMatrixByRow() const {
   updateRowMatrix_();
   return &rowMatrix_;
}

//-----------------------------------------------------------------------

const CoinPackedMatrix *
OsiVolSolverInterface::getMatrixByCol() const {
   updateColMatrix_();
   return &colMatrix_;
}

//#############################################################################
// Problem information methods (results)
//#############################################################################

std::vector<double*> OsiVolSolverInterface::getDualRays(int maxNumRays) const
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getDualRays",
		 "OsiVolSolverInterface");
  return std::vector<double*>();
}
//------------------------------------------------------------------
std::vector<double*> OsiVolSolverInterface::getPrimalRays(int maxNumRays) const
{
  // *FIXME* : must write the method -LL
  throw CoinError("method is not yet written", "getPrimalRays",
		 "OsiVolSolverInterface");
  return std::vector<double*>();
}

//#############################################################################
// Problem modifying methods (rim vectors)
//#############################################################################

//-----------------------------------------------------------------------------
void OsiVolSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  while (indexFirst < indexLast) {
    const int ind = *indexFirst;
    collower_[ind] = boundList[0];
    colupper_[ind] = boundList[1];
    ++indexFirst;
    boundList += 2;
  }
}

//-----------------------------------------------------------------------------
void OsiVolSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  if (indexLast - indexFirst < getNumRows() / 3) {
    while (indexFirst < indexLast) {
      setRowBounds(*indexFirst, boundList[0], boundList[1]);
      ++indexFirst;
      boundList += 2;
    }
  } else {
    // it's better to convert everything at once
    while (indexFirst < indexLast) {
      const int ind = *indexFirst;
      rowlower_[ind] = boundList[0];
      rowupper_[ind] = boundList[1];
      ++indexFirst;
      boundList += 2;
    }
    convertBoundsToSenses_();
  }
}

//-----------------------------------------------------------------------------
void OsiVolSolverInterface::setRowSetTypes(const int* indexFirst,
					   const int* indexLast,
					   const char* senseList,
					   const double* rhsList,
					   const double* rangeList)
{
  if (indexLast - indexFirst < getNumRows() / 3) {
    while (indexFirst < indexLast) {
      setRowType(*indexFirst++, *senseList++, *rhsList++, *rangeList++);
    }
  } else {
    // it's better to convert everything at once
    while (indexFirst < indexLast) {
      const int ind = *indexFirst++;
      rowsense_[ind] = *senseList++;
      rhs_[ind] = *rhsList++;
      rowrange_[ind] = *rangeList++;
    }
    convertSensesToBounds_();
  }
}

//#############################################################################

void
OsiVolSolverInterface::setContinuous(int index)
{
  assert(continuous_ != NULL);
  if (index < 0 || index > getNumCols()) {
    throw CoinError("Index out of bound.", "setContinuous",
		   "OsiVolSolverInterface");
  }
  continuous_[index] = true;
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::setInteger(int index)
{
  assert(continuous_ != NULL);
  if (index < 0 || index > getNumCols()) {
    throw CoinError("Index out of bound.", "setContinuous",
		   "OsiVolSolverInterface");
  }
  continuous_[index] = false;
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::setContinuous(const int* indices, int len)
{
  assert(continuous_ != NULL);
  const int colnum = getNumCols();
  int i;

  for (i = len - 1; i >= 0; --i) {
    if (indices[i] < 0 || indices[i] > colnum) {
      throw CoinError("Index out of bound.", "setContinuous",
		     "OsiVolSolverInterface");
    }
  }
  
  for (i = len - 1; i >= 0; --i) {
    continuous_[indices[i]] = true;
  }
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::setInteger(const int* indices, int len)
{
  assert(continuous_ != NULL);
  const int colnum = getNumCols();
  int i;

  for (i = len - 1; i >= 0; --i) {
    if (indices[i] < 0 || indices[i] > colnum) {
      throw CoinError("Index out of bound.", "setContinuous",
		     "OsiVolSolverInterface");
    }
  }
  
  for (i = len - 1; i >= 0; --i) {
    continuous_[indices[i]] = false;
  }
}

//#############################################################################

void
OsiVolSolverInterface::setColSolution(const double *colsol)
{
   CoinDisjointCopyN(colsol, getNumCols(), colsol_);
  // Compute the left hand side (row activity levels)
  if (isZeroOneMinusOne_) {
    colMatrixOneMinusOne_->timesMajor(colsol_, lhs_);
  } else {
    colMatrix_.times(colsol_, lhs_);
  }
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::setRowPrice(const double *rowprice)
{
   CoinDisjointCopyN(rowprice, getNumRows(), rowprice_);
   compute_rc_(rowprice_, rc_);
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################

void 
OsiVolSolverInterface::addCol(const CoinPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  const int colnum = getNumCols();
  colRimResize_(colnum + 1);
  collower_[colnum]   = collb;
  colupper_[colnum]   = colub;
  objcoeffs_[colnum]  = obj;
  continuous_[colnum] = true;
  colsol_[colnum]     = fabs(collb)<fabs(colub) ? collb : colub;
  rc_[colnum]         = 0.0;

  updateColMatrix_();
  colMatrix_.appendCol(vec);
  rowMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::addCols(const int numcols,
			       const CoinPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  if (numcols > 0) {
    const int colnum = getNumCols();
    colRimResize_(colnum + numcols);
    CoinDisjointCopyN(collb, numcols, collower_ + colnum);
    CoinDisjointCopyN(colub, numcols, colupper_ + colnum);
    CoinDisjointCopyN(obj, numcols, objcoeffs_ + colnum);
    CoinFillN(continuous_ + colnum, numcols, true);
    int c;
    for ( c=0; c<numcols; c++ ) {
      if ( fabs(collb[c]) < fabs(colub[c]) ) {
        colsol_[colnum+c] = collb[c];
      }
      else {
        colsol_[colnum+c] = colub[c];
      }
    }
    //CoinFillN(colsol_     + colnum, numcols, 0.0);
    CoinFillN(rc_         + colnum, numcols, 0.0);

    updateColMatrix_();
    colMatrix_.appendCols(numcols, cols);
    rowMatrixCurrent_ = false;
  }
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  if (num > 0) {
    int * delPos = new int[num];
    CoinDisjointCopyN(columnIndices, num, delPos);
    std::sort(delPos, delPos + num);
    const int delNum = std::unique(delPos, delPos + num) - delPos;

    const int colnum = getNumCols();
    CoinDeleteEntriesFromArray(collower_, collower_ + colnum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(colupper_, colupper_ + colnum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(objcoeffs_, objcoeffs_ + colnum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(continuous_, continuous_ + colnum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(colsol_, colsol_ + colnum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rc_, rc_ + colnum,
			       delPos, delPos + delNum);

    updateColMatrix_();
    colMatrix_.deleteCols(delNum, delPos);
    rowMatrixCurrent_ = false;
  }
}


//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  const int rownum = getNumRows();
  rowRimResize_(rownum + 1);
  rowlower_[rownum] = rowlb;
  rowupper_[rownum] = rowub;
  convertBoundToSense(rowlb, rowub,
		      rowsense_[rownum], rhs_[rownum], rowrange_[rownum]);
  rowprice_[rownum] = 0.0;
  lhs_[rownum] = 0.0;

  updateRowMatrix_();
  rowMatrix_.appendRow(vec);
  colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::addRow(const CoinPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  const int rownum = getNumRows();
  rowRimResize_(rownum + 1);
  rowsense_[rownum] = rowsen;
  rhs_[rownum] = rowrhs;
  rowrange_[rownum] = rowrng;
  convertSenseToBound(rowsen, rowrhs, rowrng,
		      rowlower_[rownum], rowupper_[rownum]);
  rowprice_[rownum] = 0.0;
  lhs_[rownum] = 0.0;

  updateRowMatrix_();
  rowMatrix_.appendRow(vec);
  colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  if (numrows > 0) {
    const int rownum = getNumRows();
    rowRimResize_(rownum + numrows);
    CoinDisjointCopyN(rowlb, numrows, rowlower_ + rownum);
    CoinDisjointCopyN(rowub, numrows, rowupper_ + rownum);
    for (int i = rownum + numrows - 1; i >= rownum; --i) {
      convertBoundToSense(rowlower_[i], rowupper_[i],
			  rowsense_[i], rhs_[i], rowrange_[i]);
    }
    CoinFillN(rowprice_ + rownum, numrows, 0.0);
    CoinFillN(lhs_      + rownum, numrows, 0.0);

    updateRowMatrix_();
    rowMatrix_.appendRows(numrows, rows);
    colMatrixCurrent_ = false;
  }
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::addRows(const int numrows,
			       const CoinPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  if (numrows > 0) {
    const int rownum = getNumRows();
    rowRimResize_(rownum + numrows);
    CoinDisjointCopyN(rowsen, numrows, rowsense_ + rownum);
    CoinDisjointCopyN(rowrhs, numrows, rhs_ + rownum);
    CoinDisjointCopyN(rowrng, numrows, rowrange_ + rownum);
    for (int i = rownum + numrows - 1; i >= rownum; --i) {
      convertSenseToBound(rowsense_[i], rhs_[i], rowrange_[i],
			  rowlower_[i], rowupper_[i]);
    }
    CoinFillN(rowprice_ + rownum, numrows, 0.0);
    CoinFillN(lhs_      + rownum, numrows, 0.0);

    updateRowMatrix_();
    rowMatrix_.appendRows(numrows, rows);
    colMatrixCurrent_ = false;
  }
}

//-----------------------------------------------------------------------------

void 
OsiVolSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  if (num > 0) {
    int * delPos = new int[num];
    CoinDisjointCopyN(rowIndices, num, delPos);
    std::sort(delPos, delPos + num);
    const int delNum = std::unique(delPos, delPos + num) - delPos;

    const int rownum = getNumRows();
    CoinDeleteEntriesFromArray(rowlower_, rowlower_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rowupper_, rowupper_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rowsense_, rowsense_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rowrange_, rowrange_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rhs_, rhs_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(rowprice_, rowprice_ + rownum,
			       delPos, delPos + delNum);
    CoinDeleteEntriesFromArray(lhs_, lhs_ + rownum,
			       delPos, delPos + delNum);

    updateRowMatrix_();
    rowMatrix_.deleteRows(delNum, delPos);
    colMatrixCurrent_ = false;

    delete[] delPos;
  }
}
   
//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

OsiVolSolverInterface::OsiVolSolverInterface () :
   rowMatrixCurrent_(true),
   rowMatrix_(),
   colMatrixCurrent_(true),
   colMatrix_(),
   isZeroOneMinusOne_(false),

   colupper_(0),
   collower_(0),
   continuous_(0),
   rowupper_(0),
   rowlower_(0),
   rowsense_(0),
   rhs_(0),
   rowrange_(0),

   objcoeffs_(0),
   objsense_(1.0),
   
   colsol_(0),
   rowprice_(0),
   rc_(0),
   lhs_(0),
   lagrangeanCost_(0.0),

   rowpriceHotStart_(0),

   maxNumrows_(0),
   maxNumcols_(0),

   volprob_()
{
   volprob_.parm.granularity = 0.0;
}

//-----------------------------------------------------------------------

OsiSolverInterface *
OsiVolSolverInterface::clone(bool copyData) const {
  return copyData ?
    new OsiVolSolverInterface(*this) :
    new OsiVolSolverInterface();
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::OsiVolSolverInterface(const OsiVolSolverInterface& x) :
   rowMatrixCurrent_(true),
   rowMatrix_(),
   colMatrixCurrent_(true),
   colMatrix_(),
   isZeroOneMinusOne_(false),

   colupper_(0),
   collower_(0),
   continuous_(0),
   rowupper_(0),
   rowlower_(0),
   rowsense_(0),
   rhs_(0),
   rowrange_(0),

   objcoeffs_(0),
   objsense_(1.0),
   
   colsol_(0),
   rowprice_(0),
   rc_(0),
   lhs_(0),
   lagrangeanCost_(0.0),

   rowpriceHotStart_(0),

   maxNumrows_(0),
   maxNumcols_(0),

   volprob_()
{
   operator=(x);
   volprob_.parm.granularity = 0.0;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface&
OsiVolSolverInterface::operator=(const OsiVolSolverInterface& rhs)
{
   if (&rhs == this)
      return *this;

   OsiSolverInterface::operator=(rhs);
   gutsOfDestructor_();

   rowMatrixCurrent_ = rhs.rowMatrixCurrent_;
   if (rowMatrixCurrent_)
      rowMatrix_ = rhs.rowMatrix_;
   colMatrixCurrent_ = rhs.colMatrixCurrent_;
   if (colMatrixCurrent_)
      colMatrix_ = rhs.colMatrix_;

   if (rhs.maxNumrows_) {
      maxNumrows_ = rhs.maxNumrows_;
      rowRimAllocator_();
      const int rownum = getNumRows();
      CoinDisjointCopyN(rhs.rowupper_, rownum, rowupper_);
      CoinDisjointCopyN(rhs.rowlower_, rownum, rowlower_);
      CoinDisjointCopyN(rhs.rowsense_, rownum, rowsense_);
      CoinDisjointCopyN(rhs.rhs_, rownum, rhs_);
      CoinDisjointCopyN(rhs.rowrange_, rownum, rowrange_);
      CoinDisjointCopyN(rhs.rowprice_, rownum, rowprice_);
      CoinDisjointCopyN(rhs.lhs_, rownum, lhs_);
   }
   if (rhs.maxNumcols_) {
      maxNumcols_ = rhs.maxNumcols_;
      colRimAllocator_();
      const int colnum = getNumCols();
      CoinDisjointCopyN(rhs.colupper_, colnum, colupper_);
      CoinDisjointCopyN(rhs.collower_, colnum, collower_);
      CoinDisjointCopyN(rhs.continuous_, colnum, continuous_);
      CoinDisjointCopyN(rhs.objcoeffs_, colnum, objcoeffs_);
      CoinDisjointCopyN(rhs.colsol_, colnum, colsol_);
      CoinDisjointCopyN(rhs.rc_, colnum, rc_);
   }
   volprob_.parm.granularity = 0.0;
   return *this;
}

//-----------------------------------------------------------------------

OsiVolSolverInterface::~OsiVolSolverInterface ()
{
   gutsOfDestructor_();
}

//#############################################################################
// Applying cuts
//#############################################################################

void
OsiVolSolverInterface::applyRowCut(const OsiRowCut& rc)
{
   const int rownum = getNumRows();
   const double lb = rc.lb();
   const double ub = rc.ub();
   rowRimResize_(rownum + 1);
   rowprice_[rownum] = 0.0;
   rowlower_[rownum] = lb;
   rowupper_[rownum] = ub;
   convertBoundToSense(lb, ub,
		       rowsense_[rownum], rhs_[rownum], rowrange_[rownum]);

   updateRowMatrix_();
   rowMatrix_.appendRow(rc.row());
   colMatrixCurrent_ = false;
}

//-----------------------------------------------------------------------

void
OsiVolSolverInterface::applyColCut(const OsiColCut& cc)
{
   int i;

   const double* lb_elem = cc.lbs().getElements();
   const int* lb_ind = cc.lbs().getIndices();
   for (i = cc.lbs().getNumElements() - 1; i >= 0; --i) {
      collower_[lb_ind[i]] = CoinMax(collower_[lb_ind[i]], lb_elem[i]);
   }
   
   const double* ub_elem = cc.ubs().getElements();
   const int* ub_ind = cc.ubs().getIndices();
   for (i = cc.ubs().getNumElements() - 1; i >= 0; --i) {
      colupper_[ub_ind[i]] = CoinMin(colupper_[ub_ind[i]], ub_elem[i]);
   }
}

//#############################################################################
