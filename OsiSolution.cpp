// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <math.h>

#include "CoinHelperFunctions.hpp"
#include "OsiSolution.hpp"
#include "OsiFactorization.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiIndexedVector.hpp"
#include "OsiWarmStartBasis.hpp"
#include <cassert>
#include <cmath>
#include <cfloat>

const bool True = true;
const bool False = false;

#include <math.h>
#include <string>
#include <stdio.h>
#include <iostream>

//#############################################################################

OsiSimplexModel::OsiSimplexModel () :

  primalTolerance_(1.0e-7),
  dualTolerance_(1.0e-7),
  optimizationDirection_(1),
  numberRefinements_(0),
  numberRows_(0),
  numberColumns_(0),
  pivotVariable_(NULL),
  factorization_(NULL),
  rowSense_(NULL),
  rowLower_(NULL),
  rowUpper_(NULL),
  objective_(NULL),
  colLower_(NULL),
  colUpper_(NULL),
  matrix_(NULL),
  columnStart_(NULL),
  columnLength_(NULL),
  row_(NULL),
  element_(NULL),
  basis_(NULL)
{
  description_.rhsSense=0;
  description_.rhsMine=0;
  description_.matrixType=0;
  description_.matrixMine=0;
  description_.columnMine=0;
  description_.basisMine=0;
  description_.factorizationMine=0;
  primalSolution_[0]=NULL;
  primalSolution_[1]=NULL;
  dualSolution_[0]=NULL;
  dualSolution_[1]=NULL;
  factorization_=NULL;
}

//-----------------------------------------------------------------------------

OsiSimplexModel::~OsiSimplexModel ()
{
  if (description_.factorizationMine) 
    delete factorization_;
  if (description_.basisMine) 
    delete basis_;
  gutsOfDelete();
}
void OsiSimplexModel::gutsOfDelete()
{
  delete [] primalSolution_[0];
  delete [] primalSolution_[1];
  delete [] dualSolution_[0];
  delete [] dualSolution_[1];
  delete [] pivotVariable_;
  if (description_.rhsMine) {
    delete [] rowSense_;
    delete [] rowLower_;
    delete [] rowUpper_;
  }
  description_.rhsMine=0;
  if (description_.columnMine) {
    delete [] colLower_;
    delete [] colUpper_;
    delete [] objective_;
  }
  description_.columnMine=0;
  if (description_.matrixMine) {
    if (description_.matrixType) {
      delete [] columnStart_;
      delete [] columnLength_;
      delete [] row_;
      delete [] element_;
    } else {
      delete matrix_;
    }
  }
  description_.matrixMine=0;
}
//#############################################################################
void OsiSimplexModel::setPrimalTolerance( double value) 
{
  if (value>0.0&&value<1.0e10)
    primalTolerance_=value;
}
void OsiSimplexModel::setDualTolerance( double value) 
{
  if (value>0.0&&value<1.0e10)
    dualTolerance_=value;
}
void OsiSimplexModel::setnumberRefinements( int value) 
{
  if (value>=0&&value<10)
    numberRefinements_=value;
}
void OsiSimplexModel::setOptimizationDirection( int value) 
{
  if (value>=-1&&value<=1)
    optimizationDirection_=value;
}
void
OsiSimplexModel::borrowModel (  const OsiPackedMatrix& matrix,
		     const double* collb, const double* colub,   
		     const double* obj,
		     const double* rowlb, const double* rowub)
{
  gutsOfDelete();
  numberRows_=matrix.getNumRows();
  numberColumns_=matrix.getNumCols();
  primalSolution_[0]=new double[numberRows_];
  primalSolution_[1]=new double[numberColumns_];
  dualSolution_[0]=new double[numberRows_];
  dualSolution_[1]=new double[numberColumns_];
  pivotVariable_=new int[numberRows_];

  CoinFillN(dualSolution_[0],numberRows_,0.0);
  CoinFillN(dualSolution_[1],numberColumns_,0.0);
  int iRow,iColumn;

  rowLower_=rowlb;
  rowUpper_=rowub;
  objective_=obj;
  colLower_=collb;
  colUpper_=colub;
  matrix_=&matrix;
  description_.rhsSense=0;
  description_.matrixType=0;
  // set default solution
  for (iRow=0;iRow<numberRows_;iRow++) {
    if (rowLower_[iRow]>0.0) {
      primalSolution_[0][iRow]=rowLower_[iRow];
    } else if (rowUpper_[iRow]<0.0) {
      primalSolution_[0][iRow]=rowUpper_[iRow];
    } else {
      primalSolution_[0][iRow]=0.0;
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (colLower_[iColumn]>0.0) {
      primalSolution_[1][iColumn]=colLower_[iColumn];
    } else if (colUpper_[iColumn]<0.0) {
      primalSolution_[1][iColumn]=colUpper_[iColumn];
    } else {
      primalSolution_[1][iColumn]=0.0;
    }
  }
}
// Factorizes a problem 
int 
OsiSimplexModel::factorize ( OsiFactorization & factorization, 
			     const OsiWarmStartBasis & basis)
{
  if (description_.factorizationMine) 
    delete factorization_;
  if (description_.basisMine) 
    delete basis_;
  description_.factorizationMine=0;
  description_.basisMine=0;
  factorization_=&factorization;
  basis_=&basis;
  
  int * rowIsBasic = new int[numberRows_];
  int * columnIsBasic = new int[numberColumns_];
  int i;
  int numberBasic=0;
  for (i=0;i<numberRows_;i++) {
    if (basis.getArtifStatus(i) == OsiWarmStartBasis::basic) {
      rowIsBasic[i]=1;
      numberBasic++;
    } else {
      rowIsBasic[i]=-1;
    }
  }
  for (i=0;i<numberColumns_;i++) {
    if (basis.getStructStatus(i) == OsiWarmStartBasis::basic) {
      columnIsBasic[i]=1;
      numberBasic++;
    } else {
      columnIsBasic[i]=-1;
    }
  }
  assert (numberBasic==numberRows_);
  int status=-1;
  while (status<0) {
    if (!description_.matrixType) {
      status =  factorization_->factorize(*matrix_,rowIsBasic,
					 columnIsBasic);
    } else {
      if (description_.matrixType==2) {
	status =  factorization_->factorize(numberRows_,numberColumns_,
					   columnStart_,columnLength_,
					   row_,element_,
					   rowIsBasic, columnIsBasic);
      } else {
	int * length = new int [numberColumns_];
	for (i=0;i<numberColumns_;i++) {
	  length[i] = columnStart_[i+1]-columnStart_[i];
	}
	status =  factorization_->factorize(numberRows_,numberColumns_,
					   columnStart_,length,
					   row_,element_,
		   rowIsBasic, columnIsBasic);

      }
    }
    if (status==-99) {
      // get more memory
      factorization_->areaFactor(1.5*factorization_->areaFactor());
    } else {
      break;
    }
  }
  if (!status) {
    // do pivot information
    numberBasic=0;
    for (i=0;i<numberRows_;i++) {
      if (basis.getArtifStatus(i) == OsiWarmStartBasis::basic) {
	pivotVariable_[rowIsBasic[i]]=i+numberRows_;
	numberBasic++;
      }
    }
    for (i=0;i<numberColumns_;i++) {
      if (basis.getStructStatus(i) == OsiWarmStartBasis::basic) {
	pivotVariable_[columnIsBasic[i]]=i;
	numberBasic++;
      }
    }
  }
  return status;
}

void
OsiSimplexModel::getRowBound(int iRow, double& lower, double& upper) const
{
  lower=-DBL_MAX;
  upper=DBL_MAX;
  if (!description_.rhsSense) {
    if (rowUpper_)
      upper=rowUpper_[iRow];
    if (rowLower_)
      lower=rowLower_[iRow];
  } else {
    double right=0.0,range=0.0;
    char sense = 'G';
    if (rowLower_)
      right = rowLower_[iRow];
    if (rowUpper_)
      range = rowUpper_[iRow];
    if (rowSense_)
      sense = rowSense_[iRow];
    switch (sense) {
    case 'E':
      lower = upper = right;
      break;
    case 'L':
      lower = -DBL_MAX;
      upper = right;
      break;
    case 'G':
      lower = right;
      upper = DBL_MAX;
      break;
    case 'R':
      lower = right - range;
      upper = right;
      break;
    case 'N':
      lower = -DBL_MAX;
      upper = DBL_MAX;
      break;
    }
  }
}
// Factorizes a problem and sets solution (nonbasic)
int 
OsiSimplexModel::factorizeAndSet ( OsiFactorization & factorization, 
			     const OsiWarmStartBasis & basis)
{
  int iRow,iColumn;
  int numberSuperBasic=0;
  int numberFree=0;
						
  for (iRow=0;iRow<numberRows_;iRow++) {
    double lower,upper;
    getRowBound(iRow, lower, upper);
    if (basis.getArtifStatus(iRow) == OsiWarmStartBasis::atUpperBound) {
      assert (upper<1.0e20);
      primalSolution_[0][iRow]=upper;
    } else if (basis.getArtifStatus(iRow) == 
	       OsiWarmStartBasis::atLowerBound) {
      assert (lower>-1.0e20);
      primalSolution_[0][iRow]=lower;
    } else if (basis.getArtifStatus(iRow) == OsiWarmStartBasis::isFree) {
      if (upper>1.0e20&&lower<-1.0e20) {
	numberFree++;
	primalSolution_[0][iRow]=0.0;
      } else {
	numberSuperBasic++;
	if (lower>0.0) {
	  primalSolution_[0][iRow]=lower;
	} else if (upper<0.0) {
	  primalSolution_[0][iRow]=upper;
	} else {
	  primalSolution_[0][iRow]=0.0;
	}
      }
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double lower=0.0,upper=DBL_MAX;
    if (colLower_)
      lower=colLower_[iColumn];
    if (colUpper_)
      upper=colUpper_[iColumn];
    if (basis.getStructStatus(iColumn) == OsiWarmStartBasis::atUpperBound) {
      assert (upper<1.0e20);
      primalSolution_[1][iColumn]=upper;
    } else if (basis.getStructStatus(iColumn) == 
	       OsiWarmStartBasis::atLowerBound) {
      assert (lower>-1.0e20);
      primalSolution_[1][iColumn]=lower;
    } else if (basis.getStructStatus(iColumn) == OsiWarmStartBasis::isFree) {
      if (upper>1.0e20&&lower<-1.0e20) {
	numberFree++;
	primalSolution_[1][iColumn]=0.0;
      } else {
	numberSuperBasic++;
	if (lower>0.0) {
	  primalSolution_[1][iColumn]=lower;
	} else if (upper<0.0) {
	  primalSolution_[1][iColumn]=upper;
	} else {
	  primalSolution_[1][iColumn]=0.0;
	}
      }
    }
  }
  int status = factorize(factorization,basis);
  if (!status) {
    if(numberFree)
      status=1;
    if(numberSuperBasic)
      status=2;
  }
  return status;
}
//#############################################################################
//#############################################################################

OsiSolution::OsiSolution () :

  OsiSimplexModel(),
  columnPrimalInfeasibility_(0.0),
  columnPrimalSequence_(-2),
  rowPrimalInfeasibility_(0.0),
  rowPrimalSequence_(-2), 
  columnDualInfeasibility_(0.0),
  columnDualSequence_(-2),
  rowDualInfeasibility_(0.0),
  rowDualSequence_(-2),
  primalToleranceToGetOptimal_(-1.0),
  dualInfeasibilityWithoutFree_(0.0),
  largeValue_(1.0e15),
  largestPrimalError_(0.0),
  largestDualError_(0.0),
  largestSolutionError_(0.0),
  gotSolution_(false)
{
}


//-----------------------------------------------------------------------------

OsiSolution::~OsiSolution ()
{
}
//#############################################################################
void OsiSolution::setLargeValue( double value) 
{
  if (value>0.0&&value<DBL_MAX)
    largeValue_=value;
}
void
OsiSolution::gutsOfSolution ( const double * rowActivities,
			      const double * columnActivities,
			      const int numberRows,
			      const int numberColumns, 
			      const int* columnStart, 
			      const int * columnLength,
			      const int* row,
			      const double* element,
			      const double* colLower, 
			      const double* colUpper,   
			      const double* objective,
			      const double* rowLower, const double* rowUpper)
{

  gotSolution_=true;
  //work space
  OsiIndexedVector workSpace;
  workSpace.reserve(numberRows_);
  double * base = new double [numberRows_];
  double * save = new double [numberRows_];
  double * previous = new double [numberRows_];

  // accumulate non basic stuff 
  int iRow,iColumn;

  for (iRow=0;iRow<numberRows_;iRow++) {
    if (basis_->getArtifStatus(iRow) != OsiWarmStartBasis::basic) {
      base[iRow]=rowActivities[iRow];
    } else {
      base[iRow]=0.0;
    }
  }
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    if (basis_->getStructStatus(iColumn) != OsiWarmStartBasis::basic) {
      int j;
      double value = columnActivities[iColumn];
      if (value) {
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  iRow=row[j];
	  base[iRow] -= value*element[j];
	}
      }
    } 
  }
  CoinDisjointCopyN ( base, numberRows_ , save);

  double lastError=DBL_MAX;
  int iRefine,numberRefinements=0;
  double * work = workSpace.denseVector();
  factorization_->updateColumn(&workSpace,base);
  for (iRefine=0;iRefine<numberRefinements+1;iRefine++) {
    // check Ax == b  (for basic part)
    for (iRow=0;iRow<numberRows_;iRow++) {
      int iPivot=pivotVariable_[iRow];
      if (iPivot>=numberColumns_) {
	// slack
	work[iPivot-numberColumns_] += base[iRow];
      } else {
	// column
	int j;
	double value = base[iRow];
	for (j=columnStart[iPivot];
	     j<columnStart[iPivot]+columnLength[iPivot];j++) {
	  int jRow=row[j];
	  work[jRow] += value*element[j];
	}
      }
    }
    largestPrimalError_=0.0;
    for (iRow=0;iRow<numberRows_;iRow++) {
      if (fabs(work[iRow]-save[iRow])>largestPrimalError_) {
	largestPrimalError_=fabs(work[iRow]-save[iRow]);
      }
      work[iRow] -= save[iRow];
    }
    if (largestPrimalError_>=lastError) {
      // restore
      double * temp = base;
      base = previous;
      previous=temp;
      break;
    }
    if (iRefine<numberRefinements) {
      // try and make better
      // save this
      double * temp = base;
      base = previous;
      previous=temp;
      double multiplier = 131072.0;
      for (iRow=0;iRow<numberRows_;iRow++) {
	base[iRow] = multiplier*work[iRow];
	work[iRow]=0.0;
      }
      lastError=largestPrimalError_;
      factorization_->updateColumn(&workSpace,base);
      multiplier = 1.0/multiplier;
      for (iRow=0;iRow<numberRows_;iRow++) {
	base[iRow] = previous[iRow] -multiplier*base[iRow];
	work[iRow]=0.0;
      }
    }
  }
  CoinFillN(work,numberRows_,0.0);
  // put solution in correct place
  for (iRow=0;iRow<numberRows_;iRow++) {
    int iPivot=pivotVariable_[iRow];
    if (iPivot>=numberColumns_) {
      // slack
      primalSolution_[0][iPivot-numberColumns_] = base[iRow];
    } else {
      // column
      primalSolution_[1][iPivot] = base[iRow];
     }
  }
  // now look at primal solution
  columnPrimalInfeasibility_=0.0;
  columnPrimalSequence_=-1;
  rowPrimalInfeasibility_=0.0;
  rowPrimalSequence_=-1;
  largestSolutionError_=0.0;
  double * solution = primalSolution_[0];
  for (iRow=0;iRow<numberRows_;iRow++) {
    double infeasibility=0.0;
    if (solution[iRow]>rowUpper[iRow]) {
      infeasibility=solution[iRow]-rowUpper[iRow];
    } else if (solution[iRow]<rowLower[iRow]) {
      infeasibility=rowLower[iRow]-solution[iRow];
    }
    if (infeasibility>rowPrimalInfeasibility_) {
      rowPrimalInfeasibility_=infeasibility;
      rowPrimalSequence_=iRow;
    }
    infeasibility = fabs(rowActivities[iRow]-solution[iRow]);
    if (infeasibility>largestSolutionError_)
      largestSolutionError_=infeasibility;
  }
  solution = primalSolution_[1];
  for (iColumn=0;iColumn<numberColumns_;iColumn++) {
    double infeasibility=0.0;
    if (solution[iColumn]>colUpper[iColumn]) {
      infeasibility=solution[iColumn]-colUpper[iColumn];
    } else if (solution[iColumn]<colLower[iColumn]) {
      infeasibility=colLower[iColumn]-solution[iColumn];
    }
    if (infeasibility>columnPrimalInfeasibility_) {
      columnPrimalInfeasibility_=infeasibility;
      columnPrimalSequence_=iColumn;
    }
    infeasibility = fabs(columnActivities[iColumn]-solution[iColumn]);
    if (infeasibility>largestSolutionError_)
      largestSolutionError_=infeasibility;
  }
  delete [] previous;
  delete [] base;
  delete [] save;
}
int OsiSolution::getSolution ( const double * rowActivities,
			       const double * columnActivities)
{
  if (!factorization_->status()) {
    // get rim stuff
    double * objective = NULL;
    double * colLower = NULL;
    double * colUpper = NULL;
    double * rowLower = NULL;
    double * rowUpper = NULL;
    
    int iRow,iColumn;
    
    if (!objective_) {
      objective = new double[numberColumns_];
      CoinFillN(objective,numberColumns_,0.0);
    }
    if (!colLower_) {
      colLower = new double[numberColumns_];
      CoinFillN(colLower,numberColumns_,0.0);
    }
    if (!colUpper_) {
      colUpper = new double[numberColumns_];
      CoinFillN(colUpper,numberColumns_,DBL_MAX);
    }
    if (!description_.rhsSense) {
      if (!rowLower_) {
	rowLower = new double[numberRows_];
	CoinFillN(rowLower,numberRows_,-DBL_MAX);
      }
      if (!rowUpper_) {
	rowUpper = new double[numberRows_];
	CoinFillN(rowUpper,numberRows_,DBL_MAX);
      }
    } else {
      rowLower = new double[numberRows_];
      rowUpper = new double[numberRows_];
      for (iRow=0;iRow<numberRows_;iRow++) {
	getRowBound(iRow, rowLower[iRow], rowUpper[iRow]);
      }
    }

    if (!description_.matrixType) {
      const int * row = matrix_->getIndices();
      const int * columnStart = matrix_->getVectorStarts();
      const int * columnLength = matrix_->getVectorLengths(); 
      const double * element = matrix_->getElements();
      gutsOfSolution ( rowActivities,
		       columnActivities,
		       numberRows_,
		       numberColumns_, 
		       columnStart, 
		       columnLength,
		       row,
		       element,
		       (colLower) ? colLower : colLower_, 
		       (colUpper) ? colUpper : colUpper_, 
		       (objective) ? objective : objective_, 
		       (rowLower) ? rowLower : rowLower_, 
		       (rowUpper) ? rowUpper : rowUpper_);
    } else {
      if (description_.matrixType==2) {
	gutsOfSolution ( rowActivities,
			 columnActivities,
			 numberRows_,
			 numberColumns_, 
			 columnStart_, 
			 columnLength_,
			 row_,
			 element_,
			 (colLower) ? colLower : colLower_, 
			 (colUpper) ? colUpper : colUpper_, 
			 (objective) ? objective : objective_, 
			 (rowLower) ? rowLower : rowLower_, 
			 (rowUpper) ? rowUpper : rowUpper_);
      } else {
	int * length = new int [numberColumns_];
	for (iColumn=0;iColumn<numberColumns_;iColumn++) {
	  length[iColumn] = columnStart_[iColumn+1]-columnStart_[iColumn];
	}
	gutsOfSolution ( rowActivities,
			 columnActivities,
			 numberRows_,
			 numberColumns_, 
			 columnStart_, 
			 columnLength_,
			 row_,
			 element_,
			 (colLower) ? colLower : colLower_, 
			 (colUpper) ? colUpper : colUpper_, 
			 (objective) ? objective : objective_, 
			 (rowLower) ? rowLower : rowLower_, 
			 (rowUpper) ? rowUpper : rowUpper_);
	delete [] length;
      }
    }
    
    delete [] objective;
    delete [] colLower;
    delete [] colUpper;
    delete [] rowLower;
    delete [] rowUpper;
  }
  return factorization_->status();
}
int OsiSolution::getSolution ( )
{
  double * rowActivities = new double[numberRows_];
  double * columnActivities = new double[numberColumns_];
  CoinDisjointCopyN ( primalSolution_[0], numberRows_ , rowActivities);
  CoinDisjointCopyN ( primalSolution_[1], numberColumns_ , columnActivities);
  int status = getSolution( rowActivities, columnActivities);
  delete [] rowActivities;
  delete [] columnActivities;
  return status;
}
