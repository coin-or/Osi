#include <cassert>

#include "OsiFactorization.hpp"
#include "OsiIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiPackedMatrix.hpp"

#include <stdio.h>
//:class OsiFactorization.  Deals with Factorization and Updates
//  OsiFactorization.  Constructor
OsiFactorization::OsiFactorization (  )
{
  gutsOfInitialize(7);
}

/// Copy constructor 
OsiFactorization::OsiFactorization ( const OsiFactorization &other)
{
  gutsOfInitialize(3);
  gutsOfCopy(other);
}
/// The real work of constructors etc
void OsiFactorization::gutsOfDestructor()
{
  delete [] elementU_ ;
  delete [] startRowU_ ;
  delete [] convertRowToColumnU_ ;
  delete [] indexRowU_ ;
  delete [] indexColumnU_ ;
  delete [] startColumnU_ ;
  delete [] elementL_ ;
  delete [] indexRowL_ ;
  delete [] startColumnL_ ;
  delete [] startColumnR_ ;
  delete [] numberInRow_ ;
  delete [] numberInColumn_ ;
  delete [] numberInColumnPlus_ ;
  delete [] pivotColumn_ ;
  delete [] pivotColumnBack_ ;
  delete [] firstCount_ ;
  delete [] nextCount_ ;
  delete [] lastCount_ ;
  delete [] permute_ ;
  delete [] permuteBack_ ;
  delete [] nextColumn_ ;
  delete [] lastColumn_ ;
  delete [] nextRow_ ;
  delete [] lastRow_ ;
  delete [] saveColumn_ ;
  delete [] markRow_ ;
  delete [] pivotRowL_ ;
  delete [] pivotRegion_ ;
  delete [] elementByRowL_ ;
  delete [] startRowL_ ;
  delete [] indexColumnL_ ;
  delete []sparse_;
#ifdef ZEROFAULT
  // make more visible
  elementU_ = NULL;
  startRowU_ = NULL;
  convertRowToColumnU_ = NULL;
  indexRowU_ = NULL;
  indexColumnU_ = NULL;
  startColumnU_ = NULL;
  elementL_ = NULL;
  indexRowL_ = NULL;
  startColumnL_ = NULL;
  startColumnR_ = NULL;
  numberInRow_ = NULL;
  numberInColumn_ = NULL;
  numberInColumnPlus_ = NULL;
  pivotColumn_ = NULL;
  pivotColumnBack_ = NULL;
  firstCount_ = NULL;
  nextCount_ = NULL;
  lastCount_ = NULL;
  permute_ = NULL;
  permuteBack_ = NULL;
  nextColumn_ = NULL;
  lastColumn_ = NULL;
  nextRow_ = NULL;
  lastRow_ = NULL;
  saveColumn_ = NULL;
  markRow_ = NULL;
  pivotRowL_ = NULL;
  pivotRegion_ = NULL;
  elementByRowL_ = NULL;
  startRowL_ = NULL;
  indexColumnL_ = NULL;
  sparse_= NULL;
#endif
}
// type - 1 bit tolerances etc, 2 rest
void OsiFactorization::gutsOfInitialize(int type)
{
  // change defines to allow for 64 bits
  assert (OSIFACTORIZATION_BITS_PER_INT==8*sizeof(int));
  if ((type&1)!=0) {
    areaFactor_ = 0.0;
    increasingRows_ = 1;
    pivotTolerance_ = 1.0e-1;
    zeroTolerance_ = 1.0e-13;
    slackValue_ = 1.0;
    messageLevel_=0;
    maximumPivots_=200;
    numberTrials_ = 1;
  }
  if ((type&2)!=0) {
    numberCompressions_ = 0;
    biggerDimension_ = 0;
    numberRows_ = 0;
    numberRowsExtra_ = 0;
    maximumRowsExtra_ = 0;
    numberColumns_ = 0;
    numberColumnsExtra_ = 0;
    maximumColumnsExtra_ = 0;
    numberGoodU_ = 0;
    numberGoodL_ = 0;
    totalElements_ = 0;
    factorElements_ = 0;
    status_ = 0;
    doForrestTomlin_=true;
    numberPivots_ = 0;
    pivotColumn_ = NULL;
    permute_ = NULL;
    permuteBack_ = NULL;
    pivotColumnBack_ = NULL;
    startRowU_ = NULL;
    numberInRow_ = NULL;
    numberInColumn_ = NULL;
    numberInColumnPlus_ = NULL;
    firstCount_ = NULL;
    nextCount_ = NULL;
    lastCount_ = NULL;
    nextColumn_ = NULL;
    lastColumn_ = NULL;
    nextRow_ = NULL;
    lastRow_ = NULL;
    saveColumn_ = NULL;
    markRow_ = NULL;
    indexColumnU_ = NULL;
    pivotRowL_ = NULL;
    pivotRegion_ = NULL;
    numberSlacks_ = 0;
    numberU_ = 0;
    lengthU_ = 0;
    lengthAreaU_ = 0;
    elementU_ = NULL;
    indexRowU_ = NULL;
    startColumnU_ = NULL;
    convertRowToColumnU_ = NULL;
    numberL_ = 0;
    baseL_ = 0;
    lengthL_ = 0;
    lengthAreaL_ = 0;
    elementL_ = NULL;
    indexRowL_ = NULL;
    startColumnL_ = NULL;
    numberR_ = 0;
    lengthR_ = 0;
    lengthAreaR_ = 0;
    elementR_ = NULL;
    indexRowR_ = NULL;
    startColumnR_ = NULL;
    elementByRowL_=NULL;
    startRowL_=NULL;
    indexColumnL_=NULL;
    // always switch off sparse
    sparseThreshold_=0;
    sparse_=NULL;
  }
  if ((type&4)!=0) {
    // we need to get 1 element arrays for any with length n+1 !!
    startColumnL_ = new OsiBigIndex [ 1 ];
    startColumnR_ = new OsiBigIndex [ 1 ];
    startRowU_ = new OsiBigIndex [ 1 ];
    numberInRow_ = new int [ 1 ];
    nextRow_ = new int [ 1 ];
    lastRow_ = new int [ 1 ];
    pivotRegion_ = new double [ 1 ];
    permuteBack_ = new int [ 1 ];
    permute_ = new int [ 1 ];
    pivotColumnBack_ = new int [ 1 ];
    startColumnU_ = new OsiBigIndex [ 1 ];
    numberInColumn_ = new int [ 1 ];
    numberInColumnPlus_ = new int [ 1 ];
    pivotColumn_ = new int [ 1 ];
    nextColumn_ = new int [ 1 ];
    lastColumn_ = new int [ 1 ];
#ifdef ZEROFAULT
    startColumnL_[0] = 0;
    startColumnR_[0] = 0;
    startRowU_[0] = 0;
    numberInRow_[0] = 0;
    nextRow_[0] = 0;
    lastRow_[0] = 0;
    pivotRegion_[0] = 0.0;
    permuteBack_[0] = 0;
    permute_[0] = 0;
    pivotColumnBack_[0] = 0;
    startColumnU_[0] = 0;
    numberInColumn_[0] = 0;
    numberInColumnPlus_[0] = 0;
    pivotColumn_[0] = 0;
    nextColumn_[0] = 0;
    lastColumn_[0] = 0;
#endif
  }
}
//Part of LP
int OsiFactorization::factorize (
				 const OsiPackedMatrix & matrix,
				 int rowIsBasic[],
				 int columnIsBasic[],
				 double areaFactor )
{
  // maybe for speed will be better to leave as many regions as possible
  gutsOfDestructor();
  gutsOfInitialize(2);
  if (areaFactor)
    areaFactor_ = areaFactor;
  const int * row = matrix.getIndices();
  const int * columnStart = matrix.getVectorStarts();
  const int * columnLength = matrix.getVectorLengths(); 
  const double * element = matrix.getElements();
  int numberRows=matrix.getNumRows();
  int numberColumns=matrix.getNumCols();
  int numberBasic = 0;
  OsiBigIndex numberElements=0;
  int numberRowBasic=0;

  // compute how much in basis

  int i;

  for (i=0;i<numberRows;i++) {
    if (rowIsBasic[i]>=0)
      numberRowBasic++;
  }

  numberBasic = numberRowBasic;

  for (i=0;i<numberColumns;i++) {
    if (columnIsBasic[i]>=0) {
      numberBasic++;
      numberElements += columnLength[i];
    }
  }
  if ( numberBasic > numberRows ) {
    return -2; // say too many in basis
  }				
  numberElements = 3 * numberBasic + 3 * numberElements + 10000;
  getAreas ( numberRows, numberBasic, numberElements,
	     2 * numberElements );
  //fill
  //copy
  numberBasic=0;
  numberElements=0;
  for (i=0;i<numberRows;i++) {
    if (rowIsBasic[i]>=0) {
      indexRowU_[numberElements]=i;
      indexColumnU_[numberElements]=numberBasic;
      elementU_[numberElements++]=slackValue_;
      numberBasic++;
    }
  }
  for (i=0;i<numberColumns;i++) {
    if (columnIsBasic[i]>=0) {
      int j;
      for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	indexRowU_[numberElements]=row[j];
	indexColumnU_[numberElements]=numberBasic;
	elementU_[numberElements++]=element[j];
      }
      numberBasic++;
    }
  }
  lengthU_ = numberElements;

  preProcess ( 0 );
  factor (  );
  numberBasic=0;
  if (status_ == 0) {
    int * permuteBack = permuteBack_;
    int * back = pivotColumnBack_;
    for (i=0;i<numberRows;i++) {
      if (rowIsBasic[i]>=0) {
	rowIsBasic[i]=permuteBack[back[numberBasic++]];
      }
    }
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	columnIsBasic[i]=permuteBack[back[numberBasic++]];
      }
    }
    if (increasingRows_>1) {
      // Set up permutation vector
      if (increasingRows_<3) {
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	CoinDisjointCopyN ( permute_, numberRows_ , pivotColumn_  );
	CoinDisjointCopyN ( permuteBack_, numberRows_ , pivotColumnBack_  );
      }
    } else {
      // Set up permutation vector
      // (we could use permute_ instead of pivotColumn (not back though))
      for (i=0;i<numberRows_;i++) {
	int k=pivotColumn_[i];
	pivotColumn_[i]=pivotColumnBack_[i];
	pivotColumnBack_[i]=k;
      }
    }
  } else if (status_ == -1) {
    // mark as basic or non basic
    for (i=0;i<numberRows_;i++) {
      if (rowIsBasic[i]>=0) {
	if (pivotColumn_[numberBasic]>=0) 
	  rowIsBasic[i]=pivotColumn_[numberBasic];
	else 
	  rowIsBasic[i]=-1;
	numberBasic++;
      }
    }
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	if (pivotColumn_[numberBasic]>=0) 
	  columnIsBasic[i]=pivotColumn_[numberBasic];
	 else 
	  columnIsBasic[i]=-1;
	numberBasic++;
      }
    }
  }

  return status_;
}
//Part of LP
int OsiFactorization::factorize ( int numberRows, int numberColumns,
				  const OsiBigIndex * columnStart,
				  const int * columnLength,
				  const int * row,
				  const double * element,
				  int rowIsBasic[],
				  int columnIsBasic[],
				  double areaFactor ,
				  double * rowScale,
				  double * columnScale)
{
  // maybe for speed will be better to leave as many regions as possible
  gutsOfDestructor();
  gutsOfInitialize(2);
  if (areaFactor)
    areaFactor_ = areaFactor;
  int numberBasic = 0;
  OsiBigIndex numberElements=0;
  int numberRowBasic=0;

  // compute how much in basis

  int i;

  for (i=0;i<numberRows;i++) {
    if (rowIsBasic[i]>=0)
      numberRowBasic++;
  }

  numberBasic = numberRowBasic;

  for (i=0;i<numberColumns;i++) {
    if (columnIsBasic[i]>=0) {
      numberBasic++;
      numberElements += columnLength[i];
    }
  }
  if ( numberBasic > numberRows ) {
    return -2; // say too many in basis
  }				
  numberElements = 3 * numberBasic + 3 * numberElements + 10000;
  getAreas ( numberRows, numberBasic, numberElements,
	     2 * numberElements );
  //fill
  //copy
  numberBasic=0;
  numberElements=0;
  for (i=0;i<numberRows;i++) {
    if (rowIsBasic[i]>=0) {
      indexRowU_[numberElements]=i;
      indexColumnU_[numberElements]=numberBasic;
      elementU_[numberElements++]=slackValue_;
      numberBasic++;
    }
  }
  if (!rowScale) {
    // no scaling
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	int j;
	for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	  indexRowU_[numberElements]=row[j];
	  indexColumnU_[numberElements]=numberBasic;
	  elementU_[numberElements++]=element[j];
	}
	numberBasic++;
      }
    }
  } else {
    // scaling
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	int j;
	double scale = columnScale[i];
	for (j=columnStart[i];j<columnStart[i]+columnLength[i];j++) {
	  int iRow = row[j];
	  indexRowU_[numberElements]=iRow;
	  indexColumnU_[numberElements]=numberBasic;
	  elementU_[numberElements++]=element[j]*scale*rowScale[iRow];
	}
	numberBasic++;
      }
    }
  }
  lengthU_ = numberElements;

  preProcess ( 0 );
  factor (  );
  numberBasic=0;
  if (status_ == 0) {
    int * permuteBack = permuteBack_;
    int * back = pivotColumnBack_;
    for (i=0;i<numberRows;i++) {
      if (rowIsBasic[i]>=0) {
	rowIsBasic[i]=permuteBack[back[numberBasic++]];
      }
    }
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	columnIsBasic[i]=permuteBack[back[numberBasic++]];
      }
    }
    if (increasingRows_>1) {
      // Set up permutation vector
      if (increasingRows_<3) {
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	CoinDisjointCopyN ( permute_, numberRows_ , pivotColumn_  );
	CoinDisjointCopyN ( permuteBack_, numberRows_ , pivotColumnBack_  );
      }
    } else {
      // Set up permutation vector
      // (we could use permute_ instead of pivotColumn (not back though))
      for (i=0;i<numberRows_;i++) {
	int k=pivotColumn_[i];
	pivotColumn_[i]=pivotColumnBack_[i];
	pivotColumnBack_[i]=k;
      }
    }
  } else if (status_ == -1) {
    // mark as basic or non basic
    for (i=0;i<numberRows_;i++) {
      if (rowIsBasic[i]>=0) {
	if (pivotColumn_[numberBasic]>=0) 
	  rowIsBasic[i]=pivotColumn_[numberBasic];
	numberBasic++;
      }
    }
    for (i=0;i<numberColumns;i++) {
      if (columnIsBasic[i]>=0) {
	if (pivotColumn_[numberBasic]>=0) 
	  columnIsBasic[i]=pivotColumn_[numberBasic];
	numberBasic++;
      }
    }
  }

  return status_;
}
//Given as triplets
int OsiFactorization::factorize (
			     int numberOfRows,
			     int numberOfColumns,
			     OsiBigIndex numberOfElements,
			     OsiBigIndex maximumL,
			     OsiBigIndex maximumU,
			     const int indicesRow[],
			     const int indicesColumn[],
			     const double elements[] ,
			     int permutation[],
			     double areaFactor)
{
  gutsOfDestructor();
  gutsOfInitialize(2);
  if (areaFactor)
    areaFactor_ = areaFactor;
  getAreas ( numberOfRows, numberOfColumns, maximumL, maximumU );
  //copy
  CoinDisjointCopyN ( indicesRow, numberOfElements, indexRowU_ );
  CoinDisjointCopyN ( indicesColumn, numberOfElements, indexColumnU_ );
  CoinDisjointCopyN ( elements, numberOfElements, elementU_ );
  lengthU_ = numberOfElements;
  preProcess ( 0 );
  factor (  );
  //say which column is pivoting on which row
  int i;
  if (status_ == 0) {
    if (increasingRows_>1) {
      int * permuteBack = permuteBack_;
      int * back = pivotColumnBack_;
      // permute so slacks on own rows etc
      for (i=0;i<numberOfColumns;i++) {
	permutation[i]=permuteBack[back[i]];
      }
      // Set up permutation vector
      if (increasingRows_<3) {
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	CoinDisjointCopyN ( permute_, numberRows_ , pivotColumn_  );
	CoinDisjointCopyN ( permuteBack_, numberRows_ , pivotColumnBack_  );
      }
    } else {
      // Set up permutation vector
      // (we could use permute_ instead of pivotColumn (not back though))
      for (i=0;i<numberOfRows;i++) {
	int k=pivotColumn_[i];
	pivotColumn_[i]=pivotColumnBack_[i];
	pivotColumnBack_[i]=k;
      }
      for (i=0;i<numberOfColumns;i++) {
	permutation[i]=i;
      }
    }
  } else if (status_ == -1) {
    // mark as basic or non basic
    for (i=0;i<numberOfColumns;i++) {
      if (pivotColumn_[i]>=0) {
	permutation[i]=pivotColumn_[i];
      } else {
	permutation[i]=-1;
      }
    }
  }

  return status_;
}
/* Two part version for flexibility
   This part creates arrays for user to fill.
   maximumL is guessed maximum size of L part of
   final factorization, maximumU of U part.  These are multiplied by
   areaFactor which can be computed by user or internally.  
   returns 0 -okay, -99 memory */
int 
OsiFactorization::factorizePart1 ( int numberOfRows,
				   int numberOfColumns,
				   OsiBigIndex numberOfElements,
				   int * indicesRow[],
				   int * indicesColumn[],
				   double * elements[],
				   double areaFactor)
{
  // maybe for speed will be better to leave as many regions as possible
  gutsOfDestructor();
  gutsOfInitialize(2);
  if (areaFactor)
    areaFactor_ = areaFactor;
  OsiBigIndex numberElements = 3 * numberOfRows + 3 * numberOfElements + 10000;
  getAreas ( numberOfRows, numberOfRows, numberElements,
	     2 * numberElements );
  // need to trap memory for -99 code
  *indicesRow = indexRowU_ ;
  *indicesColumn = indexColumnU_ ;
  *elements = elementU_ ;
  lengthU_ = numberOfElements;
  return 0;
}
/* This is part two of factorization
   Arrays belong to factorization and were returned by part 1
   If status okay, permutation has pivot rows.
   If status is singular, then basic variables have +1 and ones thrown out have -INT_MAX
   to say thrown out.
   returns 0 -okay, -1 singular, -99 memory */
int 
OsiFactorization::factorizePart2 (int permutation[],int exactNumberElements)
{
  lengthU_ = exactNumberElements;
  preProcess ( 0 );
  factor (  );
  //say which column is pivoting on which row
  int i;
  int * permuteBack = permuteBack_;
  int * back = pivotColumnBack_;
  // permute so slacks on own rows etc
  for (i=0;i<numberColumns_;i++) {
    permutation[i]=permuteBack[back[i]];
  }
  if (status_ == 0) {
    if (increasingRows_>1) {
      // Set up permutation vector
      if (increasingRows_<3) {
	// these arrays start off as copies of permute
	// (and we could use permute_ instead of pivotColumn (not back though))
	CoinDisjointCopyN ( permute_, numberRows_ , pivotColumn_  );
	CoinDisjointCopyN ( permuteBack_, numberRows_ , pivotColumnBack_  );
      }
    } else {
      // Set up permutation vector
      // (we could use permute_ instead of pivotColumn (not back though))
      for (i=0;i<numberRows_;i++) {
	int k=pivotColumn_[i];
	pivotColumn_[i]=pivotColumnBack_[i];
	pivotColumnBack_[i]=k;
      }
    }
  } else if (status_ == -1) {
    // mark as basic or non basic
    for (i=0;i<numberColumns_;i++) {
      if (pivotColumn_[i]>=0) {
	permutation[i]=pivotColumn_[i];
      } else {
	permutation[i]=-1;
      }
    }
  }

  return status_;
}

//  ~OsiFactorization.  Destructor
OsiFactorization::~OsiFactorization (  )
{
  gutsOfDestructor();
}

//  show_self.  Debug show object
void
OsiFactorization::show_self (  ) const
{
  int i;

  for ( i = 0; i < numberRows_; i++ ) {
    std::cout << "r " << i << " " << pivotColumn_[i]
      << " " << pivotColumnBack_[i]
      << " " << permute_[i]
      << " " << permuteBack_[i]
      << " " << pivotRowL_[i]
      << " " << pivotColumn_[i]
      << " " << pivotRegion_[i] << std::endl;
  }				
  for ( i = 0; i < numberRows_; i++ ) {
    std::cout << "u " << i << " " << numberInColumn_[i] << std::endl;
    int j;
    for ( j = startColumnU_[i]; j < startColumnU_[i] + numberInColumn_[i];
	  j++ ) {
      std::cout << indexRowU_[j] << " " << elementU_[j] << std::endl;
    }				
  }				
  for ( i = 0; i < numberRows_; i++ ) {
    std::cout << "l " << i << " " << startColumnL_[i + 1] -
      startColumnL_[i] << std::endl;
    int j;
    for ( j = startColumnL_[i]; j < startColumnL_[i + 1]; j++ ) {
      std::cout << indexRowL_[j] << " " << elementL_[j] << std::endl;
    }				
  }				

}

//  getAreas.  Gets space for a factorization
//called by constructors
void
OsiFactorization::getAreas ( int numberOfRows,
			 int numberOfColumns,
			 OsiBigIndex maximumL,
			 OsiBigIndex maximumU )
{
  int extraSpace = maximumPivots_;

  numberRows_ = numberOfRows;
  numberColumns_ = numberOfColumns;
  maximumRowsExtra_ = numberRows_ + extraSpace;
  numberRowsExtra_ = numberRows_;
  maximumColumnsExtra_ = numberColumns_ + extraSpace;
  numberColumnsExtra_ = numberColumns_;
  lengthAreaU_ = maximumU;
  lengthAreaL_ = maximumL;
  if ( !areaFactor_ ) {
    areaFactor_ = 1.0;
  }
  if ( areaFactor_ != 1.0 && (messageLevel_&1)) {
    std::cout<<"Increasing factorization areas by "<<areaFactor_<<std::endl;
    lengthAreaU_ *= areaFactor_;
    lengthAreaL_ *= areaFactor_;
  }				
  elementU_ = new double [ lengthAreaU_ ];
  indexRowU_ = new int [ lengthAreaU_ ];
  indexColumnU_ = new int [ lengthAreaU_ ];
  elementL_ = new double [ lengthAreaL_ ];
  indexRowL_ = new int [ lengthAreaL_ ];
  startColumnL_ = new OsiBigIndex [ numberRows_ + 1 ];
  startColumnL_[0] = 0;
  startRowU_ = new OsiBigIndex [ maximumRowsExtra_ + 1 ];
  // make sure this is valid
  startRowU_[maximumRowsExtra_]=0;
  numberInRow_ = new int [ maximumRowsExtra_ + 1 ];
  markRow_ = new int [ numberRows_ ];
  pivotRowL_ = new int [ numberRows_ + 1 ];
  nextRow_ = new int [ maximumRowsExtra_ + 1 ];
  lastRow_ = new int [ maximumRowsExtra_ + 1 ];
  permute_ = new int [ maximumRowsExtra_ + 1 ];
  pivotRegion_ = new double [ maximumRowsExtra_ + 1 ];
#ifdef ZEROFAULT
  memset(elementU_,'a',lengthAreaU_*sizeof(double));
  memset(indexRowU_,'b',lengthAreaU_*sizeof(int));
  memset(indexColumnU_,'c',lengthAreaU_*sizeof(int));
  memset(elementL_,'d',lengthAreaL_*sizeof(double));
  memset(indexRowL_,'e',lengthAreaL_*sizeof(int));
  memset(startColumnL_+1,'f',numberRows_*sizeof(OsiBigIndex));
  memset(startRowU_,'g',maximumRowsExtra_*sizeof(OsiBigIndex));
  memset(numberInRow_,'h',(maximumRowsExtra_+1)*sizeof(int));
  memset(markRow_,'i',numberRows_*sizeof(int));
  memset(pivotRowL_,'j',(numberRows_+1)*sizeof(int));
  memset(nextRow_,'k',(maximumRowsExtra_+1)*sizeof(int));
  memset(lastRow_,'l',(maximumRowsExtra_+1)*sizeof(int));
  memset(permute_,'l',(maximumRowsExtra_+1)*sizeof(int));
  memset(pivotRegion_,'m',(maximumRowsExtra_+1)*sizeof(double));
#endif
  startColumnU_ = new OsiBigIndex [ maximumColumnsExtra_ + 1 ];
  numberInColumn_ = new int [ maximumColumnsExtra_ + 1 ];
  numberInColumnPlus_ = new int [ maximumColumnsExtra_ + 1 ];
  pivotColumn_ = new int [ maximumColumnsExtra_ + 1 ];
  nextColumn_ = new int [ maximumColumnsExtra_ + 1 ];
  lastColumn_ = new int [ maximumColumnsExtra_ + 1 ];
  saveColumn_ = new int [ numberColumns_ ];
#ifdef ZEROFAULT
  memset(startColumnU_,'a',(maximumColumnsExtra_+1)*sizeof(OsiBigIndex));
  memset(numberInColumn_,'b',(maximumColumnsExtra_+1)*sizeof(int));
  memset(numberInColumnPlus_,'c',(maximumColumnsExtra_+1)*sizeof(int));
  memset(pivotColumn_,'d',(maximumColumnsExtra_+1)*sizeof(int));
  memset(nextColumn_,'e',(maximumColumnsExtra_+1)*sizeof(int));
  memset(lastColumn_,'f',(maximumColumnsExtra_+1)*sizeof(int));
  memset(saveColumn_,'g',numberColumns_*sizeof(int));
#endif
  if ( numberRows_ + numberColumns_ ) {
    if ( numberRows_ > numberColumns_ ) {
      biggerDimension_ = numberRows_;
    } else {
      biggerDimension_ = numberColumns_;
    }				
    firstCount_ = new int [ biggerDimension_ + 2 ];
    nextCount_ = new int [ numberRows_ + numberColumns_ ];
    lastCount_ = new int [ numberRows_ + numberColumns_ ];
#ifdef ZEROFAULT
    memset(firstCount_,'g',(biggerDimension_ + 2 )*sizeof(int));
    memset(nextCount_,'h',(numberRows_+numberColumns_)*sizeof(int));
    memset(lastCount_,'i',(numberRows_+numberColumns_)*sizeof(int));
#endif
  } else {
    firstCount_ = new int [ 2 ];
    nextCount_ = NULL;
    lastCount_ = NULL;
#ifdef ZEROFAULT
    memset(firstCount_,'g', 2 *sizeof(int));
#endif
    biggerDimension_ = 0;
  }				
}

//  preProcess.  PreProcesses raw triplet data
//state is 0 - triplets, 1 - some counts etc , 2 - ..
void
OsiFactorization::preProcess ( int state,
			   int possibleDuplicates )
{
  int *indexRow = indexRowU_;
  int *indexColumn = indexColumnU_;
  double *element = elementU_;
  OsiBigIndex numberElements = lengthU_;
  int *numberInRow = numberInRow_;
  int *numberInColumn = numberInColumn_;
  OsiBigIndex *startRow = startRowU_;
  OsiBigIndex *startColumn = startColumnU_;
  int numberRows = numberRows_;
  int numberColumns = numberColumns_;

  totalElements_ = numberElements;
  //state falls through to next state
  switch ( state ) {
  case 0:			//counts
    {
      CoinFillN ( numberInRow, numberRows + 1 , 0);
      CoinFillN ( numberInColumn, maximumColumnsExtra_ + 1 , 0);
      OsiBigIndex i;
      for ( i = 0; i < numberElements; i++ ) {
	int iRow = indexRow[i];
	int iColumn = indexColumn[i];

	numberInRow[iRow]++;
	numberInColumn[iColumn]++;
      }				
    }
  case -1:			//sort
  case 1:			//sort
    {
      OsiBigIndex i, k;

      i = 0;
      int iColumn;
      for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
	//position after end of Column
	i += numberInColumn[iColumn];
	startColumn[iColumn] = i;
      }				
      for ( k = numberElements - 1; k >= 0; k-- ) {
	int iColumn = indexColumn[k];

	if ( iColumn >= 0 ) {
	  double value = element[k];
	  int iRow = indexRow[k];

	  indexColumn[k] = -1;
	  while ( true ) {
	    OsiBigIndex iLook = startColumn[iColumn] - 1;

	    startColumn[iColumn] = iLook;
	    double valueSave = element[iLook];
	    int iColumnSave = indexColumn[iLook];
	    int iRowSave = indexRow[iLook];

	    element[iLook] = value;
	    indexRow[iLook] = iRow;
	    indexColumn[iLook] = -1;
	    if ( iColumnSave >= 0 ) {
	      iColumn = iColumnSave;
	      value = valueSave;
	      iRow = iRowSave;
	    } else {
	      break;
	    }			
	  }			/* endwhile */
	}			
      }				
    }
  case 2:			//move largest in column to beginning
    //and do row part
    {
      OsiBigIndex i, k;

      i = 0;
      int iRow;
      for ( iRow = 0; iRow < numberRows; iRow++ ) {
	startRow[iRow] = i;
	i += numberInRow[iRow];
      }				
      CoinFillN ( numberInRow, numberRows , 0);
      int iColumn;
      for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
	int number = numberInColumn[iColumn];

	if ( number ) {
	  OsiBigIndex first = startColumn[iColumn];
	  OsiBigIndex largest = first;
	  int iRowSave = indexRow[first];
	  double valueSave = element[first];
	  double valueLargest = fabs ( valueSave );
	  int iLook = numberInRow[iRowSave];

	  numberInRow[iRowSave] = iLook + 1;
	  indexColumn[startRow[iRowSave] + iLook] = iColumn;
	  for ( k = first + 1; k < first + number; k++ ) {
	    int iRow = indexRow[k];
	    int iLook = numberInRow[iRow];

	    numberInRow[iRow] = iLook + 1;
	    indexColumn[startRow[iRow] + iLook] = iColumn;
	    double value = element[k];
	    double valueAbs = fabs ( value );

	    if ( valueAbs > valueLargest ) {
	      valueLargest = valueAbs;
	      largest = k;
	    }			
	  }			
	  indexRow[first] = indexRow[largest];
	  element[first] = element[largest];
	  indexRow[largest] = iRowSave;
	  element[largest] = valueSave;
	}			
      }				
    }
  case 3:			//links and initialize pivots
    {
      //set markRow so no rows updated
      CoinFillN ( markRow_, numberRows_, -1 );
      int *lastRow = lastRow_;
      int *nextRow = nextRow_;
      int *lastColumn = lastColumn_;
      int *nextColumn = nextColumn_;

      CoinFillN ( firstCount_, biggerDimension_ + 2, -1 );
      CoinFillN ( pivotColumn_, numberColumns_, -1 );
      CoinFillN ( numberInColumnPlus_, maximumColumnsExtra_ + 1, 0 );
      int iRow;
      for ( iRow = 0; iRow < numberRows; iRow++ ) {
	lastRow[iRow] = iRow - 1;
	nextRow[iRow] = iRow + 1;
	int number = numberInRow[iRow];

	addLink ( iRow, number );
      }				
      lastRow[maximumRowsExtra_] = numberRows - 1;
      nextRow[maximumRowsExtra_] = 0;
      lastRow[0] = maximumRowsExtra_;
      nextRow[numberRows - 1] = maximumRowsExtra_;
      startRow[maximumRowsExtra_] = numberElements;
      int iColumn;
      for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
	lastColumn[iColumn] = iColumn - 1;
	nextColumn[iColumn] = iColumn + 1;
	int number = numberInColumn[iColumn];

	addLink ( iColumn + numberRows, number );
      }				
      lastColumn[maximumColumnsExtra_] = numberColumns - 1;
      nextColumn[maximumColumnsExtra_] = 0;
      lastColumn[0] = maximumColumnsExtra_;
      nextColumn[numberColumns - 1] = maximumColumnsExtra_;
      startColumn[maximumColumnsExtra_] = numberElements;
    }
  }				/* endswitch */
}

//Does most of factorization
int
OsiFactorization::factor (  )
{
  //sparse
  status_ = factorSparse (  );
  switch ( status_ ) {
  case 0:			//finished
    totalElements_ = 0;
    {
      if ( numberGoodU_ < numberRows_ ) {
	int i, k;
	
	int * swap = permute_;
	permute_ = nextRow_;
	nextRow_ = swap;
	for ( i = 0; i < numberRows_; i++ ) {
	  lastRow_[i] = -1;
	}			
	for ( i = 0; i < numberColumns_; i++ ) {
	  lastColumn_[i] = -1;
	}			
	for ( i = 0; i < numberGoodU_; i++ ) {
	  int goodRow = pivotRowL_[i];	//valid pivot row

	  if ( i != nextRow_[goodRow] ) {
	    std::cout << "bad1" << std::endl;
	  }			
	  int goodColumn = pivotColumn_[i];

	  lastRow_[goodRow] = goodColumn;	//will now have -1 or column sequence
	  lastColumn_[goodColumn] = goodRow;	//will now have -1 or row sequence
	}			
	nextRow_ = NULL;
	k = 0;
	//copy back and count
	for ( i = 0; i < numberRows_; i++ ) {
	  permute_[i] = lastRow_[i];
	  if ( permute_[i] < 0 ) {
	    //std::cout << i << " " <<permute_[i] << std::endl;
	  } else {
	    k++;
	  }			
	}			
	for ( i = 0; i < numberColumns_; i++ ) {
	  pivotColumn_[i] = lastColumn_[i];
	}			
	std::cout <<"Factorization has "<<numberRows_-k
		  <<" singularities"<<std::endl;
	status_ = -1;
      }				
    }
    break;
  default:
    //singular ? or some error
    std::cout << "Error " << status_ << std::endl;
    break;
  }				/* endswitch */
  //clean up
  if ( !status_ ) {
    if ( messageLevel_ & 4 )
    std::cout<<"Factorization did "<<numberCompressions_
	     <<" compressions"<<std::endl;
    if ( numberCompressions_ > 10 * numberRows_ ) {
      areaFactor_ *= 1.1;
    }
    cleanup (  );
  }				
  return status_;
}


//  pivotRowSingleton.  Does one pivot on Row Singleton in factorization
bool
OsiFactorization::pivotRowSingleton ( int pivotRow,
				  int pivotColumn )
{
  //store pivot columns (so can easily compress)
  OsiBigIndex startColumn = startColumnU_[pivotColumn];
  int numberDoColumn = numberInColumn_[pivotColumn] - 1;
  OsiBigIndex endColumn = startColumn + numberDoColumn + 1;
  int pivotRowPosition = startColumn;
  int iRow = indexRowU_[pivotRowPosition];

  while ( iRow != pivotRow ) {
    pivotRowPosition++;
    iRow = indexRowU_[pivotRowPosition];
  }				/* endwhile */
#if OSI_DEBUG
  {
    if ( pivotRowPosition >= endColumn ) {
      abort (  );
    }				
  }
#endif
  //store column in L, compress in U and take column out
  OsiBigIndex l = lengthL_;

  if ( l + numberDoColumn > lengthAreaL_ ) {
    //need another area
    std::cout << "code pivot 1" << std::endl;
    //leave gap so starts will work
    abort (  );
  }				
  pivotRowL_[numberGoodL_] = pivotRow;
  startColumnL_[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL_[numberGoodL_] = l + numberDoColumn;
  lengthL_ += numberDoColumn;
  double pivotElement = elementU_[pivotRowPosition];
  double pivotMultiplier = 1.0 / pivotElement;

  pivotRegion_[numberGoodU_] = pivotMultiplier;
  OsiBigIndex i;

  for ( i = startColumn; i < pivotRowPosition; i++ ) {
    int iRow = indexRowU_[i];

    indexRowL_[l] = iRow;
    elementL_[l] = elementU_[i] * pivotMultiplier;
    l++;
    //take out of row list
    OsiBigIndex start = startRowU_[iRow];
    int numberInRow = numberInRow_[iRow];
    OsiBigIndex end = start + numberInRow;
    OsiBigIndex where = start;

    while ( indexColumnU_[where] != pivotColumn ) {
      where++;
    }				/* endwhile */
#if OSI_DEBUG
    if ( where >= end ) {
      abort (  );
    }				
#endif
    indexColumnU_[where] = indexColumnU_[end - 1];
    numberInRow--;
    numberInRow_[iRow] = numberInRow;
    deleteLink ( iRow );
    addLink ( iRow, numberInRow );
  }				
  for ( i = pivotRowPosition + 1; i < endColumn; i++ ) {
    int iRow = indexRowU_[i];

    indexRowL_[l] = iRow;
    elementL_[l] = elementU_[i] * pivotMultiplier;
    l++;
    //take out of row list
    OsiBigIndex start = startRowU_[iRow];
    int numberInRow = numberInRow_[iRow];
    OsiBigIndex end = start + numberInRow;
    OsiBigIndex where = start;

    while ( indexColumnU_[where] != pivotColumn ) {
      where++;
    }				/* endwhile */
#if OSI_DEBUG
    if ( where >= end ) {
      abort (  );
    }				
#endif
    indexColumnU_[where] = indexColumnU_[end - 1];
    numberInRow--;
    numberInRow_[iRow] = numberInRow;
    deleteLink ( iRow );
    addLink ( iRow, numberInRow );
  }				
  numberInColumn_[pivotColumn] = 0;
  //modify linked list for pivots
  numberInRow_[pivotRow] = 0;
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  //take out this bit of indexColumnU
  int next = nextRow_[pivotRow];
  int last = lastRow_[pivotRow];

  nextRow_[last] = next;
  lastRow_[next] = last;
  nextRow_[pivotRow] = numberGoodU_;	//use for permute
  return true;
}

//  pivotColumnSingleton.  Does one pivot on Column Singleton in factorization
bool
OsiFactorization::pivotColumnSingleton ( int pivotRow,
				     int pivotColumn )
{
  //store pivot columns (so can easily compress)
  int numberDoRow = numberInRow_[pivotRow] - 1;
  OsiBigIndex startColumn = startColumnU_[pivotColumn];
  int put = 0;
  OsiBigIndex startRow = startRowU_[pivotRow];
  OsiBigIndex endRow = startRow + numberDoRow + 1;
  OsiBigIndex i;

  for ( i = startRow; i < endRow; i++ ) {
    int iColumn = indexColumnU_[i];

    if ( iColumn != pivotColumn ) {
      saveColumn_[put++] = iColumn;
    }				
  }				
  //take out this bit of indexColumnU
  int next = nextRow_[pivotRow];
  int last = lastRow_[pivotRow];

  nextRow_[last] = next;
  lastRow_[next] = last;
  nextRow_[pivotRow] = numberGoodU_;	//use for permute
  //clean up counts
  double pivotElement = elementU_[startColumn];

  pivotRegion_[numberGoodU_] = 1.0 / pivotElement;
  numberInColumn_[pivotColumn] = 0;
  //numberInColumnPlus_[pivotColumn]++;
  //move pivot row in other columns to safe zone
  for ( i = 0; i < numberDoRow; i++ ) {
    int iColumn = saveColumn_[i];

    if ( numberInColumn_[iColumn] ) {
      int number = numberInColumn_[iColumn] - 1;

      //modify linked list
      deleteLink ( iColumn + numberRows_ );
      addLink ( iColumn + numberRows_, number );
      //move pivot row element
      if ( number ) {
	OsiBigIndex start = startColumnU_[iColumn];
	OsiBigIndex pivot = start;
	int iRow = indexRowU_[pivot];

	while ( iRow != pivotRow ) {
	  pivot++;
	  iRow = indexRowU_[pivot];
	}			/* endwhile */
#if OSI_DEBUG
	{
	  OsiBigIndex end_debug = startColumnU_[iColumn] +

	    numberInColumn_[iColumn];
	  if ( pivot >= end_debug ) {
	    abort (  );
	  }			
	}
#endif
	if ( pivot != start ) {
	  //move largest one up
	  double value = elementU_[start];

	  iRow = indexRowU_[start];
	  elementU_[start] = elementU_[pivot];
	  indexRowU_[start] = indexRowU_[pivot];
	  elementU_[pivot] = elementU_[start + 1];
	  indexRowU_[pivot] = indexRowU_[start + 1];
	  elementU_[start + 1] = value;
	  indexRowU_[start + 1] = iRow;
	} else {
	  //find new largest element
	  int iRowSave = indexRowU_[start + 1];
	  double valueSave = elementU_[start + 1];
	  double valueLargest = fabs ( valueSave );
	  OsiBigIndex end = start + numberInColumn_[iColumn];
	  OsiBigIndex largest = start + 1;

	  OsiBigIndex k;
	  for ( k = start + 2; k < end; k++ ) {
	    double value = elementU_[k];
	    double valueAbs = fabs ( value );

	    if ( valueAbs > valueLargest ) {
	      valueLargest = valueAbs;
	      largest = k;
	    }			
	  }			
	  indexRowU_[start + 1] = indexRowU_[largest];
	  elementU_[start + 1] = elementU_[largest];
	  indexRowU_[largest] = iRowSave;
	  elementU_[largest] = valueSave;
	}			
      }				
      //clean up counts
      numberInColumn_[iColumn]--;
      numberInColumnPlus_[iColumn]++;
      startColumnU_[iColumn]++;
    }				
  }				
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  numberInRow_[pivotRow] = 0;
  //put in dummy pivot in L
  OsiBigIndex l = lengthL_;

  pivotRowL_[numberGoodL_] = pivotRow;
  startColumnL_[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL_[numberGoodL_] = l;
  return true;
}


//  getColumnSpace.  Gets space for one Column with given length
//may have to do compression  (returns true)
//also moves existing vector
bool
OsiFactorization::getColumnSpace ( int iColumn,
			       int extraNeeded )
{
  int number = numberInColumnPlus_[iColumn] +

    numberInColumn_[iColumn];
  OsiBigIndex space = lengthAreaU_ - startColumnU_[maximumColumnsExtra_];

  if ( space < extraNeeded + number + 2 ) {
    //compression
    int iColumn = nextColumn_[maximumColumnsExtra_];
    OsiBigIndex put = 0;

    while ( iColumn != maximumColumnsExtra_ ) {
      //move
      OsiBigIndex get;
      OsiBigIndex getEnd;

      if ( startColumnU_[iColumn] >= 0 ) {
	get = startColumnU_[iColumn]
	  - numberInColumnPlus_[iColumn];
	getEnd = startColumnU_[iColumn] + numberInColumn_[iColumn];
	startColumnU_[iColumn] = put + numberInColumnPlus_[iColumn];
      } else {
	get = -startColumnU_[iColumn];
	getEnd = get + numberInColumn_[iColumn];
	startColumnU_[iColumn] = -put;
      }				
      OsiBigIndex i;
      for ( i = get; i < getEnd; i++ ) {
	indexRowU_[put] = indexRowU_[i];
	elementU_[put] = elementU_[i];
	put++;
      }				
      iColumn = nextColumn_[iColumn];
      numberCompressions_++;
    }				/* endwhile */
    startColumnU_[maximumColumnsExtra_] = put;
    space = lengthAreaU_ - put;
    if ( extraNeeded == INT_MAX >> 1 ) {
      return true;
    }				
    if ( space < extraNeeded + number + 2 ) {
      //need more space
      //if we can allocate bigger then do so and copy
      //if not then return so code can start again
      status_ = -99;
      return false;
    }				
  }				
  OsiBigIndex put = startColumnU_[maximumColumnsExtra_];
  int next = nextColumn_[iColumn];
  int last = lastColumn_[iColumn];

  if ( extraNeeded || next != maximumColumnsExtra_ ) {
    //out
    nextColumn_[last] = next;
    lastColumn_[next] = last;
    //in at end
    last = lastColumn_[maximumColumnsExtra_];
    nextColumn_[last] = iColumn;
    lastColumn_[maximumColumnsExtra_] = iColumn;
    lastColumn_[iColumn] = last;
    nextColumn_[iColumn] = maximumColumnsExtra_;
    //move
    OsiBigIndex get = startColumnU_[iColumn]
      - numberInColumnPlus_[iColumn];

    startColumnU_[iColumn] = put + numberInColumnPlus_[iColumn];
    if ( number < 50 ) {
      int *indexRow = indexRowU_;
      double *element = elementU_;
      int i = 0;

      if ( ( number & 1 ) != 0 ) {
	element[put] = element[get];
	indexRow[put] = indexRow[get];
	i = 1;
      }				
      for ( ; i < number; i += 2 ) {
	double value0 = element[get + i];
	double value1 = element[get + i + 1];
	int index0 = indexRow[get + i];
	int index1 = indexRow[get + i + 1];

	element[put + i] = value0;
	element[put + i + 1] = value1;
	indexRow[put + i] = index0;
	indexRow[put + i + 1] = index1;
      }				
    } else {
      CoinDisjointCopyN ( &indexRowU_[get], number, &indexRowU_[put] );
      CoinDisjointCopyN ( &elementU_[get], number, &elementU_[put] );
    }				
    put += number;
    get += number;
    //add 4 for luck
    startColumnU_[maximumColumnsExtra_] = put + extraNeeded + 4;
  } else {
    //take off space
    startColumnU_[maximumColumnsExtra_] = startColumnU_[last] +
      numberInColumn_[last];
  }				
  return true;
}

//  getRowSpace.  Gets space for one Row with given length
//may have to do compression  (returns true)
//also moves existing vector
bool
OsiFactorization::getRowSpace ( int iRow,
			    int extraNeeded )
{
  int number = numberInRow_[iRow];
  OsiBigIndex space = lengthAreaU_ - startRowU_[maximumRowsExtra_];

  if ( space < extraNeeded + number + 2 ) {
    //compression
    int iRow = nextRow_[maximumRowsExtra_];
    OsiBigIndex put = 0;

    while ( iRow != maximumRowsExtra_ ) {
      //move
      OsiBigIndex get = startRowU_[iRow];
      OsiBigIndex getEnd = startRowU_[iRow] + numberInRow_[iRow];

      startRowU_[iRow] = put;
      OsiBigIndex i;
      for ( i = get; i < getEnd; i++ ) {
	indexColumnU_[put] = indexColumnU_[i];
	put++;
      }				
      iRow = nextRow_[iRow];
      numberCompressions_++;
    }				/* endwhile */
    startRowU_[maximumRowsExtra_] = put;
    space = lengthAreaU_ - put;
    if ( space < extraNeeded + number + 2 ) {
      //need more space
      //if we can allocate bigger then do so and copy
      //if not then return so code can start again
      status_ = -99;
      return false;;
    }				
  }				
  OsiBigIndex put = startRowU_[maximumRowsExtra_];
  int next = nextRow_[iRow];
  int last = lastRow_[iRow];

  //out
  nextRow_[last] = next;
  lastRow_[next] = last;
  //in at end
  last = lastRow_[maximumRowsExtra_];
  nextRow_[last] = iRow;
  lastRow_[maximumRowsExtra_] = iRow;
  lastRow_[iRow] = last;
  nextRow_[iRow] = maximumRowsExtra_;
  //move
  OsiBigIndex get = startRowU_[iRow];

  startRowU_[iRow] = put;
  while ( number ) {
    number--;
    indexColumnU_[put] = indexColumnU_[get];
    put++;
    get++;
  }				/* endwhile */
  //add 4 for luck
  startRowU_[maximumRowsExtra_] = put + extraNeeded + 4;
  return true;
}

//  cleanup.  End of factorization
void
OsiFactorization::cleanup (  )
{
  getColumnSpace ( 0, INT_MAX >> 1 );	//compress
  OsiBigIndex lastU = startColumnU_[maximumColumnsExtra_];

  //free some memory here
  delete []  saveColumn_ ;
  delete []  markRow_ ;
  delete []  firstCount_ ;
  delete []  nextCount_ ;
  delete []  lastCount_ ;
  saveColumn_ = 0;
  markRow_ = 0;
  firstCount_ = 0;
  nextCount_ = 0;
  lastCount_ = 0;
  //make column starts OK
  //for best cache behavior get in order (last pivot at bottom of space)
  //that will need thinking about
  //use nextRow for permutation  (as that is what it is)
  int i;

  int * swap = permute_;
  permute_ = nextRow_;
  nextRow_ = swap;
  //safety feature
  permute_[numberRows_] = 0;
  permuteBack_ = new int [ maximumRowsExtra_ + 1 ];
#ifdef ZEROFAULT
  memset(permuteBack_,'w',(maximumRowsExtra_+1)*sizeof(int));
#endif
  for ( i = 0; i < numberRows_; i++ ) {
    int iRow = permute_[i];

    permuteBack_[iRow] = i;
  }				
  //redo nextRow_
  int extraSpace = maximumPivots_;

  for ( i = 0; i < numberColumns_; i++ ) {
    int number = numberInColumn_[i];	//always 0?
    int processed = numberInColumnPlus_[i];
    OsiBigIndex start = startColumnU_[i] - processed;

    number += processed;
    numberInColumn_[i] = number;
    startColumnU_[i] = start;
    //full list
    numberInColumnPlus_[i] = 0;
  }				
  int numberU = 0;

  pivotColumnBack_ = new int [ maximumRowsExtra_ + 1 ];
#ifdef ZEROFAULT
  memset(pivotColumnBack_,'q',(maximumRowsExtra_+1)*sizeof(int));
#endif
  for ( i = 0; i < numberColumns_; i++ ) {
    int iColumn = pivotColumn_[i];

    pivotColumnBack_[iColumn] = i;
    if ( iColumn >= 0 ) {
      if ( !numberInColumnPlus_[iColumn] ) {
	//wanted
	if ( numberU != iColumn ) {
	  numberInColumnPlus_[iColumn] = numberU;
	} else {
	  numberInColumnPlus_[iColumn] = -1;	//already in correct place
	}			
	numberU++;
      }				
    }				
  }				
  for ( i = 0; i < numberColumns_; i++ ) {
    int number = numberInColumn_[i];	//always 0?
    int where = numberInColumnPlus_[i];

    numberInColumnPlus_[i] = -1;
    OsiBigIndex start = startColumnU_[i];

    while ( where >= 0 ) {
      //put where it should be
      int numberNext = numberInColumn_[where];	//always 0?
      int whereNext = numberInColumnPlus_[where];
      OsiBigIndex startNext = startColumnU_[where];

      numberInColumn_[where] = number;
      numberInColumnPlus_[where] = -1;
      startColumnU_[where] = start;
      number = numberNext;
      where = whereNext;
      start = startNext;
    }				/* endwhile */
  }				
  //sort - using indexColumn
  CoinFillN ( indexColumnU_, lastU, -1 );
  OsiBigIndex k = 0;
  int *numberInColumn = numberInColumn_;
  int *indexColumnU = indexColumnU_;
  OsiBigIndex *startColumn = startColumnU_;
  int *indexRowU = indexRowU_;
  double *elementU = elementU_;

  for ( i = 0; i < numberRows_; i++ ) {
    OsiBigIndex start = startColumn[i];
    OsiBigIndex end = start + numberInColumn[i];

    OsiBigIndex j;
    for ( j = start; j < end; j++ ) {
      indexColumnU[j] = k++;
    }				
  }				
  for ( i = 0; i < numberRows_; i++ ) {
    OsiBigIndex start = startColumn[i];
    OsiBigIndex end = start + numberInColumn[i];

    OsiBigIndex j;
    for ( j = start; j < end; j++ ) {
      OsiBigIndex k = indexColumnU[j];
      int iRow = indexRowU[j];
      double element = elementU[j];

      while ( k != -1 ) {
	OsiBigIndex kNext = indexColumnU[k];
	int iRowNext = indexRowU[k];
	double elementNext = elementU[k];

	indexColumnU_[k] = -1;
	indexRowU[k] = iRow;
	elementU[k] = element;
	k = kNext;
	iRow = iRowNext;
	element = elementNext;
      }				/* endwhile */
    }				
  }				
  k = 0;
  for ( i = 0; i < numberRows_; i++ ) {
    startColumnU_[i] = k;
    k += numberInColumn_[i];
  }				
  delete []  numberInColumnPlus_ ;
  numberInColumnPlus_ = 0;
  numberU_ = numberU;
  numberGoodU_ = numberU;
  numberL_ = numberGoodL_;
#if OSI_DEBUG
  for ( i = 0; i < numberRows_; i++ ) {
    if ( permute_[i] < 0 ) {
      std::cout << i << std::endl;
      abort (  );
    }				
  }				
#endif
  for ( i = numberSlacks_; i < numberU; i++ ) {
    OsiBigIndex start = startColumnU_[i];
    OsiBigIndex end = start + numberInColumn_[i];

    totalElements_ += numberInColumn_[i];
    if ( end > start || pivotRegion_[i] != 1.0 ) {
      OsiBigIndex j;
      for ( j = start; j < end; j++ ) {
	int iRow = indexRowU_[j];

	iRow = permute_[iRow];
	indexRowU_[j] = iRow;
	numberInRow_[iRow]++;
      }		
    }				
  }				
  //space for cross reference
  convertRowToColumnU_ = new OsiBigIndex [ lengthAreaU_ ];
  OsiBigIndex *convertRowToColumn = convertRowToColumnU_;
  OsiBigIndex j = 0;
  OsiBigIndex *startRow = startRowU_;

  int iRow;
  for ( iRow = 0; iRow < numberRows_; iRow++ ) {
    startRow[iRow] = j;
    j += numberInRow_[iRow];
  }				
  OsiBigIndex numberInU = j;

  CoinFillN ( numberInRow_, numberRows_ , 0);
  OsiBigIndex lowCount = 0;
  OsiBigIndex highCount = numberInU;
  int lowC = 0;
  int highC = numberRows_ - numberSlacks_;

  for ( i = numberSlacks_; i < numberRows_; i++ ) {
    OsiBigIndex start = startColumnU_[i];
    OsiBigIndex end = start + numberInColumn_[i];

    lowCount += numberInColumn_[i];
    highCount -= numberInColumn_[i];
    lowC++;
    highC--;
    double pivotValue = pivotRegion_[i];

    OsiBigIndex j;
    for ( j = start; j < end; j++ ) {
      int iRow = indexRowU_[j];
      int iLook = numberInRow_[iRow];

      numberInRow_[iRow] = iLook + 1;
      OsiBigIndex k = startRow[iRow] + iLook;

      indexColumnU_[k] = i;
      convertRowToColumn[k] = j;
      //multiply by pivot
      elementU_[j] *= pivotValue;
    }				
  }				
  for ( j = 0; j < numberRows_; j++ ) {
    lastRow_[j] = j - 1;
    nextRow_[j] = j + 1;
  }				
  nextRow_[numberRows_ - 1] = maximumRowsExtra_;
  lastRow_[maximumRowsExtra_] = numberRows_ - 1;
  nextRow_[maximumRowsExtra_] = 0;
  lastRow_[0] = maximumRowsExtra_;
  startRow[maximumRowsExtra_] = numberInU;
  OsiBigIndex *startColumnL = startColumnL_;

  int firstReal = numberRows_;

  for ( i = numberRows_ - 1; i >= 0; i-- ) {
    OsiBigIndex start = startColumnL[i];
    OsiBigIndex end = startColumnL[i + 1];

    totalElements_ += end - start;
    int pivotRow = pivotRowL_[i];

    pivotRow = permute_[pivotRow];
    pivotRowL_[i] = pivotRow;
    if ( end > start ) {
      firstReal = i;
      OsiBigIndex j;
      for ( j = start; j < end; j++ ) {
	int iRow = indexRowL_[j];

	iRow = permute_[iRow];
	indexRowL_[j] = iRow;
      }				
    }				
  }				
  baseL_ = firstReal;
  numberL_ = numberGoodL_ - firstReal;
  factorElements_ = totalElements_;
  pivotRowL_[numberGoodL_] = numberRows_;	//so loop will be clean
  //can deletepivotRowL_ as not used
  delete []  pivotRowL_ ;
  pivotRowL_ = 0;
  startColumnR_ = new OsiBigIndex [ extraSpace + 1 ];
#ifdef ZEROFAULT
  memset(startColumnR_,'z',(extraSpace + 1)*sizeof(OsiBigIndex));
#endif
  //use L for R if room
  OsiBigIndex space = lengthAreaL_ - lengthL_;
  OsiBigIndex spaceUsed = lengthL_ + lengthU_;

  extraSpace = maximumPivots_;
  int needed = ( spaceUsed + numberRows_ - 1 ) / numberRows_;

  needed = needed * 2 * maximumPivots_;
  if ( needed < 2 * numberRows_ ) {
    needed = 2 * numberRows_;
  }				
  if ( space >= needed ) {
    lengthR_ = 0;
    lengthAreaR_ = space;
    elementR_ = elementL_ + lengthL_;
    indexRowR_ = indexRowL_ + lengthL_;
  } else {
    lengthR_ = 0;
    lengthAreaR_ = space;
    elementR_ = elementL_ + lengthL_;
    indexRowR_ = indexRowL_ + lengthL_;
    if ((messageLevel_&1))
      std::cout<<"Factorization may need some increasing area space"
	       <<std::endl;
    if ( areaFactor_ ) {
      areaFactor_ *= 1.1;
    } else {
      areaFactor_ = 1.1;
    }				
  }				
  numberR_ = 0;
}

//  checkConsistency.  Checks that row and column copies look OK
void
OsiFactorization::checkConsistency (  )
{
  bool bad = false;

  int iRow;
  for ( iRow = 0; iRow < numberRows_; iRow++ ) {
    if ( numberInRow_[iRow] ) {
      OsiBigIndex startRow = startRowU_[iRow];
      OsiBigIndex endRow = startRow + numberInRow_[iRow];

      OsiBigIndex j;
      for ( j = startRow; j < endRow; j++ ) {
	int iColumn = indexColumnU_[j];
	OsiBigIndex startColumn = startColumnU_[iColumn];
	OsiBigIndex endColumn = startColumn + numberInColumn_[iColumn];
	bool found = false;

	OsiBigIndex k;
	for ( k = startColumn; k < endColumn; k++ ) {
	  if ( indexRowU_[k] == iRow ) {
	    found = true;
	    break;
	  }			
	}			
	if ( !found ) {
	  bad = true;
	  std::cout << "row " << iRow << " column " << iColumn << " Rows" << std::endl;
	}			
      }				
    }				
  }				
  int iColumn;
  for ( iColumn = 0; iColumn < numberColumns_; iColumn++ ) {
    if ( numberInColumn_[iColumn] ) {
      OsiBigIndex startColumn = startColumnU_[iColumn];
      OsiBigIndex endColumn = startColumn + numberInColumn_[iColumn];

      OsiBigIndex j;
      for ( j = startColumn; j < endColumn; j++ ) {
	int iRow = indexRowU_[j];
	OsiBigIndex startRow = startRowU_[iRow];
	OsiBigIndex endRow = startRow + numberInRow_[iRow];
	bool found = false;

	OsiBigIndex k;
	for (  k = startRow; k < endRow; k++ ) {
	  if ( indexColumnU_[k] == iColumn ) {
	    found = true;
	    break;
	  }			
	}			
	if ( !found ) {
	  bad = true;
	  std::cout << "row " << iRow << " column " << iColumn << " Columns" <<
	    std::endl;
	}			
      }				
    }				
  }				
  if ( bad ) {
    abort (  );
  }				
}
  //  pivotOneOtherRow.  When just one other row so faster
bool 
OsiFactorization::pivotOneOtherRow ( int pivotRow,
					   int pivotColumn )
{
  int numberInPivotRow = numberInRow_[pivotRow] - 1;
  OsiBigIndex startColumn = startColumnU_[pivotColumn];
  OsiBigIndex startRow = startRowU_[pivotRow];
  OsiBigIndex endRow = startRow + numberInPivotRow + 1;

  //take out this bit of indexColumnU
  int next = nextRow_[pivotRow];
  int last = lastRow_[pivotRow];

  nextRow_[last] = next;
  lastRow_[next] = last;
  nextRow_[pivotRow] = numberGoodU_;	//use for permute
  lastRow_[pivotRow] = -2;
  numberInRow_[pivotRow] = 0;
  //store column in L, compress in U and take column out
  OsiBigIndex l = lengthL_;

  if ( l + 1 > lengthAreaL_ ) {
    //need another area
    std::cout << "code pivot 1" << std::endl;
    //leave gap so starts will work
    abort (  );
  }				
  //l+=currentAreaL_->elementByColumn-elementL_;
  //OsiBigIndex lSave=l;
  pivotRowL_[numberGoodL_] = pivotRow;
  startColumnL_[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL_[numberGoodL_] = l + 1;
  lengthL_++;
  double pivotElement;
  double otherMultiplier;
  int otherRow;

  if ( indexRowU_[startColumn] == pivotRow ) {
    pivotElement = elementU_[startColumn];
    otherMultiplier = elementU_[startColumn + 1];
    otherRow = indexRowU_[startColumn + 1];
  } else {
    pivotElement = elementU_[startColumn + 1];
    otherMultiplier = elementU_[startColumn];
    otherRow = indexRowU_[startColumn];
  }				
  int numberSave = numberInRow_[otherRow];
  double pivotMultiplier = 1.0 / pivotElement;

  pivotRegion_[numberGoodU_] = pivotMultiplier;
  numberInColumn_[pivotColumn] = 0;
  otherMultiplier = otherMultiplier * pivotMultiplier;
  indexRowL_[l] = otherRow;
  elementL_[l] = otherMultiplier;
  //take out of row list
  OsiBigIndex start = startRowU_[otherRow];
  OsiBigIndex end = start + numberSave;
  OsiBigIndex where = start;

  while ( indexColumnU_[where] != pivotColumn ) {
    where++;
  }				/* endwhile */
#if OSI_DEBUG
  if ( where >= end ) {
    abort (  );
  }				
#endif
  end--;
  indexColumnU_[where] = indexColumnU_[end];
  int numberAdded = 0;
  int numberDeleted = 0;

  //pack down and move to work
  int j;

  for ( j = startRow; j < endRow; j++ ) {
    int iColumn = indexColumnU_[j];

    if ( iColumn != pivotColumn ) {
      OsiBigIndex startColumn = startColumnU_[iColumn];
      OsiBigIndex endColumn = startColumn + numberInColumn_[iColumn];
      int iRow = indexRowU_[startColumn];
      double value = elementU_[startColumn];
      double largest;
      bool foundOther = false;

      //leave room for pivot
      OsiBigIndex put = startColumn + 1;
      OsiBigIndex positionLargest = -1;
      double thisPivotValue = 0.0;
      double otherElement = 0.0;
      double nextValue = elementU_[put];;
      int nextIRow = indexRowU_[put];

      //compress column and find largest not updated
      if ( iRow != pivotRow ) {
	if ( iRow != otherRow ) {
	  largest = fabs ( value );
	  elementU_[put] = value;
	  indexRowU_[put] = iRow;
	  positionLargest = put;
	  put++;
	  OsiBigIndex i;
	  for ( i = startColumn + 1; i < endColumn; i++ ) {
	    iRow = nextIRow;
	    value = nextValue;
#ifdef ZEROFAULT
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU_[i + 1];
	      nextValue = elementU_[i + 1];
#ifdef ZEROFAULT
	    }
#endif
	    if ( iRow != pivotRow ) {
	      if ( iRow != otherRow ) {
		//keep
		indexRowU_[put] = iRow;
		elementU_[put] = value;;
		put++;
	      } else {
		otherElement = value;
		foundOther = true;
	      }			
	    } else {
	      thisPivotValue = value;
	    }			
	  }			
	} else {
	  otherElement = value;
	  foundOther = true;
	  //need to find largest
	  largest = 0.0;
	  OsiBigIndex i;
	  for ( i = startColumn + 1; i < endColumn; i++ ) {
	    iRow = nextIRow;
	    value = nextValue;
#ifdef ZEROFAULT
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU_[i + 1];
	      nextValue = elementU_[i + 1];
#ifdef ZEROFAULT
	    }
#endif
	    if ( iRow != pivotRow ) {
	      //keep
	      indexRowU_[put] = iRow;
	      elementU_[put] = value;;
	      double absValue = fabs ( value );

	      if ( absValue > largest ) {
		largest = absValue;
		positionLargest = put;
	      }			
	      put++;
	    } else {
	      thisPivotValue = value;
	    }			
	  }			
	}			
      } else {
	//need to find largest
	largest = 0.0;
	thisPivotValue = value;
	OsiBigIndex i;
	for ( i = startColumn + 1; i < endColumn; i++ ) {
	  iRow = nextIRow;
	  value = nextValue;
#ifdef ZEROFAULT
	  // doesn't matter reading uninitialized but annoys checking
	  if ( i + 1 < endColumn ) {
#endif
	    nextIRow = indexRowU_[i + 1];
	    nextValue = elementU_[i + 1];
#ifdef ZEROFAULT
	  }
#endif
	  if ( iRow != otherRow ) {
	    //keep
	    indexRowU_[put] = iRow;
	    elementU_[put] = value;;
	    double absValue = fabs ( value );

	    if ( absValue > largest ) {
	      largest = absValue;
	      positionLargest = put;
	    }			
	    put++;
	  } else {
	    otherElement = value;
	    foundOther = true;
	  }			
	}			
      }				
      //slot in pivot
      elementU_[startColumn] = thisPivotValue;
      indexRowU_[startColumn] = pivotRow;
      //clean up counts
      startColumn++;
      numberInColumn_[iColumn] = put - startColumn;
      numberInColumnPlus_[iColumn]++;
      startColumnU_[iColumn]++;
      otherElement = otherElement - thisPivotValue * otherMultiplier;
      double absValue = fabs ( otherElement );

      if ( absValue > zeroTolerance_ ) {
	if ( !foundOther ) {
	  //have we space
	  saveColumn_[numberAdded++] = iColumn;
	  int next = nextColumn_[iColumn];
	  OsiBigIndex space;

	  space = startColumnU_[next] - put - numberInColumnPlus_[next];
	  if ( space <= 0 ) {
	    //getColumnSpace also moves fixed part
	    int number = numberInColumn_[iColumn];

	    if ( !getColumnSpace ( iColumn, number + 1 ) ) {
	      return false;
	    }
	    //redo starts
	    positionLargest =
	      positionLargest + startColumnU_[iColumn] - startColumn;
	    startColumn = startColumnU_[iColumn];
	    put = startColumn + number;
	  }			
	}			
	elementU_[put] = otherElement;
	indexRowU_[put] = otherRow;
	if ( absValue > largest ) {
	  largest = absValue;
	  positionLargest = put;
	}			
	put++;
      } else {
	if ( foundOther ) {
	  numberDeleted++;
	  //take out of row list
	  OsiBigIndex where = start;

	  while ( indexColumnU_[where] != iColumn ) {
	    where++;
	  }			/* endwhile */
#if OSI_DEBUG
	  if ( where >= end ) {
	    abort (  );
	  }			
#endif
	  end--;
	  indexColumnU_[where] = indexColumnU_[end];
	}			
      }				
      numberInColumn_[iColumn] = put - startColumn;
      //move largest
      if ( positionLargest >= 0 ) {
	value = elementU_[positionLargest];
	iRow = indexRowU_[positionLargest];
	elementU_[positionLargest] = elementU_[startColumn];
	indexRowU_[positionLargest] = indexRowU_[startColumn];
	elementU_[startColumn] = value;
	indexRowU_[startColumn] = iRow;
      }				
      //linked list for column
      if ( nextCount_[iColumn + numberRows_] != -2 ) {
	//modify linked list
	deleteLink ( iColumn + numberRows_ );
	addLink ( iColumn + numberRows_, numberInColumn_[iColumn] );
      }				
    }				
  }				
  //get space for row list
  next = nextRow_[otherRow];
  OsiBigIndex space;

  space = startRowU_[next] - end;
  totalElements_ += numberAdded - numberDeleted;
  int number = numberAdded + ( end - start );

  if ( space < numberAdded ) {
    numberInRow_[otherRow] = end - start;
    if ( !getRowSpace ( otherRow, number ) ) {
      return false;
    }
    end = startRowU_[otherRow] + end - start;
  }				
  // do linked lists and update counts
  numberInRow_[otherRow] = number;
  if ( number != numberSave ) {
    deleteLink ( otherRow );
    addLink ( otherRow, number );
  }				
  for ( j = 0; j < numberAdded; j++ ) {
    indexColumnU_[end++] = saveColumn_[j];
  }				
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  return true;
}
