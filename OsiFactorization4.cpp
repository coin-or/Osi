#include "OsiFactorization.hpp"
#include "OsiIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"


//  updateColumnU.  Updates part of column (FTRANU)
void
OsiFactorization::updateColumnU ( OsiIndexedVector * regionSparse,
			      int * indices,
			      int numberIn ) const
{
  int numberNonZero = regionSparse->getNumElements (  );


  if (numberNonZero<sparseThreshold_) {
    updateColumnUSparse(regionSparse,indices,numberIn);
  } else {
    updateColumnUDensish ( regionSparse );
  }
}

//  unPermuteTranspose.  takes off permutation vector
//zeroes out region2 if desired
void
OsiFactorization::unPermuteTranspose ( OsiIndexedVector * regionSparse,
				   OsiIndexedVector * regionSparse2,
				   bool erase ) const
{
  //zero region
  regionSparse->clear (  );
  double *region = regionSparse->denseVector (  );
  double *region2 = regionSparse2->denseVector (  );
  int *regionIndex2 = regionSparse2->getIndices (  );
  int numberNonZero2 = regionSparse2->getNumElements (  );
  int *regionIndex = regionSparse->getIndices (  );

  if ( numberNonZero2 ) {
    int *permuteBack = permuteBack_;

    //permute indices and move
    if ( erase ) {
      int j;
      
      for ( j = 0; j <  numberNonZero2 ; j++ ) {
	int iRow = regionIndex2[j];
	int newRow = permuteBack[iRow];
	
	region[newRow] = region2[iRow];
	region2[iRow] = 0.0;
	regionIndex[j] = newRow;
      }			
      regionSparse2->setNumElements ( 0 );
#ifdef DEBUG
      regionSparse2->checkClean();
#endif
    } else {
      int j;
      
      for ( j = 0; j < numberNonZero2 ; j++ ) {
	int iRow = regionIndex2[j];
	int newRow = permuteBack[iRow];
	
	region[newRow] = region2[iRow];
	regionIndex[j] = newRow;
      }			
    }			
    regionSparse->setNumElements (  numberNonZero2 );
  }				
}

//  updateColumnUDensish.  Updates part of column (FTRANU)
void
OsiFactorization::updateColumnUDensish ( OsiIndexedVector * regionSparse) const
{

  double *region = regionSparse->denseVector (  );
  double tolerance = zeroTolerance_;
  OsiBigIndex *startColumn = startColumnU_;
  int *indexRow = indexRowU_;
  double *element = elementU_;
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero = 0;
  int *numberInColumn = numberInColumn_;
  int i;
  double *pivotRegion = pivotRegion_;

  for (i = numberU_-1 ; i >= numberSlacks_; i-- ) {
    double pivotValue = region[i];
    if ( fabs ( pivotValue ) > tolerance ) {
      OsiBigIndex start = startColumn[i];
      int end = start + numberInColumn[i];
      
      OsiBigIndex j;
      for ( j = start; j < end; j++ ) {
	double value = element[j];
	int iRow = indexRow[j];
	
	region[iRow] -= value * pivotValue;
      }
      region[i] = pivotValue*pivotRegion[i];
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }				
  }			

  // now do slacks
  double factor = slackValue_;
  if (factor==1.0) {
    for ( i = numberSlacks_-1; i>=0;i--) {				
      double absValue = fabs ( region[i] );
      if ( absValue ) {
	if ( absValue > tolerance ) {
	  regionIndex[numberNonZero++] = i;
	} else {
	  region[i] = 0.0;
	}				
      }				
    }
  } else {
    for ( i = numberSlacks_-1; i>=0;i--) {				
      double absValue = fabs ( region[i] );
      if ( absValue ) {
	if ( absValue > tolerance ) {
	  regionIndex[numberNonZero++] = i;
	  region [i] *= factor;
	} else {
	  region[i] = 0.0;
	}				
      }				
    }
  }
  regionSparse->setNumElements ( numberNonZero );
}
//  updateColumnU.  Updates part of column (FTRANU)
void
OsiFactorization::updateColumnUSparse ( OsiIndexedVector * regionSparse,
			      int * indices,
			      int numberIn ) const
{
  int numberNonZero = regionSparse->getNumElements (  );
  int *regionIndex = regionSparse->getIndices (  );
  double *region = regionSparse->denseVector (  );
  double tolerance = zeroTolerance_;
  OsiBigIndex *startColumn = startColumnU_;
  int *indexRow = indexRowU_;
  double *element = elementU_;
  double *pivotRegion = pivotRegion_;
  // use sparse_ as temporary area
  // mark known to be zero
  int * stack = sparse_;  /* pivot */
  int * list = stack + maximumRowsExtra_;  /* final list */
  int * next = list + maximumRowsExtra_;  /* jnext */
  char * mark = (char *) (next + maximumRowsExtra_);
  int nList, nStack;
  int i , iPivot;
#ifdef DEBUG
  for (i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif

  startColumn = startColumnU_;
  indexRow = indexRowU_;
  element = elementU_;
  // move slacks to end of stack list
  int * putLast = stack+maximumRowsExtra_;
  int * put = putLast;
  
  int *numberInColumn = numberInColumn_;
  nList = 0;
  if (!indices) 
    indices = regionIndex;
  for (i=0;i<numberIn;i++) {
    iPivot=indices[i];
    if (iPivot>=numberSlacks_) {
      nStack=1;
      stack[0]=iPivot;
      next[0]=startColumn[iPivot]+numberInColumn[iPivot]-1;
      while (nStack) {
	int kPivot,j;
	/* take off stack */
	kPivot=stack[--nStack];
	if (!mark[kPivot]) {
	  j=next[nStack];
	  if (j<startColumn[kPivot]) {
	    /* finished so mark */
	    list[nList++]=kPivot;
	    mark[kPivot]=1;
	  } else {
	    kPivot=indexRow[j];
	    /* put back on stack */
	    next[nStack++] --;
	    if (!mark[kPivot]) {
	      if (kPivot>=numberSlacks_) {
		/* and new one */
		stack[nStack]=kPivot;
		next[nStack++]=startColumn[kPivot]+numberInColumn[kPivot]-1;
	      } else {
		// slack
		mark[kPivot]=1;
		--put;
		*put=kPivot;
	      }
	    }
	  }
	}
      }
    } else if (!mark[iPivot]) {
      // slack
      --put;
      *put=iPivot;
      mark[iPivot]=1;
    }
  }
  numberNonZero=0;
  for (i=nList-1;i>=0;i--) {
    iPivot = list[i];
    mark[iPivot]=0;
    double pivotValue = region[iPivot];
    if ( fabs ( pivotValue ) > tolerance ) {
      region[iPivot] *= pivotRegion[iPivot];
      regionIndex[numberNonZero++]= iPivot;
      OsiBigIndex start = startColumn[iPivot];
      int number = numberInColumn[iPivot];
      
      OsiBigIndex j;
      for ( j = start; j < start+number; j++ ) {
	double value = element[j];
	int iRow = indexRow[j];
	
	region[iRow] -=  value * pivotValue;
      }		
    } else {
      region[iPivot]=0.0;
    }		
  }
  // slacks
  if (slackValue_==1.0) {
    for (;put<putLast;put++) {
      int iPivot = *put;
      mark[iPivot]=0;
      double pivotValue = region[iPivot];
      if ( fabs ( pivotValue ) > tolerance ) {
	regionIndex[numberNonZero++]= iPivot;
      } else {
	region[iPivot]=0.0;
      }
    }
  } else {
    for (;put<putLast;put++) {
      int iPivot = *put;
      mark[iPivot]=0;
      double pivotValue = region[iPivot];
      if ( fabs ( pivotValue ) > tolerance ) {
	region[iPivot] = pivotValue*slackValue_;
	regionIndex[numberNonZero++]= iPivot;
      } else {
	region[iPivot]=0.0;
      }
    }
  }
  regionSparse->setNumElements ( numberNonZero );
}

//  =
OsiFactorization & OsiFactorization::operator = ( const OsiFactorization & other ) {
  if (this != &other) {    
    gutsOfDestructor();
    gutsOfInitialize(3);
    gutsOfCopy(other);
  }
  return *this;
}
void OsiFactorization::gutsOfCopy(const OsiFactorization &other)
{
  elementU_ = new double [ other.lengthAreaU_ ];
  indexRowU_ = new int [ other.lengthAreaU_ ];
  indexColumnU_ = new int [ other.lengthAreaU_ ];
  convertRowToColumnU_ = new OsiBigIndex [ other.lengthAreaU_ ];
  if (other.sparseThreshold_) {
    elementByRowL_ = new double [ other.lengthAreaL_ ];
    indexColumnL_ = new int [ other.lengthAreaL_ ];
    startRowL_ = new OsiBigIndex [other.numberRows_+1];
  }
  elementL_ = new double [ other.lengthAreaL_ ];
  indexRowL_ = new int [ other.lengthAreaL_ ];
  startColumnL_ = new OsiBigIndex [ other.numberRows_ + 1 ];
  int extraSpace = other.maximumColumnsExtra_ - other.numberColumns_;

  startColumnR_ = new OsiBigIndex [ extraSpace + 1 ];
  startRowU_ = new OsiBigIndex [ other.maximumRowsExtra_ + 1 ];
  numberInRow_ = new int [ other.maximumRowsExtra_ + 1 ];
  nextRow_ = new int [ other.maximumRowsExtra_ + 1 ];
  lastRow_ = new int [ other.maximumRowsExtra_ + 1 ];
  pivotRegion_ = new double [ other.maximumRowsExtra_ + 1 ];
  permuteBack_ = new int [ other.maximumRowsExtra_ + 1 ];
  permute_ = new int [ other.maximumRowsExtra_ + 1 ];
  pivotColumnBack_ = new int [ other.maximumRowsExtra_ + 1 ];
  startColumnU_ = new OsiBigIndex [ other.maximumColumnsExtra_ + 1 ];
  numberInColumn_ = new int [ other.maximumColumnsExtra_ + 1 ];
  pivotColumn_ = new int [ other.maximumColumnsExtra_ + 1 ];
  nextColumn_ = new int [ other.maximumColumnsExtra_ + 1 ];
  lastColumn_ = new int [ other.maximumColumnsExtra_ + 1 ];
  numberTrials_ = other.numberTrials_;
  biggerDimension_ = other.biggerDimension_;
  numberSlacks_ = other.numberSlacks_;
  numberU_ = other.numberU_;
  lengthU_ = other.lengthU_;
  lengthAreaU_ = other.lengthAreaU_;
  numberL_ = other.numberL_;
  baseL_ = other.baseL_;
  lengthL_ = other.lengthL_;
  lengthAreaL_ = other.lengthAreaL_;
  numberR_ = other.numberR_;
  lengthR_ = other.lengthR_;
  lengthAreaR_ = other.lengthAreaR_;
  pivotTolerance_ = other.pivotTolerance_;
  zeroTolerance_ = other.zeroTolerance_;
  slackValue_ = other.slackValue_;
  areaFactor_ = other.areaFactor_;
  numberRows_ = other.numberRows_;
  numberRowsExtra_ = other.numberRowsExtra_;
  maximumRowsExtra_ = other.maximumRowsExtra_;
  numberColumns_ = other.numberColumns_;
  numberColumnsExtra_ = other.numberColumnsExtra_;
  maximumColumnsExtra_ = other.maximumColumnsExtra_;
  maximumPivots_=other.maximumPivots_;
  numberGoodU_ = other.numberGoodU_;
  numberGoodL_ = other.numberGoodL_;
  numberPivots_ = other.numberPivots_;
  messageLevel_ = other.messageLevel_;
  totalElements_ = other.totalElements_;
  factorElements_ = other.factorElements_;
  status_ = other.status_;
  doForrestTomlin_ = other.doForrestTomlin_;
  increasingRows_ = other.increasingRows_;
  sparseThreshold_=other.sparseThreshold_;
  OsiBigIndex space = lengthAreaL_ - lengthL_;

  lengthAreaR_ = space;
  elementR_ = elementL_ + lengthL_;
  indexRowR_ = indexRowL_ + lengthL_;
  //now copy everything
  //assuming numberRowsExtra==numberColumnsExtra
  CoinDisjointCopyN ( other.startRowU_, numberRowsExtra_ + 1, startRowU_ );
  CoinDisjointCopyN ( other.numberInRow_, numberRowsExtra_ + 1, numberInRow_ );
  CoinDisjointCopyN ( other.nextRow_, numberRowsExtra_ + 1, nextRow_ );
  CoinDisjointCopyN ( other.lastRow_, numberRowsExtra_ + 1, lastRow_ );
  CoinDisjointCopyN ( other.pivotRegion_, numberRowsExtra_ + 1, pivotRegion_ );
  CoinDisjointCopyN ( other.permuteBack_, numberRowsExtra_ + 1, permuteBack_ );
  CoinDisjointCopyN ( other.permute_, numberRowsExtra_ + 1, permute_ );
  CoinDisjointCopyN ( other.pivotColumnBack_, numberRowsExtra_ + 1, pivotColumnBack_ );
  CoinDisjointCopyN ( other.startColumnU_, numberRowsExtra_ + 1, startColumnU_ );
  CoinDisjointCopyN ( other.numberInColumn_, numberRowsExtra_ + 1, numberInColumn_ );
  CoinDisjointCopyN ( other.pivotColumn_, numberRowsExtra_ + 1, pivotColumn_ );
  CoinDisjointCopyN ( other.nextColumn_, numberRowsExtra_ + 1, nextColumn_ );
  CoinDisjointCopyN ( other.lastColumn_, numberRowsExtra_ + 1, lastColumn_ );
  CoinDisjointCopyN ( other.startColumnR_ , numberRowsExtra_ - numberColumns_ + 1,
	    startColumnR_ );  
  CoinDisjointCopyN ( other.elementR_, lengthR_, elementR_ );
  CoinDisjointCopyN ( other.indexRowR_, lengthR_, indexRowR_ );
  //extra one at end
  startColumnU_[maximumColumnsExtra_] =
    other.startColumnU_[maximumColumnsExtra_];
  nextColumn_[maximumColumnsExtra_] = other.nextColumn_[maximumColumnsExtra_];
  lastColumn_[maximumColumnsExtra_] = other.lastColumn_[maximumColumnsExtra_];
  startRowU_[maximumRowsExtra_] = other.startRowU_[maximumRowsExtra_];
  nextRow_[maximumRowsExtra_] = other.nextRow_[maximumRowsExtra_];
  lastRow_[maximumRowsExtra_] = other.lastRow_[maximumRowsExtra_];
  //row and column copies of U
  int iRow;
  for ( iRow = 0; iRow < numberRowsExtra_; iRow++ ) {
    //row
    OsiBigIndex start = startRowU_[iRow];
    int numberIn = numberInRow_[iRow];

    CoinDisjointCopyN ( other.indexColumnU_ + start, numberIn, indexColumnU_ + start );
    CoinDisjointCopyN (other.convertRowToColumnU_ + start , numberIn,
	     convertRowToColumnU_ + start );
    //column
    start = startColumnU_[iRow];
    numberIn = numberInColumn_[iRow];
    CoinDisjointCopyN ( other.indexRowU_ + start, numberIn, indexRowU_ + start );
    CoinDisjointCopyN ( other.elementU_ + start, numberIn, elementU_ + start );
  }				
  // L is contiguous
  CoinDisjointCopyN ( other.elementL_, lengthL_, elementL_ );
  CoinDisjointCopyN ( other.indexRowL_, lengthL_, indexRowL_ );
  CoinDisjointCopyN ( other.startColumnL_, numberRows_ + 1, startColumnL_ );
  if (other.sparseThreshold_) {
    makeRowCopyL();
  }
}
//  updateColumnR.  Updates part of column (FTRANR)
void
OsiFactorization::updateColumnR ( OsiIndexedVector * regionSparse ) const
{
  double *region = regionSparse->denseVector (  );
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero = regionSparse->getNumElements (  );

  if ( !numberR_ )
    return;	//return if nothing to do
  double tolerance = zeroTolerance_;

  int * indexRow = indexRowR_;
  double * element = elementR_;
  OsiBigIndex * startColumn = startColumnR_-numberRows_;

  int iRow;
  double pivotValue;

  int i;
  for ( i = numberRows_; i < numberRowsExtra_; i++ ) {
    //move using permute_ (stored in inverse fashion)
    iRow = permute_[i];
    pivotValue = region[iRow];
    //zero out pre-permuted
    region[iRow] = 0.0;

    OsiBigIndex j;
    for ( j = startColumn[i]; j < startColumn[i+1]; j ++ ) {
      double value = element[j];
      int jRow = indexRow[j];
      pivotValue = pivotValue - value * region[jRow];
    }				
    if ( fabs ( pivotValue ) > tolerance ) {
      region[i] = pivotValue;
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }				
  }				
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}

//  updateColumnTransposeR.  Updates part of column (FTRANR)
void
OsiFactorization::updateColumnTransposeR ( OsiIndexedVector * regionSparse ) const
{
  double *region = regionSparse->denseVector (  );
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero = regionSparse->getNumElements (  );
  double tolerance = zeroTolerance_;
  int i;

#ifdef DEBUG
  regionSparse->checkClean();
#endif
  if (numberNonZero) {
    int last = numberRowsExtra_-1;

    
    int *indexRow = indexRowR_;
    double *element = elementR_;
    OsiBigIndex * startColumn = startColumnR_-numberRows_;
    //move using permute_ (stored in inverse fashion)
    int putRow;
    double pivotValue;
    
    if (numberNonZero < sparseThreshold_) {
      // we can use sparse_ as temporary array
      int * spare = sparse_;
      for (i=0;i<numberNonZero;i++) {
	spare[regionIndex[i]]=i;
      }
      // still need to do in correct order (for now)
      for ( i = last ; i >= numberRows_; i-- ) {
	putRow = permute_[i];
	pivotValue = region[i];
	//zero out  old permuted
	region[i] = 0.0;
	if ( pivotValue ) {
	  OsiBigIndex j;
	  for ( j = startColumn[i]; j < startColumn[i+1]; j++ ) {
	    double value = element[j];
	    int iRow = indexRow[j];
	    double oldValue = region[iRow];
	    double newValue = oldValue - value * pivotValue;
	    if (oldValue) {
	      if (!newValue)
		newValue=1.0e-100;
	      region[iRow]=newValue;
	    } else if (fabs(newValue)>tolerance) {
	      region[iRow] = newValue;
	      spare[iRow]=numberNonZero;
	      regionIndex[numberNonZero++]=iRow;
	    }
	  }				
	  region[putRow] = pivotValue;
	  // modify list
	  int position=spare[i];
	  regionIndex[position]=putRow;
	  spare[putRow]=position;
	}				
      }				
      regionSparse->setNumElements(numberNonZero);
    } else {
      for ( i = last ; i >= numberRows_; i-- ) {
	putRow = permute_[i];
	pivotValue = region[i];
	//zero out  old permuted
	region[i] = 0.0;
	if ( pivotValue ) {
	  OsiBigIndex j;
	  for ( j = startColumn[i]; j < startColumn[i+1]; j++ ) {
	    double value = element[j];
	    int iRow = indexRow[j];
	    region[iRow] -= value * pivotValue;
	  }				
	  region[putRow] = pivotValue;
	  //putRow must have been zero before so put on list ??
	  //but can't catch up so will have to do L from end
	  //unless we update lookBtran in replaceColumn - yes
	}				
      }
    }				
  }
}
/* Updates one column (FTRAN) from region2 and permutes.
   region1 starts as zero
   If increasingRows_>1
   - returns permuted result in region1 and region2 is zero.
   otherwise
   - returns un-permuted result in region2 and region1 is zero */
int OsiFactorization::updateColumn ( OsiIndexedVector * regionSparse,
					     OsiIndexedVector * regionSparse2,
					     bool FTUpdate)
{
  //permute and move indices into index array
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero = regionSparse2->getNumElements();
  int *permute = permute_;
  int * index = regionSparse2->getIndices();
  double * region = regionSparse->denseVector();
  double * array = regionSparse2->denseVector();

  int j;
  for ( j = 0; j < numberNonZero; j ++ ) {
    int iRow = index[j];
    double value = array[iRow];
    array[iRow]=0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }				
  regionSparse->setNumElements ( numberNonZero );
  // will be negative if no room
  int number=updateColumn ( regionSparse, FTUpdate );
  if (increasingRows_>1) {
    // say region2 empty
    regionSparse2->setNumElements(0);
  } else {
    // permute back
    numberNonZero = regionSparse->getNumElements();
    int * permuteBack = pivotColumnBack_;
    for ( j = 0; j < numberNonZero; j ++ ) {
      int iRow = regionIndex[j];
      double value = region[iRow];
      region[iRow]=0.0;
      iRow = permuteBack[iRow];
      array[iRow] = value;
      index[j] = iRow;
    }			
    regionSparse2->setNumElements(numberNonZero);
  }
  return number;
}
/* Updates one column (FTRAN) to/from array
   This assumes user is thinking non-permuted
   - returns un-permuted result in array.
   region starts as zero and is zero at end */
int OsiFactorization::updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
			int number,
			bool FTUpdate ) 
{
  //permute and move indices into index array
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero;
  int *permute = permute_;
  double * region = regionSparse->denseVector();

  int j;
  for ( j = 0; j < number; j ++ ) {
    int iRow = index[j];
    double value = array[iRow];
    array[iRow]=0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }				
  regionSparse->setNumElements ( number );
  // if no room will return negative
  numberNonZero = updateColumn ( regionSparse, FTUpdate );
  // permute back
  number = regionSparse->getNumElements();
  int * permuteBack = pivotColumnBack_;
  for ( j = 0; j < number; j ++ ) {
    int iRow = regionIndex[j];
    double value = region[iRow];
    region[iRow]=0.0;
    iRow = permuteBack[iRow];
    array[iRow] = value;
    index[j] = iRow;
  }			
  regionSparse->setNumElements(0);
  return numberNonZero;
}

// const version
int OsiFactorization::updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
				    int number) const
{
  //permute and move indices into index array
  int *regionIndex = regionSparse->getIndices (  );
  int numberNonZero;
  int *permute = permute_;
  double * region = regionSparse->denseVector();

  int j;
  for ( j = 0; j < number; j ++ ) {
    int iRow = index[j];
    double value = array[iRow];
    array[iRow]=0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }				
  regionSparse->setNumElements ( number );
  // if no room will return negative
  numberNonZero = updateColumn ( regionSparse );
  // permute back
  number = regionSparse->getNumElements();
  int * permuteBack = pivotColumnBack_;
  for ( j = 0; j < number; j ++ ) {
    int iRow = regionIndex[j];
    double value = region[iRow];
    region[iRow]=0.0;
    iRow = permuteBack[iRow];
    array[iRow] = value;
    index[j] = iRow;
  }			
  regionSparse->setNumElements(0);
  return numberNonZero;
}
//  makes a row copy of L
void
OsiFactorization::makeRowCopyL ( )
{
  if (!sparseThreshold_)
    sparseThreshold_=(numberRows_+9)/10;
  //sparseThreshold_=99999;
  // allow for stack, list, next and char map of mark
  int nRowIndex = (maximumRowsExtra_+sizeof(int)-1)/
    sizeof(char);
  delete []sparse_;
  sparse_ = new int [ 3*maximumRowsExtra_ + nRowIndex ];
  // zero out mark
  memset(sparse_+3*maximumRowsExtra_,0,maximumRowsExtra_*sizeof(char));
  delete []elementByRowL_;
  delete []startRowL_;
  delete []indexColumnL_;
  elementByRowL_=new double[lengthAreaL_];
  startRowL_=new OsiBigIndex[numberRows_+1];
  indexColumnL_=new int[lengthAreaL_];
  // counts
  CoinFillN(startRowL_,numberRows_,0);
  int i;
  for (i=baseL_;i<baseL_+numberL_;i++) {
    OsiBigIndex j;
    for (j=startColumnL_[i];j<startColumnL_[i+1];j++) {
      int iRow = indexRowL_[j];
      startRowL_[iRow]++;
    }
  }
  // convert count to lasts
  OsiBigIndex count=0;
  for (i=0;i<numberRows_;i++) {
    int numberInRow=startRowL_[i];
    count += numberInRow;
    startRowL_[i]=count;
  }
  startRowL_[numberRows_]=count;
  // now insert
  for (i=baseL_+numberL_-1;i>=baseL_;i--) {
    OsiBigIndex j;
    for (j=startColumnL_[i];j<startColumnL_[i+1];j++) {
      int iRow = indexRowL_[j];
      OsiBigIndex start = startRowL_[iRow]-1;
      startRowL_[iRow]=start;
      elementByRowL_[start]=elementL_[j];
      indexColumnL_[start]=i;
    }
  }
}
//  get sparse threshold
int
OsiFactorization::sparseThreshold ( ) const
{
  return sparseThreshold_;
}

//  set sparse threshold
void
OsiFactorization::sparseThreshold ( int value ) 
{
  if (value>0&&sparseThreshold_) {
    sparseThreshold_=value;
  } else if (!value&&sparseThreshold_) {
    // delete copy
    sparseThreshold_=0;
    delete []elementByRowL_;
    delete []startRowL_;
    delete []indexColumnL_;
    elementByRowL_=NULL;
    startRowL_=NULL;
    indexColumnL_=NULL;
    delete []sparse_;
    sparse_=NULL;
  } else if (value>0&&!sparseThreshold_) {
    sparseThreshold_=value;
    makeRowCopyL();
  }
}
void OsiFactorization::maximumPivots (  int value )
{
  if (value>0) {
    maximumPivots_=value;
  }
}
void OsiFactorization::messageLevel (  int value )
{
  if (value>0&&value<16) {
    messageLevel_=value;
  }
}
void OsiFactorization::pivotTolerance (  double value )
{
  if (value>0.0&&value<=1.0) {
    pivotTolerance_=value;
  }
}
void OsiFactorization::zeroTolerance (  double value )
{
  if (value>0.0&&value<1.0) {
    zeroTolerance_=value;
  }
}
void OsiFactorization::slackValue (  double value )
{
  if (value>=0.0) {
    slackValue_=1.0;
  } else {
    slackValue_=-1.0;
  }
}
