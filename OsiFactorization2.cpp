#include "OsiFactorization.hpp"
#include "OsiIndexedVector.hpp"

//  pivot.  Does one pivot in factorization
//will have several versions depending on markRow
bool
OsiFactorization::pivot ( int pivotRow,
			  int pivotColumn,
			  int pivotRowPosition,
			  int pivotColumnPosition,
			  double work[],
			  unsigned int workArea2[],
			  int increment,
			  int increment2,
			  int markRow[] )
#define LARGEINTEGER INT_MAX
#include "OsiFactorization2.ipp"
bool 
OsiFactorization::pivot ( int pivotRow,
			  int pivotColumn,
			  int pivotRowPosition,
			  int pivotColumnPosition,
			  double work[],
			  unsigned int workArea2[],
			  int increment,
			  int increment2,
			  short markRow[] )
#undef LARGEINTEGER
#define LARGEINTEGER 32767
#include "OsiFactorization2.ipp"
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
    cout << "code pivot 1" << endl;
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
#if DEBUG_OSI
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
#ifdef CHECKING
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU_[i + 1];
	      nextValue = elementU_[i + 1];
#ifdef CHECKING
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
#ifdef CHECKING
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU_[i + 1];
	      nextValue = elementU_[i + 1];
#ifdef CHECKING
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
#ifdef CHECKING
	  // doesn't matter reading uninitialized but annoys checking
	  if ( i + 1 < endColumn ) {
#endif
	    nextIRow = indexRowU_[i + 1];
	    nextValue = elementU_[i + 1];
#ifdef CHECKING
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
#if DEBUG_OSI
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
