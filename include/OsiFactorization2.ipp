{
  int *indexColumnU = indexColumnU_;
  OsiBigIndex *startColumnU = startColumnU_;
  int *numberInColumn = numberInColumn_;
  double *elementU = elementU_;
  int *indexRowU = indexRowU_;
  OsiBigIndex *startRowU = startRowU_;
  int *numberInRow = numberInRow_;
  double *elementL = elementL_;
  int *indexRowL = indexRowL_;
  int *saveColumn = saveColumn_;
  int *nextRow = nextRow_;
  int *lastRow = lastRow_;

  //store pivot columns (so can easily compress)
  int numberInPivotRow = numberInRow[pivotRow] - 1;
  OsiBigIndex startColumn = startColumnU[pivotColumn];
  int numberInPivotColumn = numberInColumn[pivotColumn] - 1;
  OsiBigIndex endColumn = startColumn + numberInPivotColumn + 1;
  int put = 0;
  OsiBigIndex startRow = startRowU[pivotRow];
  OsiBigIndex endRow = startRow + numberInPivotRow + 1;

  if ( pivotColumnPosition < 0 ) {
    OsiBigIndex i;
    for ( i = startRow; i < endRow; i++ ) {
      int iColumn = indexColumnU[i];

      if ( iColumn != pivotColumn ) {
	saveColumn[put] = iColumn;
	put++;
      } else {
	pivotColumnPosition = i;
      }				
    }				
  } else {
    OsiBigIndex i;

    for ( i = startRow; i < pivotColumnPosition; i++ ) {
      int iColumn = indexColumnU[i];

      saveColumn[put] = iColumn;
      put++;
    }				
    for ( i = pivotColumnPosition + 1; i < endRow; i++ ) {
      int iColumn = indexColumnU[i];

      saveColumn[put] = iColumn;
      put++;
    }				
  }				
  //take out this bit of indexColumnU
  int next = nextRow[pivotRow];
  int last = lastRow[pivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  nextRow[pivotRow] = numberGoodU_;	//use for permute
  lastRow[pivotRow] = -2;
  numberInRow[pivotRow] = 0;
  //store column in L, compress in U and take column out
  OsiBigIndex l = lengthL_;

  if ( l + numberInPivotColumn > lengthAreaL_ ) {
    //need another area
    cout << "code pivot 1" << endl;
    //leave gap so starts will work
    exit ( 99 );
  }				
  //l+=currentAreaL_->elementByColumn-elementL;
  OsiBigIndex lSave = l;

  pivotRowL_[numberGoodL_] = pivotRow;
  startColumnL_[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL_[numberGoodL_] = l + numberInPivotColumn;
  lengthL_ += numberInPivotColumn;
  if ( pivotRowPosition < 0 ) {
    OsiBigIndex i;
    for ( i = startColumn; i < endColumn; i++ ) {
      int iRow = indexRowU[i];

      if ( iRow != pivotRow ) {
	indexRowL[l] = iRow;
	elementL[l] = elementU[i];
	markRow[iRow] = l - lSave;
	l++;
	//take out of row list
	OsiBigIndex start = startRowU[iRow];
	OsiBigIndex end = start + numberInRow[iRow];
	OsiBigIndex where = start;

	while ( indexColumnU[where] != pivotColumn ) {
	  where++;
	}			/* endwhile */
#if DEBUG_OSI
	if ( where >= end ) {
	  abort (  );
	}			
#endif
	indexColumnU[where] = indexColumnU[end - 1];
	numberInRow[iRow]--;
      } else {
	pivotRowPosition = i;
      }				
    }				
  } else {
    OsiBigIndex i;

    for ( i = startColumn; i < pivotRowPosition; i++ ) {
      int iRow = indexRowU[i];

      markRow[iRow] = l - lSave;
      indexRowL[l] = iRow;
      elementL[l] = elementU[i];
      l++;
      //take out of row list
      OsiBigIndex start = startRowU[iRow];
      OsiBigIndex end = start + numberInRow[iRow];
      OsiBigIndex where = start;

      while ( indexColumnU[where] != pivotColumn ) {
	where++;
      }				/* endwhile */
#if DEBUG_OSI
      if ( where >= end ) {
	abort (  );
      }				
#endif
      indexColumnU[where] = indexColumnU[end - 1];
      numberInRow[iRow]--;
    }				
    for ( i = pivotRowPosition + 1; i < endColumn; i++ ) {
      int iRow = indexRowU[i];

      markRow[iRow] = l - lSave;
      indexRowL[l] = iRow;
      elementL[l] = elementU[i];
      l++;
      //take out of row list
      OsiBigIndex start = startRowU[iRow];
      OsiBigIndex end = start + numberInRow[iRow];
      OsiBigIndex where = start;

      while ( indexColumnU[where] != pivotColumn ) {
	where++;
      }				/* endwhile */
#if DEBUG_OSI
      if ( where >= end ) {
	abort (  );
      }				
#endif
      indexColumnU[where] = indexColumnU[end - 1];
      numberInRow[iRow]--;
    }				
  }				
  markRow[pivotRow] = LARGEINTEGER;
  //compress pivot column (move pivot to front including saved)
  double pivotElement = elementU[pivotRowPosition];
  double pivotMultiplier = 1.0 / pivotElement;

  pivotRegion_[numberGoodU_] = pivotMultiplier;
  numberInColumn[pivotColumn] = 0;
  //use end of L for temporary space
  int *indexL = &indexRowL[lSave];
  double *multipliersL = &elementL[lSave];

  //adjust
  int j;

  for ( j = 0; j < numberInPivotColumn; j++ ) {
    multipliersL[j] = pivotMultiplier * multipliersL[j];
  }				
  //zero out fill
  OsiBigIndex iErase;
  for ( iErase = 0; iErase < increment2 * numberInPivotRow;
	iErase++ ) {
    workArea2[iErase] = 0;
  }				
  OsiBigIndex added = numberInPivotRow * numberInPivotColumn;
  unsigned int *temp2 = workArea2;

  //pack down and move to work
  int jColumn;
  for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
    int iColumn = saveColumn[jColumn];
    OsiBigIndex startColumn = startColumnU[iColumn];
    OsiBigIndex endColumn = startColumn + numberInColumn[iColumn];
    int iRow = indexRowU[startColumn];
    double value = elementU[startColumn];
    double largest;
    OsiBigIndex put = startColumn;
    OsiBigIndex positionLargest = -1;
    double thisPivotValue = 0.0;

    //compress column and find largest not updated
    bool checkLargest;
    int mark = markRow[iRow];

    if ( mark < 0 ) {
      largest = fabs ( value );
      positionLargest = put;
      put++;
      checkLargest = false;
    } else {
      //need to find largest
      largest = 0.0;
      checkLargest = true;
      if ( mark != LARGEINTEGER ) {
	//will be updated
	work[mark] = value;
	int word = mark >> SHIFT_PER_INT;
	int bit = mark & MASK_PER_INT;

	temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	added--;
      } else {
	thisPivotValue = value;
      }				
    }				
    OsiBigIndex i;
    for ( i = startColumn + 1; i < endColumn; i++ ) {
      iRow = indexRowU[i];
      value = elementU[i];
      int mark = markRow[iRow];

      if ( mark < 0 ) {
	//keep
	indexRowU[put] = iRow;
	elementU[put] = value;;
	if ( checkLargest ) {
	  double absValue = fabs ( value );

	  if ( absValue > largest ) {
	    largest = absValue;
	    positionLargest = put;
	  }			
	}			
	put++;
      } else if ( mark != LARGEINTEGER ) {
	//will be updated
	work[mark] = value;;
	int word = mark >> SHIFT_PER_INT;
	int bit = mark & MASK_PER_INT;

	temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	added--;
      } else {
	thisPivotValue = value;
      }				
    }				
    //slot in pivot
    elementU[put] = elementU[startColumn];
    indexRowU[put] = indexRowU[startColumn];
    if ( positionLargest == startColumn ) {
      positionLargest = put;	//follow if was largest
    }				
    put++;
    elementU[startColumn] = thisPivotValue;
    indexRowU[startColumn] = pivotRow;
    //clean up counts
    startColumn++;
    numberInColumn[iColumn] = put - startColumn;
    numberInColumnPlus_[iColumn]++;
    startColumnU[iColumn]++;
    //how much space have we got
    int next = nextColumn_[iColumn];
    OsiBigIndex space;

    space = startColumnU[next] - put - numberInColumnPlus_[next];
    //assume no zero elements
    if ( numberInPivotColumn > space ) {
      //getColumnSpace also moves fixed part
      if ( !getColumnSpace ( iColumn, numberInPivotColumn ) ) {
	return false;
      }				
      //redo starts
      positionLargest = positionLargest + startColumnU[iColumn] - startColumn;
      startColumn = startColumnU[iColumn];
      put = startColumn + numberInColumn[iColumn];
    }				
    double tolerance = zeroTolerance_;

    for ( j = 0; j < numberInPivotColumn; j++ ) {
      value = work[j] - thisPivotValue * multipliersL[j];
      double absValue = fabs ( value );

      if ( absValue > tolerance ) {
	work[j] = 0.0;
	elementU[put] = value;
	indexRowU[put] = indexL[j];
	if ( absValue > largest ) {
	  largest = absValue;
	  positionLargest = put;
	}			
	put++;
      } else {
	work[j] = 0.0;
	added--;
	int word = j >> SHIFT_PER_INT;
	int bit = j & MASK_PER_INT;

	if ( temp2[word] & ( 1 << bit ) ) {
	  //take out of row list
	  iRow = indexL[j];
	  OsiBigIndex start = startRowU[iRow];
	  OsiBigIndex end = start + numberInRow[iRow];
	  OsiBigIndex where = start;

	  while ( indexColumnU[where] != iColumn ) {
	    where++;
	  }			/* endwhile */
#if DEBUG_OSI
	  if ( where >= end ) {
	    abort (  );
	  }			
#endif
	  indexColumnU[where] = indexColumnU[end - 1];
	  numberInRow[iRow]--;
	} else {
	  //make sure won't be added
	  int word = j >> SHIFT_PER_INT;
	  int bit = j & MASK_PER_INT;

	  temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	}			
      }				
    }				
    numberInColumn[iColumn] = put - startColumn;
    //move largest
    if ( positionLargest >= 0 ) {
      value = elementU[positionLargest];
      iRow = indexRowU[positionLargest];
      elementU[positionLargest] = elementU[startColumn];
      indexRowU[positionLargest] = indexRowU[startColumn];
      elementU[startColumn] = value;
      indexRowU[startColumn] = iRow;
    }				
    //linked list for column
    if ( nextCount_[iColumn + numberRows_] != -2 ) {
      //modify linked list
      deleteLink ( iColumn + numberRows_ );
      addLink ( iColumn + numberRows_, numberInColumn[iColumn] );
    }				
    temp2 += increment2;
  }				
  //get space for row list
  unsigned int *putBase = workArea2;
  int bigLoops = numberInPivotColumn >> SHIFT_PER_INT;
  int i = 0;

  // do linked lists and update counts
  while ( bigLoops ) {
    bigLoops--;
    int bit;
    for ( bit = 0; bit < BITS_PER_INT; i++, bit++ ) {
      unsigned int *putThis = putBase;
      int iRow = indexL[i];

      //get space
      int number = 0;
      int jColumn;

      for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
	unsigned int test = *putThis;

	putThis += increment2;
	test = 1 - ( ( test >> bit ) & 1 );
	number += test;
      }				
      int next = nextRow[iRow];
      OsiBigIndex space;

      space = startRowU[next] - startRowU[iRow];
      number += numberInRow[iRow];
      if ( space < number ) {
	if ( !getRowSpace ( iRow, number ) ) {
	  return false;
	}
      }				
      // now do
      putThis = putBase;
      next = nextRow[iRow];
      number = numberInRow[iRow];
      OsiBigIndex end = startRowU[iRow] + number;
      int saveIndex = indexColumnU[startRowU[next]];

      //add in
      for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
	unsigned int test = *putThis;

	putThis += increment2;
	test = 1 - ( ( test >> bit ) & 1 );
	indexColumnU[end] = saveColumn[jColumn];
	end += test;
      }				
      //put back next one in case zapped
      indexColumnU[startRowU[next]] = saveIndex;
      markRow[iRow] = -1;
      number = end - startRowU[iRow];
      numberInRow[iRow] = number;
      deleteLink ( iRow );
      addLink ( iRow, number );
    }				
    putBase++;
  }				/* endwhile */
  int bit;

  for ( bit = 0; i < numberInPivotColumn; i++, bit++ ) {
    unsigned int *putThis = putBase;
    int iRow = indexL[i];

    //get space
    int number = 0;
    int jColumn;

    for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
      unsigned int test = *putThis;

      putThis += increment2;
      test = 1 - ( ( test >> bit ) & 1 );
      number += test;
    }				
    int next = nextRow[iRow];
    OsiBigIndex space;

    space = startRowU[next] - startRowU[iRow];
    number += numberInRow[iRow];
    if ( space < number ) {
      if ( !getRowSpace ( iRow, number ) ) {
	return false;
      }
    }				
    // now do
    putThis = putBase;
    next = nextRow[iRow];
    number = numberInRow[iRow];
    OsiBigIndex end = startRowU[iRow] + number;
    int saveIndex;

#ifdef CHECKING
    if ( next != maximumRowsExtra_ )
#endif
      saveIndex = indexColumnU[startRowU[next]];

    //add in
    for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
      unsigned int test = *putThis;

      putThis += increment2;
      test = 1 - ( ( test >> bit ) & 1 );

      indexColumnU[end] = saveColumn[jColumn];
      end += test;
    }				
    //put back next one in case zapped
#ifdef CHECKING
    if ( next != maximumRowsExtra_ )
#endif
      indexColumnU[startRowU[next]] = saveIndex;
    markRow[iRow] = -1;
    number = end - startRowU[iRow];
    numberInRow[iRow] = number;
    deleteLink ( iRow );
    addLink ( iRow, number );
  }				
  markRow[pivotRow] = -1;
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  totalElements_ += added;
  return true;
}
