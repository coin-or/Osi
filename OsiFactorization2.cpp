#include "OsiFactorization.hpp"
#include "OsiIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
//  factorSparse.  Does sparse phase of factorization
//return code is <0 error, 0= finished
int
OsiFactorization::factorSparse (  )
{
  int status = 0;
  int *indexRow = indexRowU_;
  int *indexColumn = indexColumnU_;
  double *element = elementU_;
  int count = 1;
  double *workArea = new double [ numberRows_ ];

  CoinFillN ( workArea, numberRows_ , 0.0);
  //get space for bit work area
  OsiBigIndex workSize = 1000;
  unsigned int *workArea2 = ( unsigned int * ) new int [ workSize ];
  int lastColumnInBlock;

  lastColumnInBlock = numberColumns_;
  int larger;

  if ( numberRows_ < numberColumns_ ) {
    larger = numberColumns_;
  } else {
    larger = numberRows_;
  }
  //do slacks first
  int pivotColumn;
  for ( pivotColumn = 0; pivotColumn < lastColumnInBlock;
	pivotColumn++ ) {
    if ( numberInColumn_[pivotColumn] == 1 ) {
      OsiBigIndex start = startColumnU_[pivotColumn];
      double value = element[start];

      if ( value == slackValue_ && numberInColumnPlus_[pivotColumn] == 0 ) {
	//treat as slack
	int iRow = indexRowU_[start];

	totalElements_--;
	if ( !pivotColumnSingleton ( iRow, pivotColumn ) ) {
	  return -99;
	}
	pivotColumn_[numberGoodU_] = pivotColumn;
	numberGoodU_++;
      }
    }
  }
  numberSlacks_ = numberGoodU_;
  while ( count <= biggerDimension_ ) {
    OsiBigIndex minimumCount = INT_MAX;
    OsiBigIndex minimumCost = INT_MAX;

    count = 1;
    bool stopping = false;
    int pivotRow = -1;
    int pivotColumn = -1;
    int pivotRowPosition = -1;
    int pivotColumnPosition = -1;
    int look = firstCount_[count];
    int trials = 0;

    while ( !stopping ) {
      if ( count == 1 && firstCount_[1] >= 0 ) {
	//do column singletons first to put more in U
	while ( look >= 0 ) {
	  if ( look >= numberRows_ ) {
	    int iColumn = look - numberRows_;

#if OSI_DEBUG
	    if ( numberInColumn_[iColumn] != count ) {
	      abort (  );
	    }
#endif
	    OsiBigIndex start = startColumnU_[iColumn];
	    int iRow = indexRow[start];

	    pivotRow = iRow;
	    pivotRowPosition = start;
	    pivotColumn = iColumn;
	    pivotColumnPosition = -1;
	    stopping = true;
	    look = -1;
	    break;
	  } else {
	    look = nextCount_[look];
	  }
	}			/* endwhile */
	if ( !stopping ) {
	  //back to singletons
	  look = firstCount_[1];
	}
      }
      int *nextCount = nextCount_;
      int *numberInRow = numberInRow_;
      int *numberInColumn = numberInColumn_;
      OsiBigIndex *startRow = startRowU_;
      OsiBigIndex *startColumn = startColumnU_;
      double pivotTolerance = pivotTolerance_;
      int numberTrials = numberTrials_;
      int numberRows = numberRows_;

      while ( look >= 0 ) {
	if ( look < numberRows_ ) {
	  int iRow = look;

#if OSI_DEBUG
	  if ( numberInRow[iRow] != count ) {
	    abort (  );
	  }
#endif
	  look = nextCount[look];
	  bool rejected = false;
	  OsiBigIndex start = startRow[iRow];
	  OsiBigIndex end = start + count;

	  OsiBigIndex i;
	  for ( i = start; i < end; i++ ) {
	    int iColumn = indexColumn[i];
	    OsiBigIndex cost = ( count - 1 ) * numberInColumn[iColumn];

	    if ( cost < minimumCost ) {
	      OsiBigIndex where = startColumn[iColumn];
	      double minimumValue = element[where];

	      minimumValue = fabs ( minimumValue ) * pivotTolerance;
	      while ( indexRow[where] != iRow ) {
		where++;
	      }			/* endwhile */
#if OSI_DEBUG
	      {
		OsiBigIndex end_debug = startColumn[iColumn] +
		  numberInColumn[iColumn];

		if ( where >= end_debug ) {
		  abort (  );
		}
	      }
#endif
	      double value = element[where];

	      value = fabs ( value );
	      if ( value >= minimumValue ) {
		minimumCost = cost;
		minimumCount = numberInColumn[iColumn];
		pivotRow = iRow;
		pivotRowPosition = -1;
		pivotColumn = iColumn;
		pivotColumnPosition = i;
		if ( minimumCount < count ) {
		  stopping = true;
		  look = -1;
		  break;
		}
	      } else if ( pivotRow == -1 ) {
		rejected = true;
	      }
	    }
	  }
	  trials++;
	  if ( trials >= numberTrials && pivotRow >= 0 ) {
	    stopping = true;
	    look = -1;
	    break;
	  }
	  if ( rejected ) {
	    //take out for moment
	    //eligible when row changes
	    deleteLink ( iRow );
	    addLink ( iRow, biggerDimension_ + 1 );
	  }
	} else {
	  int iColumn = look - numberRows;

#if OSI_DEBUG
	  if ( numberInColumn[iColumn] != count ) {
	    abort (  );
	  }
#endif
	  look = nextCount[look];
	  OsiBigIndex start = startColumn[iColumn];
	  OsiBigIndex end = start + numberInColumn[iColumn];
	  double minimumValue = element[start];

	  minimumValue = fabs ( minimumValue ) * pivotTolerance;
	  OsiBigIndex i;
	  for ( i = start; i < end; i++ ) {
	    double value = element[i];

	    value = fabs ( value );
	    if ( value >= minimumValue ) {
	      int iRow = indexRow[i];
	      OsiBigIndex cost = ( count - 1 ) * numberInRow[iRow];

	      if ( cost < minimumCost ) {
		minimumCost = cost;
		minimumCount = numberInRow[iRow];
		pivotRow = iRow;
		pivotRowPosition = i;
		pivotColumn = iColumn;
		pivotColumnPosition = -1;
		if ( minimumCount <= count + 1 ) {
		  stopping = true;
		  look = -1;
		  break;
		}
	      }
	    }
	  }
	  trials++;
	  if ( trials >= numberTrials && pivotRow >= 0 ) {
	    stopping = true;
	    look = -1;
	    break;
	  }
	}
      }				/* endwhile */
      //end of this - onto next
      if ( !stopping ) {
	count++;
	if ( count <= biggerDimension_ ) {
	  look = firstCount_[count];
	} else {
	  stopping = true;
	}
      } else {
	if ( pivotRow >= 0 ) {
	  int numberDoRow = numberInRow_[pivotRow] - 1;
	  int numberDoColumn = numberInColumn_[pivotColumn] - 1;

	  totalElements_ -= ( numberDoRow + numberDoColumn + 1 );
	  if ( numberDoColumn > 0 ) {
	    if ( numberDoRow > 0 ) {
	      if ( numberDoColumn > 1 ) {
		//  if (1) {
		//need to adjust more for cache and SMP
		//allow at least 4 extra
		int increment = numberDoColumn + 1 + 4;

		if ( increment & 15 ) {
		  increment = increment & ( ~15 );
		  increment += 16;
		}
		int increment2 =

		  ( increment + OSIFACTORIZATION_BITS_PER_INT - 1 ) >> OSIFACTORIZATION_SHIFT_PER_INT;
		OsiBigIndex size = increment2 * numberDoRow;

		if ( size > workSize ) {
		  workSize = size;
		  delete []  workArea2 ;
		  workArea2 = ( unsigned int * ) new int [ workSize ];
		}
		bool goodPivot;

		if ( larger < 32766 ) {
		  //branch out to best pivot routine 
		  goodPivot = pivot ( pivotRow, pivotColumn,
				      pivotRowPosition, pivotColumnPosition,
				      workArea, workArea2, increment,
				      increment2, ( short * ) markRow_ ,
				      32767);
		} else {
		  //might be able to do better by permuting
		  goodPivot = pivot ( pivotRow, pivotColumn,
				      pivotRowPosition, pivotColumnPosition,
				      workArea, workArea2, increment,
				      increment2, ( int * ) markRow_ ,
				      INT_MAX);
		}
		if ( !goodPivot ) {
		  return -99;
		}
	      } else {
		if ( !pivotOneOtherRow ( pivotRow, pivotColumn ) ) {
		  return -99;
		}
	      }
	    } else {
	      if ( !pivotRowSingleton ( pivotRow, pivotColumn ) ) {
		return -99;
	      }
	    }
	  } else {
	    if ( !pivotColumnSingleton ( pivotRow, pivotColumn ) ) {
	      return -99;
	    }
	  }
	  pivotColumn_[numberGoodU_] = pivotColumn;
	  numberGoodU_++;
	}
      }
    }				/* endwhile */
#if OSI_DEBUG==2
    checkConsistency (  );
#endif
  }				/* endwhile */
  delete []  workArea ;
  delete [] workArea2 ;
  return status;
}
