// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

/* 
   Authors
   
   John Forrest

 */
#ifndef OsiFactorization_H
#define OsiFactorization_H

#include <iostream>
#include <string>

class OsiPackedMatrix;
class OsiIndexedVector;
/** This deals with Factorization and Updates

    This class started with a parallel simplex code I was writing in the
    mid 90's.  The need for parallelism led to many complications and
    I have simplified as much as I could to get back to this.

    I was aiming at problems where I might get speed-up so I was looking at dense
    problems or ones with structure.  This led to permuting input and output
    vectors and to increasing the number of rows each rank-one update.  This is 
    still in as a minor overhead.

    I have also put in handling for hyper-sparsity.  I have taken out
    all outer loop unrolling, dense matrix handling and most of the
    book-keeping for slacks.  Also I always use FTRAN approach to updating
    even if factorization fairly dense.  All these could improve performance.

    I will know if anyone is interested in developing this if she or he
    does a PFI update.  I don't know how to :-)

    I blame some of the coding peculiarities on the history of the code
    but mostly it is just because I can't do elegant code (or useful
    comments).

    I am assuming that 32 bits is enough for number of rows or columns, but OsiBigIndex
    may be redefined to get 64 bits.
 */

typedef int OsiBigIndex;


class OsiFactorization {
   friend void OsiFactorizationUnitTest( const std::string & mpsDir );

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    OsiFactorization (  );
  /// Copy constructor 
  OsiFactorization ( const OsiFactorization &other);

  /// Destructor
   ~OsiFactorization (  );
  /// Debug show object (shows one representation)
  void show_self (  ) const;
  /// = copy
    OsiFactorization & operator = ( const OsiFactorization & other );
  //@}

  /**@name Do factorization */
  //@{
  /** When part of LP - given by basic variables.
  Actually does factorization.
  Arrays passed in have non negative value to say basic.
  If status is okay, basic variables have pivot row - this is only needed
  if increasingRows_ >1.
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  int factorize ( const OsiPackedMatrix & matrix, 
		  int rowIsBasic[], int columnIsBasic[] , 
		  double areaFactor = 0.0 );
  /// Same but with OsiPackedMatrix components split out
  int factorize ( int numberRows, int numberColumns,
		  const OsiBigIndex * columnStart,
		  const int * columnLength,
		  const int * row,
		  const double * element,
		  int rowIsBasic[], int columnIsBasic[] , 
		  double areaFactor = 0.0 );
  /** When given as triplets.
  Actually does factorization.  maximumL is guessed maximum size of L part of
  final factorization, maximumU of U part.  These are multiplied by
  areaFactor which can be computed by user or internally.  
  Arrays are copied in.  I could add flag to delete arrays to save a 
  bit of memory.
  If status okay, permutation has pivot rows - this is only needed
  if increasingRows_ >1 (if not then 0,1,2..
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -99 memory */
  int factorize ( int numberRows,
		  int numberColumns,
		  OsiBigIndex numberElements,
		  OsiBigIndex maximumL,
		  OsiBigIndex maximumU,
		  const int indicesRow[],
		  const int indicesColumn[], const double elements[] ,
		  int permutation[],
		  double areaFactor = 0.0);
  /** Two part version for maximum flexibility
      This part creates arrays for user to fill.
      estimateNumberElements is safe estimate of number
      returns 0 -okay, -99 memory */
  int factorizePart1 ( int numberRows,
		       int numberColumns,
		       OsiBigIndex estimateNumberElements,
		       int * indicesRow[],
		       int * indicesColumn[],
		       double * elements[],
		  double areaFactor = 0.0);
  /** This is part two of factorization
      Arrays belong to factorization and were returned by part 1
      If status okay, permutation has pivot rows - this is only needed
      if increasingRows_ >1 (if not then 0,1,2..
      If status is singular, then basic variables have pivot row
      and ones thrown out have -1
      returns 0 -okay, -1 singular, -99 memory */
  int factorizePart2 (int permutation[],int exactNumberElements);
  //@}

  /**@name general stuff such as permutation or status */
  //@{ 
  /// Returns status
  inline int status (  ) const {
    return status_;
  };
  /// Returns number of pivots since factorization
  inline int pivots (  ) const {
    return numberPivots_;
  };
  /// Returns address of permute region
  inline int *permute (  ) const {
    return permute_;
  };
  /// Returns address of pivotColumn region (also used for permuting)
  inline int *pivotColumn (  ) const {
    return pivotColumn_;
  };
  /// Returns address of permuteBack region
  inline int *permuteBack (  ) const {
    return permuteBack_;
  };
  /// Returns address of pivotColumnBack region (also used for permuting)
  inline int *pivotColumnBack (  ) const {
    return pivotColumnBack_;
  };
  /// Number of Rows after iterating
  inline int numberRowsExtra (  ) const {
    return numberRowsExtra_;
  };
  /// Number of Rows after factorization
  inline int numberRows (  ) const {
    return numberRows_;
  };
  /// Maximum of Rows after iterating
  inline int maximumRowsExtra (  ) const {
    return maximumRowsExtra_;
  };
  /// Total number of columns in factorization
  inline int numberColumns (  ) const {
    return numberColumns_;
  };
  /// Length of FT vector
  inline int numberForrestTomlin (  ) const {
    return numberInColumn_[numberColumnsExtra_];
  };
  /// Number of good columns in factorization
  inline int numberGoodColumns (  ) const {
    return numberGoodU_;
  };
  /// Whether larger areas needed
  inline double areaFactor (  ) const {
    return areaFactor_;
  };
  inline void areaFactor ( double value ) {
    areaFactor_=value;
  };
  /// Whether rows increase after pivoting
  inline bool increasingRows (  ) const {
    return increasingRows_ > 1;
  };
  /** 0 - no increasing rows - no nothing (not coded)
      1 - no permutation (i.e. basis order in is pivot order), 
      2 user wants slacks pivoting on own rows,
      3 user needs to know everything as row are really increasing */
  inline void increasingRows ( int value  ) {
    increasingRows_ = value;
  };
  /// Level of detail of messages
  inline int messageLevel (  ) const {
    return messageLevel_ ;
  };
  void messageLevel (  int value );
  /// Maximum number of pivots between factorizations
  inline int maximumPivots (  ) const {
    return maximumPivots_ ;
  };
  void maximumPivots (  int value );
  /// Pivot tolerance
  inline double pivotTolerance (  ) const {
    return pivotTolerance_ ;
  };
  void pivotTolerance (  double value );
  /// Zero tolerance
  inline double zeroTolerance (  ) const {
    return zeroTolerance_ ;
  };
  void zeroTolerance (  double value );
  /// Whether slack value is +1 or -1
  inline double slackValue (  ) const {
    return slackValue_ ;
  };
  void slackValue (  double value );
  /// Returns maximum absolute value in factorization
  double maximumCoefficient() const;
  //@}

  /**@name some simple stuff */
  //@{

  /// Returns number in U area
  inline OsiBigIndex numberElementsU (  ) const {
    return lengthU_;
  };
  /// Returns length of U area
  inline OsiBigIndex lengthAreaU (  ) const {
    return lengthAreaU_;
  };
  /// Returns number in L area
  inline OsiBigIndex numberElementsL (  ) const {
    return lengthL_;
  };
  /// Returns length of L area
  inline OsiBigIndex lengthAreaL (  ) const {
    return lengthAreaL_;
  };
  //@}

  /**@name rank one updates which do exist */
  //@{

  /** Replaces one Column in basis,
      If checkBeforeModifying is true will do all accuracy checks
      before modifying factorization.  Whether to set this depends on
      speed considerations.  You could just do this on first iteration
      after factorization and thereafter re-factorize
      returns 0=OK, 1=Probably OK, 2=singular, 3=no room */
  int replaceColumn ( int pivotRow,
		      double pivotCheck,
		      int numberElements,
		      int index[], double array[] ,
		      bool checkBeforeModifying=false);

  /** Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
      If checkBeforeModifying is true will do all accuracy checks
      before modifying factorization.  Whether to set this depends on
      speed considerations.  You could just do this on first iteration
      after factorization and thereafter re-factorize
   partial update already in U */
  int replaceColumn ( OsiIndexedVector * regionSparse,
		      int pivotRow,
		      double pivotCheck ,
		      bool checkBeforeModifying=false);
  /// Throws away incoming column (I am not sure really needed)
  void throwAwayColumn (  );
  //@}

  /**@name various uses of factorization (return code number elements) 
   which user may want to know about */
  //@{
  /** Updates one column (FTRAN) from region2
      number returned is negative if no room
      region1 starts as zero and is zero at end */
  int updateColumn ( OsiIndexedVector * regionSparse,
			OsiIndexedVector * regionSparse2,
			bool FTUpdate = false ) ;
  /** Updates one column (FTRAN) to/from array 
      number returned is negative if no room
      indices contains list of nonzeros and is updated
      region starts as zero and is zero at end */
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
			int number,
			bool FTUpdate = false ) ;
  /** Updates one column (FTRAN) to/from array 
      number returned is negative if no room
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. rhs)
      region starts as zero and is zero at end */
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			bool FTUpdate = false ) ;
  /** These two versions has same effect as above with FTUpdate==false
      so number returned is always >=0 */
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
			int number) const;
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[] ) const;
  /** Updates one column transpose (BTRAN)
      indices contains list of nonzeros and is updated.
      returns number of nonzeros */
  int updateColumnTranspose ( OsiIndexedVector * regionSparse,
				 double array[], //unpacked
				 int index[],
				 int number) const;
  /** Updates one column transpose (BTRAN)
      ** For large problems you should ALWAYS know where the nonzeros
      are, so please try and migrate to previous method after you
      have got code working using this simple method - thank you!
      (the only exception is if you know input is dense e.g. dense objective)
      returns number of nonzeros */
  int updateColumnTranspose ( OsiIndexedVector * regionSparse,
				 double array[] ) const;
  /** Updates one column (BTRAN) from region2
      region1 starts as zero and is zero at end */
  int updateColumnTranspose ( OsiIndexedVector * regionSparse,
			      OsiIndexedVector * regionSparse2) const;
  /** makes a row copy of L for speed and to allow very sparse problems */
  void goSparse();
  /**  get sparse threshold */
  int sparseThreshold ( ) const;
  /**  set sparse threshold */
  void sparseThreshold ( int value );
  //@}
  /// *** Below this user may not want to know about

  /**@name various uses of factorization (return code number elements) 
   which user may not want to know about (left over from my LP code) */
  //@{
  /// Updates one column (FTRAN) already permuted
  int updateColumn ( OsiIndexedVector * regionSparse,
			bool FTUpdate = false ) ;
  /// This version has same effect as above with FTUpdate==false
  int updateColumn ( OsiIndexedVector * regionSparse) const;

  /// Updates one column transpose (BTRAN) already permuted
  int updateColumnTranspose ( OsiIndexedVector * regionSparse ) const;

  /** Takes off permutation vector (only needed if increasingRows_>1),
      zeroes out region2 if wanted. */
  void unPermuteTranspose ( OsiIndexedVector * regionSparse,
			    OsiIndexedVector * regionSparse2,
			    bool erase = true ) const;

  //@}

  /**@name various updates - none of which have been written! */
  //@{

  /** Adds given elements to Basis and updates factorization,
      can increase size of basis. Returns rank */
  int add ( OsiBigIndex numberElements,
	       int indicesRow[],
	       int indicesColumn[], double elements[] );

  /** Adds one Column to basis,
      can increase size of basis. Returns rank */
  int addColumn ( OsiBigIndex numberElements,
		     int indicesRow[], double elements[] );

  /** Adds one Row to basis,
      can increase size of basis. Returns rank */
  int addRow ( OsiBigIndex numberElements,
		  int indicesColumn[], double elements[] );

  /// Deletes one Column from basis, returns rank
  int deleteColumn ( int Row );
  /// Deletes one Row from basis, returns rank
  int deleteRow ( int Row );

  /** Replaces one Row in basis,
      returns 0=OK, 1=Probably OK, 2=singular */
  int replaceRow ( OsiBigIndex numberElements,
		      int indicesColumn[], double elements[] );
  //@}

private:

  /**@name used by factorization */
  /// Gets space for a factorization, called by constructors
  void getAreas ( int numberRows,
		  int numberColumns,
		  OsiBigIndex maximumL,
		  OsiBigIndex maximumU );

  /** PreProcesses raw triplet data.
      state is 0 - triplets, 1 - some counts etc , 2 - .. */
  void preProcess ( int state,
		    int possibleDuplicates = -1 );
  /// Does most of factorization
  int factor (  );
  /// Does sparse phase of factorization
  /// return code is <0 error, 0= finished
  int factorSparse (  );

  /// Pivots when just one other row so faster?
  bool pivotOneOtherRow ( int pivotRow,
			  int pivotColumn );
  /// Does one pivot on Row Singleton in factorization
  bool pivotRowSingleton ( int pivotRow,
			   int pivotColumn );
  /// Does one pivot on Column Singleton in factorization
  bool pivotColumnSingleton ( int pivotRow,
			      int pivotColumn );

  /** Gets space for one Column with given length,
   may have to do compression  (returns True if successful),
   also moves existing vector,
   extraNeeded is over and above present */
  bool getColumnSpace ( int iColumn,
			int extraNeeded );

  /** Gets space for one Row with given length,
  may have to do compression  (returns True if successful),
  also moves existing vector */
  bool getRowSpace ( int iRow, int extraNeeded );

  /** Gets space for one Row with given length while iterating,
  may have to do compression  (returns True if successful),
  also moves existing vector */
  bool getRowSpaceIterate ( int iRow,
			    int extraNeeded );
  /// Checks that row and column copies look OK
  void checkConsistency (  );
  /// Adds a link in chain of equal counts
  inline void addLink ( int index, int count ) {
    int next = firstCount_[count];
      lastCount_[index] = -2 - count;
    if ( next < 0 ) {
      //first with that count
      firstCount_[count] = index;
      nextCount_[index] = -1;
    } else {
      firstCount_[count] = index;
      nextCount_[index] = next;
      lastCount_[next] = index;
  }};
  /// Deletes a link in chain of equal counts
  inline void deleteLink ( int index ) {
    int next = nextCount_[index];
    int last = lastCount_[index];
    if ( last >= 0 ) {
      nextCount_[last] = next;
    } else {
      int count = -last - 2;

      firstCount_[count] = next;
    }				
    if ( next >= 0 ) {
      lastCount_[next] = last;
    }				
    nextCount_[index] = -2;
    lastCount_[index] = -2;
    return;
  };
  /// Cleans up at end of factorization
  void cleanup (  );

  /// Updates part of column (FTRANL)
  void updateColumnL ( OsiIndexedVector * region ) const;

  /// Updates part of column (FTRANR)
  void updateColumnR ( OsiIndexedVector * region ) const;

  /// Updates part of column (FTRANU)
  void updateColumnU ( OsiIndexedVector * region, int * indices,
			 int numberIn ) const;

  /// Updates part of column (FTRANU) when sparse
  void updateColumnUSparse ( OsiIndexedVector * regionSparse,
			     int * indices,
			     int numberIn ) const;
  /// Updates part of column (FTRANU)
  void updateColumnUDensish ( OsiIndexedVector * regionSparse) const;

  /** Updates part of column transpose (BTRANU),
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeU ( OsiIndexedVector * region) const;

  /// Updates part of column transpose (BTRANR)
  void updateColumnTransposeR ( OsiIndexedVector * region ) const;

  /// Updates part of column transpose (BTRANL)
  void updateColumnTransposeL ( OsiIndexedVector * region ) const;

  /** Returns accuracy status of replaceColumn
      returns 0=OK, 1=Probably OK, 2=singular */
  int checkPivot(double saveFromU, double oldPivot) const;
  /// The real work of constructors etc
  void gutsOfDestructor();
  /// 1 bit - tolerances etc, 2 rest (could be more granular)
  void gutsOfInitialize(int type);
  void gutsOfCopy(const OsiFactorization &other);


  /********************************* START LARGE TEMPLATE ********/
// this should have defines so will work in 64 bit mode
#define OSIFACTORIZATION_BITS_PER_INT 32
#define OSIFACTORIZATION_SHIFT_PER_INT 5
#define OSIFACTORIZATION_MASK_PER_INT 0x1f
  template <class T>  bool
  pivot ( int pivotRow,
	  int pivotColumn,
	  int pivotRowPosition,
	  int pivotColumnPosition,
	  double work[],
	  unsigned int workArea2[],
	  int increment,
	  int increment2,
	  T markRow[] ,
	  int largeInteger)
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
    std::cout << "code pivot 1" << std::endl;
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
  markRow[pivotRow] = largeInteger;
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
      if ( mark != largeInteger ) {
	//will be updated
	work[mark] = value;
	int word = mark >> OSIFACTORIZATION_SHIFT_PER_INT;
	int bit = mark & OSIFACTORIZATION_MASK_PER_INT;

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
      } else if ( mark != largeInteger ) {
	//will be updated
	work[mark] = value;;
	int word = mark >> OSIFACTORIZATION_SHIFT_PER_INT;
	int bit = mark & OSIFACTORIZATION_MASK_PER_INT;

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
	int word = j >> OSIFACTORIZATION_SHIFT_PER_INT;
	int bit = j & OSIFACTORIZATION_MASK_PER_INT;

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
	  int word = j >> OSIFACTORIZATION_SHIFT_PER_INT;
	  int bit = j & OSIFACTORIZATION_MASK_PER_INT;

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
  int bigLoops = numberInPivotColumn >> OSIFACTORIZATION_SHIFT_PER_INT;
  int i = 0;

  // do linked lists and update counts
  while ( bigLoops ) {
    bigLoops--;
    int bit;
    for ( bit = 0; bit < OSIFACTORIZATION_BITS_PER_INT; i++, bit++ ) {
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

  /********************************* END LARGE TEMPLATE ********/
  //@}
////////////////// data //////////////////
private:

  /**@name data */
  //@{
  /// Pivot tolerance
  double pivotTolerance_;
  /// Zero tolerance
  double zeroTolerance_;
  /// Whether slack value is  +1 or -1
  double slackValue_;
  /// How much to multiply areas by
  double areaFactor_;
  /// Number of Rows in factorization
  int numberRows_;
  /// Number of Rows after iterating
  int numberRowsExtra_;
  /// Maximum number of Rows after iterating
  int maximumRowsExtra_;
  /// Number of Columns in factorization
  int numberColumns_;
  /// Number of Columns after iterating
  int numberColumnsExtra_;
  /// Maximum number of Columns after iterating
  int maximumColumnsExtra_;
  /// Number factorized in U (not row singletons)
  int numberGoodU_;
  /// Number factorized in L
  int numberGoodL_;
  /// Maximum number of pivots before factorization
  int maximumPivots_;
  /// Number pivots since last factorization
  int numberPivots_;
  /// Number of elements in U (to go)
  ///       or while iterating total overall
  OsiBigIndex totalElements_;
  /// Number of elements after factorization
  OsiBigIndex factorElements_;
  /// Pivot order for each Column
  int *pivotColumn_;
  /// Permutation vector for pivot row order
  int *permute_;
  /// DePermutation vector for pivot row order
  int *permuteBack_;
  /// Inverse Pivot order for each Column
  int *pivotColumnBack_;
  /// Status of factorization
  int status_;

  /** 0 - no increasing rows - no permutations,
   1 - no increasing rows but permutations 
   2 - increasing rows */
  int increasingRows_;

  /// Detail in messages
  int messageLevel_;

  /// Number of trials before rejection
  int numberTrials_;
  /// Start of each Row as pointer
  OsiBigIndex *startRowU_;

  /// Number in each Row
  int *numberInRow_;

  /// Number in each Column
  int *numberInColumn_;

  /// Number in each Column including pivoted
  int *numberInColumnPlus_;

  /** First Row/Column with count of k,
      can tell which by offset - Rows then Columns */
  int *firstCount_;

  /// Next Row/Column with count
  int *nextCount_;

  /// Previous Row/Column with count
  int *lastCount_;

  /// Next Column in memory order
  int *nextColumn_;

  /// Previous Column in memory order
  int *lastColumn_;

  /// Next Row in memory order
  int *nextRow_;

  /// Previous Row in memory order
  int *lastRow_;

  /// Columns left to do in a single pivot
  int *saveColumn_;

  /// Marks rows to be updated
  int *markRow_;

  /// Larger of row and column size
  int biggerDimension_;

  /// Base address for U (may change)
  int *indexColumnU_;

  /// Pivots for L
  int *pivotRowL_;

  /// Inverses of pivot values
  double *pivotRegion_;

  /// Number of slacks at beginning of U
  int numberSlacks_;

  /// Number in U
  int numberU_;

  /// Base of U is always 0
  //int baseU_;

  /// Length of U
  OsiBigIndex lengthU_;

  /// Length of area reserved for U
  OsiBigIndex lengthAreaU_;

/// Elements of U
  double *elementU_;

/// Row indices of U
  int *indexRowU_;

/// Start of each column in U
  OsiBigIndex *startColumnU_;

/// Converts rows to columns in U 
  OsiBigIndex *convertRowToColumnU_;

  /// Number in L
  int numberL_;

/// Base of L
  int baseL_;

  /// Length of L
  OsiBigIndex lengthL_;

  /// Length of area reserved for L
  OsiBigIndex lengthAreaL_;

  /// Elements of L
  double *elementL_;

  /// Row indices of L
  int *indexRowL_;

  /// Start of each column in L
  OsiBigIndex *startColumnL_;

  /// Number in R
  int numberR_;

  /// Length of R stuff
  OsiBigIndex lengthR_;

  /// length of area reserved for R
  OsiBigIndex lengthAreaR_;

  /// Elements of R
  double *elementR_;

  /// Row indices for R
  int *indexRowR_;

  /// Start of columns for R
  OsiBigIndex *startColumnR_;

  /// Number of compressions doen
  OsiBigIndex numberCompressions_;

  /// true if Forrest Tomlin update, false if PFI (dummy)
  bool doForrestTomlin_;

  /// Below this use sparse technology - if 0 then no L row copy
  int sparseThreshold_;

  /// Start of each row in L
  OsiBigIndex * startRowL_;

  /// Index of column in row for L
  int * indexColumnL_;

  /// Elements in L (row copy)
  double * elementByRowL_;

  /// Sparse regions
  mutable int * sparse_;

  //@}
};
#endif
