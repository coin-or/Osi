#ifndef OsiFactorization_H
#define OsiFactorization_H

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

public:

  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
    OsiFactorization (  );
  /// Copy constructor 
  OsiFactorization ( const OsiFactorization &other);

  /// Destructor
   ~OsiFactorization (  );
  /// Debug show object
  void show_self (  ) const;
  /// = copy
    OsiFactorization & operator = ( const OsiFactorization & other );
  //@}

  /**@name Do factorization */
  //@{
  /** When part of LP - given by basic variables.
  Actually does factorization.
  Arrays passed in have non negative value to say basic.
  If status is okay, basic variables have pivot row.
  If status is singular, then basic variables have +1 and ones thrown out have -INT_MAX
  to say thrown out.
  returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  int factorize ( OsiPackedMatrix & matrix, 
		  int rowIsBasic[], int columnIsBasic[] , 
		  double areaFactor = 0.0 );

  /** When given as triplets.
  Actually does factorization.  maximumL is guessed maximum size of L part of
  final factorization, maximumU of U part.  These are multiplied by
  areaFactor which can be computed by user or internally.  
  Arrays are copied in.  I could add flag to delete arrays to save a 
  bit of memory.
  If status okay, array has pivot rows.
  If status is singular, then basic variables have +1 and ones thrown out have -INT_MAX
  to say thrown out.
  returns 0 -okay, -1 singular, -99 memory */
  int factorize ( int numberOfRows,
		  int numberOfColumns,
		  OsiBigIndex numberOfElements,
		  OsiBigIndex maximumL,
		  OsiBigIndex maximumU,
		  const int indicesRow[],
		  const int indicesColumn[], const double elements[] ,
		  int permutation[],
		  double areaFactor = 0.0);
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
  /// 0 - no permutation (not coded), 1 user need not bother,2 user needs to know
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
      returns 0=OK, 1=Probably OK, 2=singular, 3=no room */
  int replaceColumn ( int pivotRow,
			 double pivotCheck,
			 int numberOfElements,
			 int index[], double array[] );

  /** Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   partial update already in U */
  int replaceColumn ( OsiIndexedVector * regionSparse,
			 int pivotRow,
			 double pivotCheck );
  /// Throws away incoming column
  void throwAwayColumn (  );
  //@}

  /**@name various uses of factorization (return code number elements) */
  //@{
  /** Updates one column (FTRAN) from region2
      This assumes user is thinking non-permuted
      - returns un-permuted result in region2.
      region1 starts as zero and is zero at end */
  int updateColumn ( OsiIndexedVector * regionSparse,
			OsiIndexedVector * regionSparse2,
			bool FTUpdate = false ) ;
  /** Updates one column (FTRAN) to/from array
      This assumes user is thinking non-permuted
      - returns un-permuted result in array.
      number returned is negative if no room
      region starts as zero and is zero at end */
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
			int number,
			bool FTUpdate = false ) ;
  /// This version has same effect as above with FTUpdate==false
  int updateColumn ( OsiIndexedVector * regionSparse,
			double array[], //unpacked
			int index[],
			int number) const;
  /// Updates one column (FTRAN) already permuted
  int updateColumn ( OsiIndexedVector * regionSparse,
			bool FTUpdate = false ) ;
  /// This version has same effect as above with FTUpdate==false
  int updateColumn ( OsiIndexedVector * regionSparse) const;

  /// Updates one column transpose (BTRAN) already permuted
  int updateColumnTranspose ( OsiIndexedVector * regionSparse ) const;

  /// Updates one column transpose (BTRAN) - assumes user is thinking non-permuted
  int updateColumnTranspose ( OsiIndexedVector * regionSparse,
				 double array[], //unpacked
				 int index[],
				 int number) const;


  /** Takes off permutation vector (only needed if increasingRows_>1),
      zeroes out region2 if wanted. */
  void unPermuteTranspose ( OsiIndexedVector * regionSparse,
			    OsiIndexedVector * regionSparse2,
			    bool erase = true ) const;

  /** makes a row copy of L for speed */
  void makeRowCopyL();
  /**  get sparse threshold */
  int sparseThreshold ( ) const;
  /**  set sparse threshold */
  void sparseThreshold ( int value );
  //@}

  /**@name various updates - none of which have been written! */
  //@{

  /** Adds given elements to Basis and updates factorization,
      can increase size of basis. Returns rank */
  int add ( OsiBigIndex numberOfElements,
	       int indicesRow[],
	       int indicesColumn[], double elements[] );

  /** Adds one Column to basis,
      can increase size of basis. Returns rank */
  int addColumn ( OsiBigIndex numberOfElements,
		     int indicesRow[], double elements[] );

  /** Adds one Row to basis,
      can increase size of basis. Returns rank */
  int addRow ( OsiBigIndex numberOfElements,
		  int indicesColumn[], double elements[] );

  /// Deletes one Column from basis, returns rank
  int deleteColumn ( int Row );
  /// Deletes one Row from basis, returns rank
  int deleteRow ( int Row );

  /** Replaces one Row in basis,
      returns 0=OK, 1=Probably OK, 2=singular */
  int replaceRow ( OsiBigIndex numberOfElements,
		      int indicesColumn[], double elements[] );
  //@}

private:

  /**@name used by factorization */
  /// Gets space for a factorization, called by constructors
  void getAreas ( int numberOfRows,
		  int numberOfColumns,
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

  /** Does one pivot in factorization,
  pivot routines return false on lack of memory.
   This version is for large nucleus */
  bool pivot ( int pivotRow,
	       int pivotColumn,
	       int pivotRowPosition,
	       int pivotColumnPosition,
	       double work[],
	       unsigned int workArea2[],
	       int increment,
	       int increment2, int markRow[] );

  /** Does one pivot in factorization,
  pivot routines return false on lack of memory.
   This version is for small nucleus */
  bool pivot ( int pivotRow,
	       int pivotColumn,
	       int pivotRowPosition,
	       int pivotColumnPosition,
	       double work[],
	       unsigned int workArea2[],
	       int increment,
	       int increment2, short markRow[] );
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

  /// The real work of constructors etc
  void gutsOfDestructor();
  /// 1 bit - tolerances etc, 2 rest (could be more granular)
  void gutsOfInitialize(int type);
  void gutsOfCopy(const OsiFactorization &other);

  //@}
////////////////// data //////////////////
protected:

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
// this should have defines so will work in 64 bit mode
#define BITS_PER_INT 32
#define SHIFT_PER_INT 5
#define MASK_PER_INT 0x1f
#endif
