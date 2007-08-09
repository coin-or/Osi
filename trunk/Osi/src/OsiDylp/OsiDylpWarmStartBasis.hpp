#ifndef OsiDylpWarmStartBasis_H
#define OsiDylpWarmStartBasis_H

/*! \legal
  Copyright (C) 2003 -- 2007
  Lou Hafer, International Business Machines Corporation and others.
  All Rights Reserved.
*/

/*! \file OsiDylpWarmStartBasis.hpp
    \brief Declaration of the warm start class for dylp.
*/

/*
  sccs: @(#)OsiDylpWarmStartBasis.hpp	1.5	09/16/04
  cvs: $Id$
*/

#include "CoinWarmStartBasis.hpp"

#define DYLP_INTERNAL
extern "C" {
#include "dylp.h"
}



/*! \class OsiDylpWarmStartBasis
    \brief The dylp warm start class

  This derived class is necessary because dylp by default works with a subset
  of the full constraint system. The warm start object needs to contain a
  list of the active constraints in addition to the status information
  included in CoinWarmStartBasis. It is also convenient to include the solver
  phase in the warm start object.
*/


class OsiDylpWarmStartBasis : public CoinWarmStartBasis

{ public:

/*! \name Methods to get and set basis information.

  Methods for structural and artificial variables are inherited from
  CoinWarmStartBasis.  Constraint status is coded using the
  CoinWarmStartBasis::Status codes. Active constraints are coded as
  atLowerBound, inactive as isFree.
*/
//@{

  /*! \brief Return the number of active constraints */

  int numberActiveConstraints() const ;


  /*! \brief Return the status of the specified constraint. */

  inline Status getConStatus (int i) const

  { const int st = (constraintStatus_[i>>2] >> ((i&3)<<1)) & 3 ;
    return (static_cast<CoinWarmStartBasis::Status>(st)) ; }


  /*! \brief Set the status of the specified constraint */

  inline void setConStatus (int i, Status st)

  { char &st_byte = constraintStatus_[i>>2] ;
    st_byte &= ~(3 << ((i&3)<<1)) ;
    st_byte |= (st << ((i&3)<<1)) ; }


  /*! \brief Return the status array for constraints. */

  inline char *getConstraintStatus () { return (constraintStatus_) ; }

  /** \c const overload for
    \link OsiDylpWarmStartBasis::getConstraintStatus()
    	  getConstraintStatus()
    \endlink
  */

  inline const char *getConstraintStatus () const
  
  { return (constraintStatus_) ; }


  /*! \brief Set the lp phase for this basis */

  inline void setPhase (dyphase_enum phase) { phase_ = phase ; }

  /*! \brief Get the lp phase for this basis. */

  inline dyphase_enum getPhase () const { return (phase_) ; }

//@}

/*! \name Basis `diff' methods */
//@{

  /*! \brief Generate a `diff' that can convert oldBasis to this basis.  */

  CoinWarmStartDiff *generateDiff (const CoinWarmStart *const oldCWS) const ;

  /*! \brief Apply \p diff to this basis */

  void applyDiff (const CoinWarmStartDiff *const cwsdDiff) ;

//@}

/*! \name Methods to modify the warm start object */
//@{

  /*! \brief Set basis capacity; existing basis is discarded */

  void setSize (int ns, int na) ;

  /*! \brief Set basis capacity; existing basis is maintained */

  void resize (int numRows, int numCols) ;

  /*! \brief Delete a set of rows from the basis

    \warning
    This routine assumes that the set of indices to be deleted is sorted
    in ascending order and is free from duplicates. Use deleteRows if
    this is not guaranteed.

    \warning
    The resulting basis is guaranteed valid only if all deleted
    constraints are slack (hence the associated logicals are basic).
  */

  void compressRows (int tgtCnt, const int *tgts) ;

  /*! \brief Delete a set of rows from the basis
  
    \warning
    The resulting basis is guaranteed valid only if all deleted
    constraints are slack (hence the associated logicals are basic).
  */

  void deleteRows (int number, const int *which) ;

  /** \brief Merge entries from a source basis into this basis.

    \warning
    It's the client's responsibility to ensure validity of the merged basis,
    if that's important to the application.

    The vector xferCols (xferRows) specifies runs of entries to be taken from
    the source basis and placed in this basis. Each entry is a CoinTriple,
    with first specifying the starting source index of a run, second
    specifying the starting destination index, and third specifying the run
    length.
  */

  virtual void mergeBasis(const CoinWarmStartBasis *src,
			  const XferVec *xferRows,
			  const XferVec *xferCols) ;

//@}

/*! \name Constructors, destructors, and related functions */
//@{

  /*! \brief Default constructor (empty object) */

  OsiDylpWarmStartBasis () ;

  /*! \brief Constructs a warm start object with the specified status arrays */

  OsiDylpWarmStartBasis (int ns, int na, const char *sStat,
			 const char *aStat, const char *cStat = 0) ;

  /*! \brief Construct an OsiDylpWarmStartBasis from a CoinWarmStartBasis */

  OsiDylpWarmStartBasis (const CoinWarmStartBasis &cwsb) ;

  /*! \brief Copy constructor */

  OsiDylpWarmStartBasis (const OsiDylpWarmStartBasis &ws) ;

  /*! \brief `Virtual constructor' */

  CoinWarmStart *clone () const ;

  /*! \brief Destructor */

  ~OsiDylpWarmStartBasis () ;

  /*! \brief Assignment */

  OsiDylpWarmStartBasis& operator= (const OsiDylpWarmStartBasis &rhs) ;

  /*! \brief Assign the status vectors to be the warm start information */

  void assignBasisStatus
    (int ns, int na, char *&sStat, char *&aStat, char *&cStat) ;

  /*! \brief Assign the status vectors to be the warm start information */

  void assignBasisStatus
    (int ns, int na, char *&sStat, char *&aStat) ;


//@}

/*! \name Miscellaneous methods */
//@{

  /*! \brief Prints in readable format (for debug) */

  void print () const ;

  /*! \brief Performs basis consistency checks (for debug) */

  void checkBasis () const ;

//@}

  private:

/*! \name Constraint status private data members */
//@{

  dyphase_enum phase_ ;		///< dylp phase

  char *constraintStatus_ ;	///< vector of constraint status information

//@}

} ;



/*! \class OsiDylpWarmStartBasisDiff
    \brief A `diff' between two OsiDylpWarmStartBasis objects

  This class exists in order to hide from the world the details of
  calculating and representing a `diff' between two OsiDylpWarmStartBasis
  objects. For convenience, assignment, cloning, and deletion are visible to
  the world, and default and copy constructors are visible to derived classes.
  Knowledge of the rest of this structure, and of generating and
  applying diffs, is restricted to the functions
  OsiDylpWarmStartBasis::generateDiff() and OsiDylpWarmStartBasis::applyDiff().

  The actual data structure is a pair of unsigned int vectors, #diffNdxs_ and
  #diffVals_, and a CoinWarmStartBasisDiff object.

  \todo This is a pretty generic structure, and vector diff is a pretty generic
	activity. We should be able to convert this to a template.

  \todo Using unsigned int as the data type for the diff vectors might help
	to contain the damage when this code is inevitably compiled for 64 bit
	architectures. But the notion of int as 4 bytes is hardwired into
	CoinWarmStartBasis, so changes are definitely required.
*/

class OsiDylpWarmStartBasisDiff : public CoinWarmStartBasisDiff
{ public:

  /*! \brief `Virtual constructor' */
  virtual CoinWarmStartDiff *clone() const
  { OsiDylpWarmStartBasisDiff *odwsbd =  new OsiDylpWarmStartBasisDiff(*this) ;
    return (dynamic_cast<CoinWarmStartDiff *>(odwsbd)) ; }

  /*! \brief Assignment */
  virtual
  OsiDylpWarmStartBasisDiff &operator= (const OsiDylpWarmStartBasisDiff &rhs) ;

  /*! \brief Destructor */
  virtual ~OsiDylpWarmStartBasisDiff()
  { delete[] condiffNdxs_ ;
    delete[] condiffVals_ ; }

  private:

  friend CoinWarmStartDiff *OsiDylpWarmStartBasis::generateDiff
    (const CoinWarmStart *const oldCWS) const ;
  friend void OsiDylpWarmStartBasis::applyDiff
    (const CoinWarmStartDiff *const diff) ;

  /*! \brief Standard constructor */
  OsiDylpWarmStartBasisDiff (int sze, const unsigned int *const diffNdxs,
			     const unsigned int *const diffVals,
			     const CoinWarmStartBasisDiff *const cwsbd) ;
  
  /*! \brief Default constructor */
  OsiDylpWarmStartBasisDiff ()
    : CoinWarmStartBasisDiff(),
      consze_(0),
      condiffNdxs_(0),
      condiffVals_(0)
  { /* intentionally left blank */ } 

  /*! \brief Copy constructor
  
    For convenience when copying objects containing OsiDylpWarmStartBasisDiff
    objects. But consider whether you should be using #clone() to retain
    polymorphism.
  */
  OsiDylpWarmStartBasisDiff (const OsiDylpWarmStartBasisDiff &odwsbd) ;
    
  /* Data members */

  /*! \brief Number of entries (and allocated capacity), in units of \c int. */
  int consze_ ;

  /*! \brief Array of diff indices for constraint status */

  unsigned int *condiffNdxs_ ;

  /*! \brief Array of diff values for constraint status */

  unsigned int *condiffVals_ ;

} ;



#endif // OsiDylpWarmStartBasis_H
