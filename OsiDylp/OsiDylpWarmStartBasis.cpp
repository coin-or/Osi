/*! \legal
  Copyright (C) 2002, 2003
  Lou Hafer, International Business Machines Corporation and others.
  All Rights Reserved.
*/

#ifdef COIN_USE_DYLP

#ifdef _MSC_VER

/* Turn off compiler warning about long names */

#  pragma warning(disable:4786)

#endif // _MSC_VER

/* Cut name lengths for readability. */

#define ODWSB OsiDylpWarmStartBasis
#define CWSB CoinWarmStartBasis

/*!
   \file OsiDylpWarmStartBasis.cpp

   \brief Implementation of warm start object for dylp.

  This file contains the implementation of a warm start object for dylp. It
  extends the CoinWarmStartBasis class by adding an explicit list of active
  constraints. When operating in dynamic mode, dylp needs to be told which
  constraints are active in the basis.
  
  For ease of use, constraint status is handled like variable status. There's
  an array, one entry per constraint, coded using
  CoinWarmStartBasis::Status.  Inactive constraints are marked with isFree,
  active constraints with atLowerBound. Note that active/inactive is not
  equivalent to tight/loose.  The default behaviour in dynamic mode is to
  purge constraints which are strictly loose, but dylp can be instructed to
  also purge tight constraints when the associated dual variable is zero.
*/

namespace {
  char sccsid[] = "%W%	%G%" ;
  char cvsid[] = "$Id$" ;
}

#include <vector>
#include <cassert>
#include <iostream>
#include "OsiDylpWarmStartBasis.hpp"

using std::vector ;

/*
  A macro to standardize the size of the status arrays. CWSB::Status is the
  status enum, and currently it occupies 2 bits, packed 4 per char (note the
  implicit assumption that a char is 1 byte).  This macro calculates the
  number of bytes needed to hold ns status entires, rounded up to the nearest
  int.
*/

#define STATPERBYTE 4
#define STATBYTES(zz_ns_zz) \
  (((zz_ns_zz+sizeof(int)*STATPERBYTE-1)/(sizeof(int)*STATPERBYTE))* \
   sizeof(int))
  

/*! \defgroup ODWSBConstructorsDestructors ODWSB Constructors and Destructors
    \brief ODWSB Constructors and destructors

  ODWSB provides a default constructor, several variations on a copy
  constructor, a default destructor, and two variations on assignment.
  There are also routines to change the capacity of the warm start object.
*/
//@{

/*! Creates an empty basis. */

ODWSB::OsiDylpWarmStartBasis ()

  : CWSB(),
    phase_(dyINV),
    constraintStatus_(0)

{ /* intentionally left blank */ }


/*!
  A true copy --- no data structures are shared with the original.
*/

ODWSB::OsiDylpWarmStartBasis (const OsiDylpWarmStartBasis &ws)

  : CWSB(ws),
    phase_(ws.phase_),
    constraintStatus_(0)

{ int constatsze = STATBYTES(getNumArtificial()) ;
  constraintStatus_ = new char[constatsze] ;
  memcpy(constraintStatus_,ws.constraintStatus_,constatsze) ; }


/*!
  `Virtual constructor' method, for when we need to duplicate the ODWSB object
  from code that doesn't know it has an ODSWB object.
*/

CoinWarmStart *ODWSB::clone () const

{ const ODWSB *odwsb_orig = dynamic_cast<const ODWSB *>(this) ;
  ODWSB *odwsb_new = 0 ;
  if (odwsb_orig)
  { odwsb_new = new OsiDylpWarmStartBasis(*odwsb_orig) ; }
  return (dynamic_cast<CoinWarmStart *>(odwsb_new)) ; }

/*!
  Construct by copying a set of status arrays (parameters preserved). If
  no constraint status is provided, the default is to declare all constraints
  active. The phase is set to dyPRIMAL1, which is safe for any basis.

  Consider assignBasisStatus(int,int,char*&,char*&) if the object should
  assume ownership.
*/

ODWSB::OsiDylpWarmStartBasis
  (int ns, int na, const char *sStat, const char *aStat, const char *cStat)

  : CWSB(ns,na,sStat,aStat),
    phase_(dyPRIMAL1),
    constraintStatus_(0)

{ int constatsze = STATBYTES(na) ;
  constraintStatus_ = new char[constatsze] ;
/*
  If a status array is given, copy it in. If there's no array provided, assume
  all constraints are active. Set up a byte with the correct pattern and use it
  to initialize the constraint status.
*/
  if (cStat)
  { memcpy(constraintStatus_,cStat,constatsze) ; }
  else
  { char byteActive = 0 ;
    int i ;
    for (i = 0 ; i <= 3 ; i++) setStatus(&byteActive,i,CWSB::atLowerBound) ;
    memset(constraintStatus_,byteActive,constatsze) ; } }

/*!
  Assignment of structure (parameter preserved)
*/

OsiDylpWarmStartBasis& ODWSB::operator= (const OsiDylpWarmStartBasis &rhs)

{ if (this != &rhs)
  { CWSB::operator=(rhs) ;
    phase_ = rhs.phase_ ;
    delete[] constraintStatus_ ;
    int constatsze = STATBYTES(getNumArtificial()) ;
    constraintStatus_ = new char[constatsze] ;
    memcpy(constraintStatus_,rhs.constraintStatus_,constatsze) ; }
  
  return *this ; }


/*!
  Assignment of status arrays (parameters destroyed). The phase is set to
  dyPRIMAL1, which is safe for any basis.

  In this method the OsiDylpWarmStartBasis object assumes ownership of the
  arrays and upon return the argument pointers will be NULL.
  If copying is desirable, use the
  \link OsiDylpWarmStartBasis(int,int,const char*,const char*,const char*)
	array constructor \endlink
  or the
  \link operator=(const OsiDylpWarmStartBasis&)
	assignment operator \endlink.

  \note
  The pointers passed to this method will be
  freed using delete[], so they must be created using new[].
*/

void ODWSB::assignBasisStatus
  (int ns, int na, char *&sStat, char *&aStat, char *&cStat)

{ CWSB::assignBasisStatus(ns,na,sStat,aStat) ;
  phase_ = dyPRIMAL1 ;
  delete[] constraintStatus_ ;
  constraintStatus_ = cStat ;
  cStat = 0 ; }

/*!
  Assignment of status arrays (parameters destroyed). When no constraint status
  is provided, the default is to create an array which indicates all
  constraints are active. The phase is set to dyPRIMAL1, which is safe for
  any basis.

  In this method the OsiDylpWarmStartBasis object assumes ownership of the
  arrays and upon return the argument pointers will be NULL.
  If copying is desirable, use the
  \link OsiDylpWarmStartBasis(int,int,const char*,const char*,const char*)
	array constructor \endlink
  or the
  \link operator=(const OsiDylpWarmStartBasis&)
	assignment operator \endlink.

  \note
  The pointers passed to this method will be
  freed using delete[], so they must be created using new[].
*/
void ODWSB::assignBasisStatus (int ns, int na, char *&sStat, char *&aStat)

{ int constatsze = STATBYTES(na) ;
  char byteActive = 0 ;
  int i ;

  CWSB::assignBasisStatus(ns,na,sStat,aStat) ;
  phase_ = dyPRIMAL1 ;

  delete[] constraintStatus_ ;
  constraintStatus_ = new char[constatsze] ;
  for (i = 0 ; i <= 3 ; i++) setStatus(&byteActive,i,CWSB::atLowerBound) ;
  memset(constraintStatus_,byteActive,constatsze) ; }

ODWSB::~OsiDylpWarmStartBasis () { delete[] constraintStatus_ ; }

/*!
  This routine sets the capacity of the warm start object.
  Any existing basis information is lost.
*/

void ODWSB::setSize (int ns, int na)

{ CWSB::setSize(ns,na) ;
  phase_ = dyINV ;
  delete[] constraintStatus_ ;
  int constatsze = STATBYTES(na) ;
  constraintStatus_ = new char[constatsze] ;
  assert(((int) CWSB::isFree == 0)) ;
  memset(constraintStatus_,0,constatsze) ; }

/*!
  This routine sets the capacity of the warm start object.
  Any existing basis information is retained. If the new size is smaller
  than the existing basis, the existing basis is truncated to fit.
  If the new size is larger than the existing basis, additional structural
  variables are given the status nonbasic at lower bound, additional artificial
  variables are given the status basic, and additional constraints
  are given the status active.
*/

void ODWSB::resize (int numRows, int numCols)

{ int concnt = getNumArtificial() ;
  
  CWSB::resize(numRows,numCols) ;

  if (numRows != concnt)
  { int oldsze = STATBYTES(concnt) ;
    int newsze = STATBYTES(numRows) ;
    char *newStat = new char[newsze] ;
/*
  If the new capacity is smaller, truncate the existing basis.
*/
    if (oldsze > newsze)
    { memcpy(newStat,constraintStatus_,newsze) ; }
/*
  If the new capacity is larger, copy the existing basis, then mark the
  additional constraints as active.
*/
    else
    { char byteActive = 0 ;
      int i ;
      for (i = 0 ; i <= 3 ; i++) setStatus(&byteActive,i,CWSB::atLowerBound) ;
      memcpy(newStat,constraintStatus_,oldsze*sizeof(char)) ;
      memset(newStat+oldsze,byteActive,newsze-oldsze) ; }

    delete [] constraintStatus_ ;
    constraintStatus_ = newStat ; } }


/*!
   Removal of a tight constraint with a nonbasic logical implies that some
   basic variable must be kicked out of the basis. There's no cheap, efficient
   way to choose this variable, so the choice is left to the client.
*/
/*
  To elaborate on the above comment: When a tight constraint is removed, new
  extreme points are exposed. There's a reasonable likelihood they'll be
  desireable points (since the tight constraint was most likely tight because
  it was blocking movement in a desireable direction). In essence, we want to
  take the next primal pivot and see what comes tight. That's a fair bit of
  work, and on balance best left to the client to decide if it's appropriate,
  or if something else should be done.
*/

void ODWSB::deleteRows (int number, const int *which)

{ int oldconcnt = getNumArtificial() ;
  int i,k ;

  CWSB::deleteRows(number,which) ;

  int delcnt = 0 ;
  vector<bool> deleted = vector<bool>(oldconcnt) ;
  for (k = 0 ; k < number ; k++)
  { i = which[k] ;
    if (i >= 0 && i < oldconcnt && !deleted[i])
    { delcnt++ ;
      deleted[i] = true ; } }

  int newsze = STATBYTES(oldconcnt-delcnt) ;
  char *newStat = new char[newsze] ;

  k = 0 ;
  for (i = 0 ; i < oldconcnt ; i++)
  { Status status = getConStatus(i);
    if (!deleted[i])
    { setStatus(newStat,k,status) ;
      k++ ; } }
  
  delete [] constraintStatus_ ;
  constraintStatus_ = newStat ; }

/*
  The good news is that CWSB::deleteColumns works just fine for ODWSB.
*/

//@}

/*! \defgroup BasisInfo Methods to get and set basis information
    \brief Methods to get and set basis information
*/

//@{

/*!
  Count the number of active constraints by scanning the constraint status
  array.
*/

int ODWSB::numberActiveConstraints () const

{ int i ;
  int concnt = getNumArtificial() ;
  int actcnt = 0 ;

  for (i = 0 ; i < concnt ; i++)
    if (getStatus(constraintStatus_,i) == CWSB::atLowerBound) actcnt++ ;
  
  return (actcnt) ; }


/*! \defgroup MiscUtil Miscellaneous Utilities
    \brief Miscellaneous utility functions
*/

//@{

/*!
  Print the warm start object in a (more or less) human-readable form.
  Intended for debugging.
*/

void ODWSB::print () const

{ char conlett[] = {'I','?','?','A'} ;
  char statlett[] = {'F','B','U','L'} ;
  int i ;

  std::cout << "ODWSB: " ;
  std::cout << getNumArtificial() << " constraints (" <<
	       numberActiveConstraints() << " active), " ;
  std::cout << getNumStructural() << " variables." << std::endl ;

  std::cout << "Rows: " ;
  for (i = 0 ; i < getNumArtificial() ; i++)
  { std::cout << conlett[getConStatus(i)] ; }
  std::cout << std::endl << "      " ;
  for (i = 0 ; i < getNumArtificial() ; i++)
  { std::cout << statlett[getArtifStatus(i)] ; }

  std::cout << std::endl << "Cols: " ;
  for (i = 0 ; i < getNumStructural() ; i++)
  { std::cout << statlett[getStructStatus(i)] ; }

  std::cout << std::endl ;
  
  return ; }

#endif // COIN_USE_DYLP
