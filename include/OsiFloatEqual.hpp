// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiFloatEqual_H
#define OsiFloatEqual_H

#include <algorithm>
#include <cmath>

#if defined(_MSC_VER)
# include<float.h>
# define isnan   _isnan
# define finite  _finite
#endif

#ifdef sun
  extern "C" int finite(double),isnan(double) ;
#endif


/*! \file OsiFloatEqual.hpp
    \brief Function objects for testing equality of real numbers.

  Two objects are provided; one tests for equality to an absolute tolerance,
  one to a scaled tolerance. The tests will handle IEEE floating point, but
  note that infinity == infinity. Mathematicians are rolling in their graves,
  but this matches the behaviour for the common practice of using DBL_MAX
  (numeric_limits<double>::max(), or similar large finite number) as infinity.

  <p>
  Example usage:
  <pre>
    double d1 = 3.14159 ;
    double d2 = d1 ;
    double d3 = d1+.0001 ;

    OsiAbsFltEq eq1 ;
    OsiAbsFltEq eq2(.001) ;

    assert(  eq1(d1,d2) ) ;
    assert( !eq1(d1,d3) ) ;
    assert(  eq2(d1,d3) ) ;
  </pre>
  OsiRelFltEq follows the same pattern.
*/

/*! \brief Equality to an absolute tolerance

  Operands are considered equal if their difference is within an epsilon ;
  the test does not consider the relative magnitude of the operands.
*/

class OsiAbsFltEq
{
  public:

  //! Compare function

  inline bool operator() (const double f1, const double f2) const

  { if (isnan(f1) || isnan(f2)) return false ;
    if (f1 == f2) return true ;
    return (fabs(f1-f2) < epsilon_) ; } ;

  /*! \name Constructors and destructors */
  //@{

  //! Default constructor

  OsiAbsFltEq () : epsilon_(1.e-10) {} ;

  //! Alternate constructor with epsilon as a parameter

  OsiAbsFltEq (const double epsilon) : epsilon_(epsilon) {} ;

  //! Destructor

  virtual ~OsiAbsFltEq () {} ;

  //! Copy constructor

  OsiAbsFltEq (const OsiAbsFltEq& src) : epsilon_(src.epsilon_) {} ;

  //! Assignment

  OsiAbsFltEq& operator= (const OsiAbsFltEq& rhs)

  { if (this != &rhs) epsilon_ = rhs.epsilon_ ;
    return (*this) ; } ;

  //@}

  private:  

  /*! \name Private member data */
  //@{

  //! Equality tolerance.

  double epsilon_ ;

  //@}

} ;



/*! \brief Equality to a scaled tolerance

  Operands are considered equal if their difference is within a scaled
  epsilon calculated as epsilon_*(1+max(|f1|,|f2|)).
*/

class OsiRelFltEq
{
  public:

  //! Compare function

  inline bool operator() (const double f1, const double f2) const

  { if (isnan(f1) || isnan(f2)) return false ;
    if (f1 == f2) return true ;
    if (!finite(f1) || !finite(f2)) return false ;

    double tol = (fabs(f1)>fabs(f2))?fabs(f1):fabs(f2) ;

    return (fabs(f1-f2) <= epsilon_*(1+tol)) ; }

  /*! \name Constructors and destructors */
  //@{

  //! Default constructor

  OsiRelFltEq () : epsilon_(1.e-10) {} ;

  //! Alternate constructor with epsilon as a parameter

  OsiRelFltEq (const double epsilon) : epsilon_(epsilon) {} ;

  //! Destructor

  virtual ~OsiRelFltEq () {} ;

  //! Copy constructor

  OsiRelFltEq (const OsiRelFltEq & src) : epsilon_(src.epsilon_) {} ;

  //! Assignment

  OsiRelFltEq& operator= (const OsiRelFltEq& rhs)

  { if (this != &rhs) epsilon_ = rhs.epsilon_ ;
    return (*this) ; } ;

  //@}

private: 

  /*! \name Private member data */
  //@{

  //! Base equality tolerance

  double epsilon_ ;

  //@}

} ;

#endif
