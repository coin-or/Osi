// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cfloat>

#include "OsiRowCut.hpp"
#include "OsiFloatEqual.hpp"


//-------------------------------------------------------------------
// Row sense, rhs, range
//-------------------------------------------------------------------
char OsiRowCut::sense() const
{
  if      ( lb_ == ub_ )                        return 'E';
  else if ( lb_ == -DBL_MAX && ub_ == DBL_MAX ) return 'N';
  else if ( lb_ == -DBL_MAX )                   return 'L';
  else if ( ub_ == DBL_MAX )                    return 'G';
  else                                          return 'R';
}

double OsiRowCut::rhs() const
{
  if      ( lb_ == ub_ )                        return ub_;
  else if ( lb_ == -DBL_MAX && ub_ == DBL_MAX ) return 0.0;
  else if ( lb_ == -DBL_MAX )                   return ub_;
  else if ( ub_ == DBL_MAX )                    return lb_;
  else                                          return ub_;
}

double OsiRowCut::range() const
{
  if      ( lb_ == ub_ )                        return 0.0;
  else if ( lb_ == -DBL_MAX && ub_ == DBL_MAX ) return 0.0;
  else if ( lb_ == -DBL_MAX )                   return 0.0;
  else if ( ub_ == DBL_MAX )                    return 0.0;
  else                                          return ub_ - lb_;
}

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiRowCut::OsiRowCut ()
:
OsiCut(),
row_(),
lb_(-/*std::numeric_limits<double>::max()*/DBL_MAX),
ub_( /*std::numeric_limits<double>::max()*/DBL_MAX)
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiRowCut::OsiRowCut (
                  const OsiRowCut & source)
:
OsiCut(source),
row_(source.row_),
lb_(source.lb_),
ub_(source.ub_)
{  
  // Nothing to do here
}


//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiRowCut * OsiRowCut::clone() const
{  return (new OsiRowCut(*this));}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiRowCut::~OsiRowCut ()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiRowCut &
OsiRowCut::operator=(
                      const OsiRowCut& rhs)
{
  if (this != &rhs) {
    OsiCut::operator=(rhs);
    row_=rhs.row_;
    lb_=rhs.lb_;
    ub_=rhs.ub_;
  }
  return *this;
}
