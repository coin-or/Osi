// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cfloat>

#include "OsiRowCut.hpp"

#ifndef OSI_INLINE_ROWCUT_METHODS

//  //-------------------------------------------------------------------
//  // Set/Get lower & upper bounds
//  //-------------------------------------------------------------------
double OsiRowCut::lb() const { return lb_; }
void OsiRowCut::setLb(double lb) { lb_ = lb; }
double OsiRowCut::ub() const { return ub_; }
void OsiRowCut::setUb(double ub) { ub_ = ub; }

//-------------------------------------------------------------------
// Set row elements
//------------------------------------------------------------------- 
void OsiRowCut::setRow(int size, 
		       const int * colIndices, const double * elements)
{
  row_.setVector(size,colIndices,elements);
}
void OsiRowCut::setRow( const CoinPackedVector & v )
{
  row_ = v;
}

//-------------------------------------------------------------------
// Get the row
//-------------------------------------------------------------------
const CoinPackedVector & OsiRowCut::row() const 
{ 
  return row_; 
}

//----------------------------------------------------------------
// == operator 
//-------------------------------------------------------------------
bool
OsiRowCut::operator==(const OsiRowCut& rhs) const
{
  if ( this->OsiCut::operator!=(rhs) ) return false;
  if ( row() != rhs.row() ) return false;
  if ( lb() != rhs.lb() ) return false;
  if ( ub() != rhs.ub() ) return false;
  return true;
}
bool
OsiRowCut::operator!=(const OsiRowCut& rhs) const
{
  return !( (*this)==rhs );
}


//----------------------------------------------------------------
// consistent & infeasible 
//-------------------------------------------------------------------
bool OsiRowCut::consistent() const
{
  const CoinPackedVector & r=row();
  r.duplicateIndex("consistent", "OsiRowCut");
  if ( r.getMinIndex() < 0 ) return false;
  return true;
}
bool OsiRowCut::consistent(const OsiSolverInterface& im) const
{  
  const CoinPackedVector & r=row();
  if ( r.getMaxIndex() >= im.getNumCols() ) return false;

  return true;
}
bool OsiRowCut::infeasible(const OsiSolverInterface &im) const
{
  if ( lb() > ub() ) return true;

  return false;
}

#endif

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
