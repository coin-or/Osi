// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiColCut.hpp"

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiColCut::OsiColCut() :
   OsiCut(),
   lbs_(),
   ubs_()
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiColCut::OsiColCut(const OsiColCut & source) :
   OsiCut(source),
   lbs_(source.lbs_),
   ubs_(source.ubs_)
{  
  // Nothing to do here
}


//----------------------------------------------------------------
// Clone
//----------------------------------------------------------------
OsiColCut * OsiColCut::clone() const
{ return (new OsiColCut(*this)); }

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiColCut::~OsiColCut ()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiColCut &
OsiColCut::operator=(const OsiColCut& rhs)
{
  if (this != &rhs) {
    
    OsiCut::operator=(rhs);
    lbs_=rhs.lbs_;
    ubs_=rhs.ubs_;
  }
  return *this;
}
