// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiSolverInterface.hpp"
#include "OsiSimplexInterface.hpp"

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiSimplexInterface::OsiSimplexInterface ()
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiSimplexInterface::OsiSimplexInterface (
                  const OsiSimplexInterface & rhs)
{  
  // nothing to do here
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiSimplexInterface::~OsiSimplexInterface ()
{
  // nothing to do here
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiSimplexInterface &
OsiSimplexInterface::operator=(const OsiSimplexInterface& rhs)
{
  if (this != &rhs) {
  }
  return *this;
}




