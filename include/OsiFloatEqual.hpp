// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiFloatEqual_H
#define OsiFloatEqual_H

//#include <stdlib.h>
//#include <stdio.h>
#include <algorithm>
//#include <utilities>
//#include <limits>
//#include <numeric>
#include <cmath>



//------------------------------------------------------
//
// Function objects for testing equality of floats
//
//------------------------------------------------------

/** Function object for testing doubles for equality.

Operands are considered equal if their difference is 
within an epsilon.
This test does not consider the relative magnitude of 
the operands.<br>
Example usage:
<pre>
   double d1=3.14159;
   double d2=d1;
   double d3=d1+.0001;
   OsiAbsFltEq eq1;
   OsiAbsFltEq eq2(.001);
   assert(  eq1(d1,d2) );
   assert( !eq1(d1,d3) );
   assert(  eq2(d1,d3) );
</pre>
*/
class OsiAbsFltEq
{
public:
  /// Compare function
  inline bool operator()(double f1, double f2) const
  { return fabs(f1-f2) < epsilon_; };

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  OsiAbsFltEq():epsilon_(1.e-10){}; 

  /// Alternate constructor. Allows setting epsilon
  OsiAbsFltEq(double epsilon):epsilon_(epsilon){}; 

  /// Destructor
  virtual ~OsiAbsFltEq(){};   

  /// Copy constructor
  OsiAbsFltEq(const OsiAbsFltEq & src):epsilon_(src.epsilon_){}; 

  /// Assignment
  OsiAbsFltEq & operator=(const OsiAbsFltEq& rhs)
  {
    if (this != &rhs) {	
      epsilon_ = rhs.epsilon_;
    }
    return *this;
  };
  //@}
private:  
  /**@name Private member data */
  //@{
  /// Maximum allowed difference between operands for them to be considered equal
  double epsilon_;
  //@}
};


/** Function object for testing doubles for eqaulity.

This test considers the relative magnitude of the operands.
Two operands f1 & f2 are considered equal if<br>
abs(f1-f2) &lt;= epsilon*(max(abs(f1),abs(f2))+1)<br>
Example usage:<br>
<pre>
   double d1=3.14159;
   double d2=d1;
   double d3=d1+.0001;
   OsiRelFltEq eq1;
   assert(  eq1(d1,d2) );
   assert( !eq1(d1,d3) );
</pre>
*/
class OsiRelFltEq
{
public:
  /// Compare function
  inline bool operator()(double f1, double f2) const
  { 
    //double tol = std::max(std::abs(f1),std::abs(f2)) + 1.;
    //return std::abs(f1-f2) <= epsilon_*tol;
    double tol = fabs(f1) > fabs(f2) ? fabs(f1) : fabs(f2);
    tol++;
    return fabs(f1-f2) <= epsilon_*tol; 
  };

  /**@name Constructors and destructors */
  //@{
  /// Default constructor
  OsiRelFltEq():epsilon_(1.e-10){}; 

  /// Alternate constructor. Allows setting epsilon.
  OsiRelFltEq(double epsilon):epsilon_(epsilon){}; 

  /// Destructor
  virtual ~OsiRelFltEq(){};   

  /// Copy constructor
  OsiRelFltEq(const OsiRelFltEq & src):epsilon_(src.epsilon_){}; 

  /// Assignment
  OsiRelFltEq & operator=(const OsiRelFltEq& rhs)
  {
    if (this != &rhs) {	
      epsilon_ = rhs.epsilon_;
    }
    return *this;
  };
  //@}

private: 
  /**@name Private member data */
  //@{
  /// Epsilon
  double epsilon_;
  //@}
};

#endif
