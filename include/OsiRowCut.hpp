// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiRowCut_H
#define OsiRowCut_H


#include "OsiCollections.hpp"
#include "OsiPackedVector.hpp"
#include "OsiCut.hpp"

/** Row Cut Class

A row cut has:
  <ul>
  <li>a lower bound<br>
  <li>an upper bound<br>
  <li>a vector of row elements
  </ul>
*/
class OsiRowCut : public OsiCut {
   friend void OsiRowCutUnitTest(const OsiSolverInterface * baseSiP,    
				 const std::string & mpsDir);

public:
  
  /**@name Row bounds */
  //@{
    /// Get lower bound
    inline double lb() const;
    /// Set lower bound
    inline void setLb(double lb);
    /// Get upper bound
    inline double ub() const;
    /// Set upper bound
    inline void setUb(double ub);
  //@}

  /**@name Row rhs, sense, range */
  //@{
    /// Get sense ('E', 'G', 'L', 'N', 'R')
    char sense() const;
    /// Get right-hand side
    double rhs() const;
    /// Get range (ub - lb for 'R' rows, 0 otherwise)
    double range() const;
  //@}

  //-------------------------------------------------------------------
  /**@name Row elements  */
  //@{
    /// Set row elements
    inline void setRow( 
      int size, 
      const int * colIndices, 
      const double * elements );
    /// Set row elements from a packed vector
    inline void setRow( const OsiPackedVector & v );
    /// Get row elements
    inline const OsiPackedVector & row() const;
  //@}

  /**@name Comparison operators  */
  //@{
    /** equal - true if lower bound, upper bound, row elements,
        and OsiCut are equal.
    */
    inline virtual bool operator==(const OsiRowCut& rhs) const; 
    /// not equal
    inline virtual bool operator!=(const OsiRowCut& rhs) const; 
  //@}
  
    
  //----------------------------------------------------------------
  /**@name Sanity checks on cut */
  //@{
    /** Returns true if the cut is consistent.
        This checks to ensure that:
        <ul>
        <li>The row element vector does not have duplicate indices
        <li>The row element vector indices are >= 0
        </ul>
    */
    inline virtual bool consistent() const; 

    /** Returns true if cut is consistent with respect to the solver
        interface's model.
        This checks to ensure that
        <ul>
        <li>The row element vector indices are < the number of columns
            in the model
        </ul>
    */
    inline virtual bool consistent(const OsiSolverInterface& im) const;

    /** Returns true if the row cut itself is infeasible and cannot be satisfied.       
        This checks whether
        <ul>
        <li>the lower bound is strictly greater than the
            upper bound.
        </ul>
    */
    inline virtual bool infeasible(const OsiSolverInterface &im) const;
  //@}


  /**@name Constructors and destructors */
  //@{
    /// Assignment operator
    OsiRowCut & operator=( const OsiRowCut& rhs);
  
    /// Copy constructor 
    OsiRowCut ( const OsiRowCut &);  

    /// Clone
    virtual OsiRowCut * clone() const;
  
    /// Default Constructor 
    OsiRowCut ();
  
    /// Destructor 
    virtual ~OsiRowCut ();
  //@}

private:
  
 
  /**@name Private member data */
  //@{
    /// Row elements
    OsiPackedVector row_;
    /// Row lower bound
    double lb_;
    /// Row upper bound
    double ub_;
  //@}
};


//-------------------------------------------------------------------
// Set/Get lower & upper bounds
//-------------------------------------------------------------------
double OsiRowCut::lb() const { return lb_; }
void OsiRowCut::setLb(double lb) { lb_ = lb; }
double OsiRowCut::ub() const { return ub_; }
void OsiRowCut::setUb(double ub) { ub_ = ub; }

//-------------------------------------------------------------------
// Set row elements
//------------------------------------------------------------------- 
void OsiRowCut::setRow( 
                          int size, 
                          const int * colIndices, 
                          const double * elements )
{
  row_.setVector(size,colIndices,elements);
}
void OsiRowCut::setRow( const OsiPackedVector & v )
{
  row_ = v;
}

//-------------------------------------------------------------------
// Get the row
//-------------------------------------------------------------------
const OsiPackedVector & OsiRowCut::row() const 
{ 
  return row_; 
}

//----------------------------------------------------------------
// == operator 
//-------------------------------------------------------------------
bool
OsiRowCut::operator==(
                       const OsiRowCut& rhs) const
{
  if ( this->OsiCut::operator!=(rhs) ) return false;
  if ( row() != rhs.row() ) return false;
  if ( lb() != rhs.lb() ) return false;
  if ( ub() != rhs.ub() ) return false;
  return true;
}
bool
OsiRowCut::operator!=(
                       const OsiRowCut& rhs) const
{
  return !( (*this)==rhs );
}


//----------------------------------------------------------------
// consistent & infeasible 
//-------------------------------------------------------------------
bool OsiRowCut::consistent() const
{
  const OsiPackedVector & r=row();
  r.duplicateIndex("consistent", "OsiRowCut");
  if ( r.getMinIndex() < 0 ) return false;
  return true;
}
bool OsiRowCut::consistent(const OsiSolverInterface& im) const
{  
  const OsiPackedVector & r=row();
  if ( r.getMaxIndex() >= im.getNumCols() ) return false;

  return true;
}
bool OsiRowCut::infeasible(const OsiSolverInterface &im) const
{
  if ( lb() > ub() ) return true;

  return false;
}

//#############################################################################
/** A function that tests the methods in the OsiRowCut class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
OsiRowCutUnitTest(const OsiSolverInterface * baseSiP,    
		  const std::string & mpsDir);

#endif
