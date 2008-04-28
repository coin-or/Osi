// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiRowCutDebugger_H
#define OsiRowCutDebugger_H

#include <string>

#include "OsiCuts.hpp"
#include "OsiSolverInterface.hpp"

/** Validate Row Cut Generator */
class OsiRowCutDebugger {
  friend void OsiRowCutDebuggerUnitTest(const OsiSolverInterface * siP,    
					const std::string & mpsDir);

public:
  
  /**@name Validate Row Cuts */
  //@{
  /** If we are on the path to the optimal integer solution then
      check if any generated cuts cut off the optimal solution!

      If so then print offending cuts and return non-zero code

      Up to user to check if on optimalPath (using function of same name).
      This is normally handled by rowCutDebugger() in OsiSolverInterface.

      Return number of invalid cuts.
  */
  virtual int validateCuts(const OsiCuts & cs, int first, int last) const;

  /// check one cut. Return true if cut is invalid
  virtual bool invalidCut(const OsiRowCut & rowcut) const;

  /// Return optimal solution
  inline const double * optimalSolution() const
  { return optimalSolution_;}

  /// Return number of columns in optimal solution
  inline int numberColumns() const { return (numberColumns_) ; }
  /// Return value of optimal solution
  inline double optimalValue() const
  { return optimalValue_;}
  //@}

  /**@name Activate Debugger */
  //@{
  /** Activate debugger using name of model.
      It may or may not work if problem presolved.
      Returns true if debugger activated.
  */
  bool activate(const OsiSolverInterface & si, const char * model);
  /** Activate debugger using full solution array.
      Only integer values need to be correct.
      Up to user to get it correct.
      Returns true if debugger activated (i.e. solution was valid).
  */
  bool activate(const OsiSolverInterface & si, const double * solution);
  /// Redo solution after preprocessing
  void redoSolution(int numberColumns,const int * originalColumns);
  /// Print optimal solution (returns -1 bad debug, 0 on optimal, 1 not)
  int printOptimalSolution(const OsiSolverInterface & si) const;
  //@}

  /**@name Test if on Optimal Path */
  //@{
  /** Returns whether still on optimal path.
      This should normally be invoked from 
      OsiSolverInterface::rowCutDebugger()
  */
  bool onOptimalPath(const OsiSolverInterface & si) const;
  //@}

  /**@name Test if debugger active */
  //@{
  /// Returns true if debugger is active 
  bool active() const;
  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor - no checking 
  OsiRowCutDebugger ();

  /** Constructor with name of model.
      It may or may not work if problem presolved
  */
  OsiRowCutDebugger (const OsiSolverInterface & si, const char * model);

  // Constructor with full solution (only integers need be correct)
  OsiRowCutDebugger (const OsiSolverInterface & si, const double * solution);
 
  /// Copy constructor 
  OsiRowCutDebugger (
    const OsiRowCutDebugger &);

  /// Assignment operator 
  OsiRowCutDebugger &
    operator=(
    const OsiRowCutDebugger& rhs);
  
  /// Destructor 
  virtual
    ~OsiRowCutDebugger ();
  //@}
      
private:
  
 // Private member methods


  /**@name Private methods */
  //@{

  //@}

  // Private member data

  /**@name Private member data */
  //@{
  /// Value of optimal solution
  double optimalValue_;
  /// number of columns in problem
  int numberColumns_;
  /// Whether integer or not
  bool * integerVariable_;
  /// Optimal column solution
  double * optimalSolution_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the OsiRowCut class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
void
OsiRowCutDebuggerUnitTest(const OsiSolverInterface * siP,    
			  const std::string & mpsDir);
  
#endif
