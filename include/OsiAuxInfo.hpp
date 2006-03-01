// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiAuxInfo_H
#define OsiAuxInfo_H

class OsiSolverInterface;

//#############################################################################
/** This class allows for a more structured use of algorithmic tweaking to
    an OsiSolverInterface.  It is designed to replace the simple use of
    appData_ pointer.

    This has been done to make it easier to use NonLinear solvers and other
    exotic beasts in a branch and bound mode.  After this class definition
    there is one for a derived class for just such a purpose.

*/

class OsiAuxInfo {
public:
  // Default Constructor 
  OsiAuxInfo (void * appData = NULL);

  // Copy Constructor 
  OsiAuxInfo (const OsiAuxInfo & rhs);
  // Destructor
  virtual ~OsiAuxInfo();
  
  /// Clone
  virtual OsiAuxInfo * clone() const;
  /// Assignment operator 
  OsiAuxInfo & operator=(const OsiAuxInfo& rhs);
  
  /// Get application data
  inline void * getApplicationData() const
  { return appData_;};
protected:
    /// Pointer to user-defined data structure
    void * appData_;
};
//#############################################################################
/** This class allows for the use of more exotic solvers e.g. Non-Linear or Volume.

    You can derive from this although at present I can't see the need.
*/

class OsiBabSolver : public OsiAuxInfo {
public:
  // Default Constructor 
  OsiBabSolver (int solverType=0);

  // Copy Constructor 
  OsiBabSolver (const OsiBabSolver & rhs);
  // Destructor
  virtual ~OsiBabSolver();
  
  /// Clone
  virtual OsiAuxInfo * clone() const;
  /// Assignment operator 
  OsiBabSolver & operator=(const OsiBabSolver& rhs);
  
  /// Update solver 
  inline void setSolver(OsiSolverInterface * solver)
  { solver_ = solver;};

  /** returns 0 if no heuristic solution, 1 if valid solution
      with better objective value than one passed in
      Sets solution values if good, sets objective value 
      Frees internal solution
  */
  int solution(double & objectiveValue,
		       double * newSolution);
  /// Set solution (and replace solver pointer if not null)
  void setSolution(const OsiSolverInterface * solver=NULL);

  /** Sets solver type
      0 - normal LP solver
      1 - DW - may also return heuristic solutions
      2 - NLP solver or similar - can't compute objective value just from solution
          check solver to see if feasible and what objective value is
          - may also return heuristic solution
      3 - NLP solver or similar - can't compute objective value just from solution
          check this (rather than solver) to see if feasible and what objective value is.
          Using Outer Approximation so called lp based
          - may also return heuristic solution
  */
  inline void setSolverType(int value)
  { solverType_=value;};
  /** gets solver type
      0 - normal LP solver
      1 - DW - may also return heuristic solutions
      2 - NLP solver or similar - can't compute objective value just from solution
          check this (rather than solver) to see if feasible and what objective value is
          - may also return heuristic solution
      3 - NLP solver or similar - can't compute objective value just from solution
          check this (rather than solver) to see if feasible and what objective value is.
          Using Outer Approximation so called lp based
          - may also return heuristic solution
  */
  inline int solverType() const
  { return solverType_;};
  /** Return true if getting solution may add cuts so hot start etc will
      be obsolete */
  inline bool solutionAddsCuts() const
  { return solverType_==3;};
  /** Returns true if can use solver objective or feasible values,
      otherwise use mipBound etc */
  inline bool solverAccurate() const
  { return solverType_==0||solverType_==2;};
  /// Returns true if can use reduced costs for fixing
  inline bool reducedCostsAccurate() const
  { return solverType_==0;};
  /// Get objective  (well mip bound)
  double mipBound() const;
  /// Returns true if node feasible
  bool mipFeasible() const;
  /// Set mip bound (only used for some solvers)
  inline void setMipBound(double value)
  { mipBound_ = value;};
  /// Says whether to get solution from this
  inline bool getSolutionFromInfo() const
  { return solverType_==3 && bestSolution_;};
  /// Says whether we want to try cuts at all
  inline bool tryCuts() const
  { return solverType_!=2;};
  /// Says whether we have a warm start (so can do strong branching)
  inline bool warmStart() const
  { return solverType_!=2;};
protected:
  /// Solver to use for getting/setting solutions etc
  const OsiSolverInterface * solver_;
  /** Solver type
      0 - normal LP solver
      1 - DW - may also return heuristic solutions
      2 - NLP solver or similar - can't compute objective value just from solution
          check this (rather than solver) to see if feasible and what objective value is
          - may also return heuristic solution
      3 - NLP solver or similar - can't compute objective value just from solution
          check this (rather than solver) to see if feasible and what objective value is.
          Using Outer Approximation so called lp based
          - may also return heuristic solution
  */
  int solverType_;
  /// Objective value of best solution (if there is one) (minimization)
  double bestObjectiveValue_;
  /// Best integer feasible solution
  double * bestSolution_;
  /// Current lower bound on solution ( if > 1.0e50 infeasible)
  double mipBound_;
};

#endif
