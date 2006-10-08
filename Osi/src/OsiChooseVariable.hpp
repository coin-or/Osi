// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiChooseVariable_H
#define OsiChooseVariable_H

#include <string>
#include <vector>

#include "CoinWarmStartBasis.hpp"
#include "OsiBranchingObject.hpp"

class OsiSolverInterface;
class OsiHotInfo;

/** This class chooses a variable to branch on

    The base class just chooses the variable and direction without strong branching but it 
    has information which would normally be used by strong branching e.g. to re-enter
    having fixed a variable but using same candidates for strong branching.

    The flow is :
    a) initialize the process.  This decides on strong branching list
       and stores indices of all infeasible objects  
    b) do strong branching on list.  If list is empty then just
       choose one candidate and return without strong branching.  If not empty then
       go through list and return best.  However we may find that the node is infeasible
       or that we can fix a variable.  If so we return and it is up to user to call
       again (after fixing a variable).
*/

class OsiChooseVariable  {
 
public:
    
  /// Default Constructor 
  OsiChooseVariable ();

  /// Constructor from solver (so we can set up arrays etc)
  OsiChooseVariable (const OsiSolverInterface * solver);

  /// Copy constructor 
  OsiChooseVariable (const OsiChooseVariable &);
   
  /// Assignment operator 
  OsiChooseVariable & operator= (const OsiChooseVariable& rhs);

  /// Clone
  virtual OsiChooseVariable * clone() const;

  /// Destructor 
  virtual ~OsiChooseVariable ();

  /** Sets up strong list and clears all if initialize is true.
      Returns number of infeasibilities. */
  int setupList ( OsiBranchingInformation *info, bool initialize);
  /** Choose a variable
      Returns - 
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from whichObject() and whichWay()
     We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
     Best object can be found from bestObjectIndex_
  */
  int chooseVariable( OsiBranchingInformation *info);
  /**  This is a utility function which does strong branching on
       a list of objects and stores the results in OsiHotInfo.objects.
       On entry the object sequence is stored in the OsiHotInfo object
       and maybe more.
       It returns the number of objects inspected.
  */
  int doStrongBranching( OsiBranchingInformation *info, int numberToDo,
			 OsiHotInfo * results);
  /// Given a candidate (bestObjectIndex_) fill in useful information e.g. estimates
  void updateInformation( OsiBranchingInformation *info);
  /// Objective value for feasible solution
  inline double goodObjectiveValue() const
  { return goodObjectiveValue_;};
  /// Estimate of up change or change on chosen if n-way
  inline double upChange() const
  { return upChange_;};
  /// Estimate of down change or max change on other possibilities if n-way
  inline double downChange() const
  { return downChange_;};
  /// Good solution - deleted by finalize
  inline const double * goodSolution() const
  { return goodSolution_;};
  /// Index of chosen object
  inline int bestObject() const
  { return bestObjectIndex_;};
  /// Set index of chosen object
  inline void setBestObject(int value)
  { bestObjectIndex_ = value;};
  /// Get the number of objects unsatisfied at this node - accurate on first pass
  inline int numberUnsatisfied() const
  {return numberUnsatisfied_;};
  /// Number of objects to choose for strong branching
  inline int numberStrong() const
  { return numberStrong_;};
  /// Set number of objects to choose for strong branching
  inline void setNumberStrong(int value)
  { numberStrong_ = value;};
  /// Number left on strong list
  inline int numberOnList() const
  { return numberOnList_;};
  /// Number of strong branches actually done 
  inline int numberStrongDone() const
  { return numberStrongDone_;};
  /// List of candidates
  inline const int * candidates() const
  { return list_;};
  /// Trust results from strong branching for changing bounds
  inline bool trustStrongForBound() const
  { return trustStrongForBound_;};
  /// Set trust results from strong branching for changing bounds
  inline void setTrustStrongForBound(bool yesNo)
  { trustStrongForBound_ = yesNo;};
  /// Trust results from strong branching for valid solution
  inline bool trustStrongForSolution() const
  { return trustStrongForSolution_;};
  /// Set trust results from strong branching for valid solution
  inline void setTrustStrongForSolution(bool yesNo)
  { trustStrongForSolution_ = yesNo;};
  /// Set solver and redo arrays
  void setSolver (const OsiSolverInterface * solver);
  /** Return status - 
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
     We can pick up branch from whichObject() and whichWay()
     We can pick up a forced branch (can change bound) from whichForcedObject() and whichForcedWay()
     If we have a solution then we can pick up from goodObjectiveValue() and goodSolution()
  */
  inline int status() const
  { return status_;};
  inline void setStatus(int value)
  { status_ = value;};


protected:
  // Data
  /// Objective value for feasible solution
  double goodObjectiveValue_;
  /// Estimate of up change or change on chosen if n-way
  double upChange_;
  /// Estimate of down change or max change on other possibilities if n-way
  double downChange_;
  /// Good solution - deleted by finalize
  double * goodSolution_;
  /// List of candidates
  int * list_;
  /// Useful array (for sorting etc)
  double * useful_;
  /// Pointer to solver
  const OsiSolverInterface * solver_;
  /* Status -
     -1 Node is infeasible
     0  Normal termination - we have a candidate
     1  All looks satisfied - no candidate
     2  We can change the bound on a variable - but we also have a strong branching candidate
     3  We can change the bound on a variable - but we have a non-strong branching candidate
     4  We can change the bound on a variable - no other candidates
  */
  int status_;
  /// Index of chosen object
  int bestObjectIndex_;
  /// The number of objects unsatisfied at this node.
  int numberUnsatisfied_;
  /// Number of objects to choose for strong branching
  int numberStrong_;
  /// Number left on strong list
  int numberOnList_;
  /// Number of strong branches actually done 
  int numberStrongDone_;
  /// List of unsatisfied objects - first numberOnList_ for strong branching
  /// Trust results from strong branching for changing bounds
  bool trustStrongForBound_;
  /// Trust results from strong branching for valid solution
  bool trustStrongForSolution_;
};

/** This class contains the result of strong branching on a variable
    When created it stores enough information for strong branching
*/

class OsiHotInfo  {
 
public:
    
  /// Default Constructor 
  OsiHotInfo ();

  /// Constructor from useful information
  OsiHotInfo ( const OsiBranchingInformation *info, int whichObject);

  /// Copy constructor 
  OsiHotInfo (const OsiHotInfo &);
   
  /// Assignment operator 
  OsiHotInfo & operator= (const OsiHotInfo& rhs);

  /// Clone
  virtual OsiHotInfo * clone() const;

  /// Destructor 
  virtual ~OsiHotInfo ();

  /** Fill in useful information after strong branch 
      for 2 way 0 is down, 1 is up
  */
  void updateInformation( const OsiSolverInterface * solver, int way);
  /// Original objective value
  inline double originalObjectiveValue() const
  { return originalObjectiveValue_;};
  /// Up change  - invalid if n-way
  inline double upChange() const
  { assert (!changes_); return change_[1];};
  /// Down change  - invalid if n-way
  inline double downChange() const
  { assert (!changes_); return change_[0];};
  /// Change on way k
  inline double change(int k) const
  { if (!changes_) return change_[k]; else return changes_[k];};

  /// Up iteration count  - invalid if n-way
  inline int upIterationCount() const
  { assert (!iterationCounts_); return iterationCount_[1];};
  /// Down iteration count  - invalid if n-way
  inline int downIterationCount() const
  { assert (!iterationCounts_); return iterationCount_[0];};
  /// Iteration count on way k
  inline int iterationCount(int k) const
  { if (!iterationCounts_) return iterationCount_[k]; else return iterationCounts_[k];};

  /// Up status  - invalid if n-way
  inline int upStatus() const
  { assert (!statuses_); return status_[1];};
  /// Down status  - invalid if n-way
  inline int downStatus() const
  { assert (!statuses_); return status_[0];};
  /// Status on way k
  inline int status(int k) const
  { if (!statuses_) return status_[k]; else return statuses_[k];};

protected:
  // Data
  /// Original objective value
  double originalObjectiveValue_;
  /// Objective changes for 2-way
  double change_[2];
  /// Iteration counts for 2-way
  int iterationCount_[2];
  /** Status for 2-way
      -1 - not done
      0 - feasible and finished
      1 -  infeasible
      2 - not finished
  */
  int status_[2];
  /// Objective changes for n-way
  double * changes_;
  /// Iteration counts for n-way
  int * iterationCounts_;
  /** Status for n-way
      -1 - not done
      0 - feasible and finished
      1 -  infeasible
      2 - not finished
  */
  int * statuses_;
  
};


#endif
