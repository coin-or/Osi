// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiClpSolverInterface.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "CoinHelperFunctions.hpp"
#include "ClpSimplex.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiOsiMessage.hpp"


#include  <time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
static double totalTime=0.0;
static double cpuTime()
{
  double cpu_temp;
#if defined(_MSC_VER)
  unsigned int ticksnow;        /* clock_t is same as int */
  
  ticksnow = (unsigned int)clock();
  
  cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else
  struct rusage usage;
  getrusage(RUSAGE_SELF,&usage);
  cpu_temp = usage.ru_utime.tv_sec;
  cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif
  return cpu_temp;
}

//#############################################################################
// Solve methods
//#############################################################################
void OsiClpSolverInterface::initialSolve()
{
  ClpSimplex solver;
  double time1 = cpuTime();
  solver.borrowModel(*modelPtr_);
  solver.scaling();
#if 0
  ClpDualRowSteepest steep;
  solver.setDualRowPivotAlgorithm(steep);
  solver.dual();
  basis_ = solver.getBasis();
  lastAlgorithm_=2; // dual
#else
  ClpPrimalColumnSteepest steep;
  solver.setPrimalColumnPivotAlgorithm(steep);
  solver.primal();
  basis_ = solver.getBasis();
  lastAlgorithm_=1; // primal
#endif
  solver.returnModel(*modelPtr_);
  time1 = cpuTime()-time1;
  totalTime += time1;
  std::cout<<time1<<" seconds - total "<<totalTime<<std::endl;
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::resolve()
{

  ClpSimplex solver;
  solver.borrowModel(*modelPtr_);
  solver.setBasis(basis_);
  solver.scaling();
  // For moment set log level to 0 for less output
  // Really there should be Osi call to do this
  solver.setLogLevel(0);

  ClpDualRowSteepest steep;
  solver.setDualRowPivotAlgorithm(steep);
  solver.dual();
  basis_ = solver.getBasis();
  lastAlgorithm_=2; // dual
  solver.returnModel(*modelPtr_);
}

//#############################################################################
// Parameter related methods
//#############################################################################

bool
OsiClpSolverInterface::setIntParam(OsiIntParam key, int value)
{
  return modelPtr_->setIntParam(key,value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setDblParam(OsiDblParam key, double value)
{
  return modelPtr_->setDblParam(key,value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::setStrParam(OsiStrParam key, const std::string & value)
{
  return modelPtr_->setStrParam(key,value);
}


//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getIntParam(OsiIntParam key, int& value) const 
{
  return modelPtr_->getIntParam(key,value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getDblParam(OsiDblParam key, double& value) const
{
  return modelPtr_->getDblParam(key,value);
}

//-----------------------------------------------------------------------------

bool
OsiClpSolverInterface::getStrParam(OsiStrParam key, std::string & value) const
{
  return modelPtr_->getStrParam(key,value);
}


//#############################################################################
// Methods returning info on how the solution process terminated
//#############################################################################

bool OsiClpSolverInterface::isAbandoned() const
{
  // *THINK*:
  // in *our* current setup if there are to many numerical difficulties, or
  // more precisely there is a severe error condition then CLP just outputs a
  // message and exits. We could change this, but will do it later.
  return false;
}

bool OsiClpSolverInterface::isProvenOptimal() const
{

  const int stat = modelPtr_->status();
  return (stat == 0);
}

bool OsiClpSolverInterface::isProvenPrimalInfeasible() const
{

  const int stat = modelPtr_->status();
  if (stat != 1)
     return false;
  return true;
}

bool OsiClpSolverInterface::isProvenDualInfeasible() const
{
  const int stat = modelPtr_->status();
  return stat == 2;
}

bool OsiClpSolverInterface::isPrimalObjectiveLimitReached() const
{
  double limit = 0.0;
  getDblParam(OsiPrimalObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = modelPtr_->objectiveValue();
  const int maxmin = modelPtr_->optimizationDirection();

  switch (lastAlgorithm_) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
   case 2: // dual simplex
     if (modelPtr_->status() == 0) // optimal
	return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
     return false;
   case 1: // primal simplex
     return maxmin > 0 ? (obj < limit) /*minim*/ : (obj > limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiClpSolverInterface::isDualObjectiveLimitReached() const
{

  double limit = 0.0;
  getDblParam(OsiDualObjectiveLimit, limit);
  if (limit > 1e30) {
    // was not ever set
    return false;
  }
   
  const double obj = modelPtr_->objectiveValue();
  const int maxmin = modelPtr_->optimizationDirection();

  switch (lastAlgorithm_) {
   case 0: // no simplex was needed
     return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
   case 1: // primal simplex
     if (modelPtr_->status() == 0) // optimal
	return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
     return false;
   case 2: // dual simplex
     if (modelPtr_->status() != 0 )
	// over dual limit
	return true;
     return maxmin > 0 ? (obj > limit) /*minim*/ : (obj < limit) /*maxim*/;
  }
  return false; // fake return
}

bool OsiClpSolverInterface::isIterationLimitReached() const
{
  const int stat = modelPtr_->status();
  return (stat == 3);
}

//#############################################################################
// WarmStart related methods
//#############################################################################

OsiWarmStart* OsiClpSolverInterface::getWarmStart() const
{

  return new OsiWarmStartBasis(basis_);
}

//-----------------------------------------------------------------------------

bool OsiClpSolverInterface::setWarmStart(const OsiWarmStart* warmstart)
{

  const OsiWarmStartBasis* ws =
    dynamic_cast<const OsiWarmStartBasis*>(warmstart);

  if (! ws)
    return false;
  basis_ = OsiWarmStartBasis(*ws);
  return true;

}

//#############################################################################
// Hotstart related methods (primarily used in strong branching)
//#############################################################################

void OsiClpSolverInterface::markHotStart()
{
  // *TEST*
  delete ws_;
  ws_ = dynamic_cast<OsiWarmStartBasis*>(getWarmStart());
  modelPtr_->getIntParam(OsiMaxNumIteration,itlimOrig_);
  int itlim;
  OsiSolverInterface::getIntParam(OsiMaxNumIterationHotStart, itlim);
  modelPtr_->setIntParam(OsiMaxNumIteration,itlim);

}

void OsiClpSolverInterface::solveFromHotStart()
{
  setWarmStart(ws_);
  resolve();
}

void OsiClpSolverInterface::unmarkHotStart()
{

  modelPtr_->setIntParam(OsiMaxNumIteration,itlimOrig_);
  delete ws_;
  ws_ = NULL;
}


// Trivial class for Branch and Bound

class OsiNode  {
  
public:
    
  // Default Constructor 
  OsiNode ();

  // Constructor from current state (and list of integers)
  // Also chooses branching variable (if none set to -1)
  OsiNode (OsiSolverInterface &model,
	   int numberIntegers, int * integer);
  
  // Copy constructor 
  OsiNode ( const OsiNode &);
   
  // Assignment operator 
  OsiNode & operator=( const OsiNode& rhs);

  // Destructor 
  ~OsiNode ();
  
  // Public data
  // Basis (should use tree, but not as wasteful as bounds!)
  OsiWarmStartBasis basis_;
  // Objective value
  double objectiveValue_;
  // Branching variable (0 is first integer)
  int variable_;
  // Way to branch - -1 down (first), 1 up, -2 down (second), 2 up (second)
  int way_;
  // Number of integers (for length of arrays)
  int numberIntegers_;
  // Current value
  double value_;
  // Now I must use tree
  // Bounds stored in full (for integers)
  int * lower_;
  int * upper_;
};


OsiNode::OsiNode() :
  basis_(OsiWarmStartBasis()),
  objectiveValue_(1.0e100),
  variable_(-100),
  way_(-1),
  numberIntegers_(0),
  value_(0.5),
  lower_(NULL),
  upper_(NULL)
{
}
OsiNode::OsiNode(OsiSolverInterface & model,
		 int numberIntegers, int * integer)
{
  const OsiWarmStartBasis* ws =
    dynamic_cast<const OsiWarmStartBasis*>(model.getWarmStart());

  assert (ws!=NULL); // make sure not volume
  basis_ = OsiWarmStartBasis(*ws);
  variable_=-1;
  way_=-1;
  numberIntegers_=numberIntegers;
  value_=0.0;
  if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
    objectiveValue_ = model.getObjSense()*model.getObjValue();
  } else {
    objectiveValue_ = 1.0e100;
    lower_ = NULL;
    upper_ = NULL;
    return; // node cutoff
  }
  lower_ = new int [numberIntegers_];
  upper_ = new int [numberIntegers_];
  assert (upper_!=NULL);
  const double * lower = model.getColLower();
  const double * upper = model.getColUpper();
  const double * solution = model.getColSolution();
  int i;
  // Hard coded integer tolerance
#define INTEGER_TOLERANCE 1.0e-6
  // Number of strong branching candidates
#define STRONG_BRANCHING 5
#ifdef STRONG_BRANCHING
  double upMovement[STRONG_BRANCHING];
  double downMovement[STRONG_BRANCHING];
  double solutionValue[STRONG_BRANCHING];
  int chosen[STRONG_BRANCHING];
  int iSmallest=0;
  // initialize distance from integer
  for (i=0;i<STRONG_BRANCHING;i++) {
    upMovement[i]=0.0;
    chosen[i]=-1;
  }
#endif
  variable_=-1;
  // This has hard coded integer tolerance
  double mostAway=INTEGER_TOLERANCE;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    lower_[i]=(int)lower[iColumn];
    upper_[i]=(int)upper[iColumn];
    double value = solution[iColumn];
    value = max(value,(double) lower_[i]);
    value = min(value,(double) upper_[i]);
    double nearest = floor(value+0.5);
    if (fabs(value-nearest)>mostAway) {
#ifdef STRONG_BRANCHING
      double away = fabs(value-nearest);
      if (away>upMovement[iSmallest]) {
	//add to list
	upMovement[iSmallest]=away;
	solutionValue[iSmallest]=value;
	chosen[iSmallest]=i;
	int j;
	iSmallest=-1;
	double smallest = 1.0;
	for (j=0;j<STRONG_BRANCHING;j++) {
	  if (upMovement[j]<smallest) {
	    smallest=upMovement[j];
	    iSmallest=j;
	  }
	}
      }
#else
      mostAway=fabs(value-nearest);
      variable_=i;
      value_=value;
      if (value<=nearest)
	way_=1; // up
      else
	way_=-1; // down
#endif
    }
  }
#ifdef STRONG_BRANCHING
  int numberStrong=0;
  for (i=0;i<STRONG_BRANCHING;i++) {
    if (chosen[i]>=0) { 
      numberStrong ++;
      variable_ = chosen[i];
    }
  }
  if (numberStrong==1) {
    // just one - makes it easy
    int iColumn = integer[variable_];
    double value = solution[iColumn];
    value = max(value,(double) lower_[variable_]);
    value = min(value,(double) upper_[variable_]);
    double nearest = floor(value+0.5);
    value_=value;
    if (value<=nearest)
      way_=1; // up
    else
      way_=-1; // down
  } else if (numberStrong) {
    // more than one - choose
    bool chooseOne=true;
    model.markHotStart();
    for (i=0;i<STRONG_BRANCHING;i++) {
      int iInt = chosen[i];
      if (iInt>=0) {
	int iColumn = integer[iInt];
	double value = solutionValue[i]; // value of variable in original
	double objectiveChange;
	value = max(value,(double) lower_[iInt]);
	value = min(value,(double) upper_[iInt]);

	// try down

	model.setColUpper(iColumn,floor(value));
	model.solveFromHotStart();
	model.setColUpper(iColumn,upper_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	downMovement[i]=objectiveChange;

	// try up

	model.setColLower(iColumn,ceil(value));
	model.solveFromHotStart();
	model.setColLower(iColumn,lower_[iInt]);
	if (model.isProvenOptimal()&&!model.isDualObjectiveLimitReached()) {
	  objectiveChange = model.getObjSense()*model.getObjValue()
	    - objectiveValue_;
	} else {
	  objectiveChange = 1.0e100;
	}
	upMovement[i]=objectiveChange;
	
	/* Possibilities are:
	   Both sides feasible - store
	   Neither side feasible - set objective high and exit
	   One side feasible - change bounds and resolve
	*/
	bool solveAgain=false;
	if (upMovement[i]<1.0e100) {
	  if(downMovement[i]<1.0e100) {
	    // feasible - no action
	  } else {
	    // up feasible, down infeasible
	    solveAgain = true;
	    model.setColLower(iColumn,ceil(value));
	  }
	} else {
	  if(downMovement[i]<1.0e100) {
	    // down feasible, up infeasible
	    solveAgain = true;
	    model.setColUpper(iColumn,floor(value));
	  } else {
	    // neither side feasible
	    objectiveValue_=1.0e100;
	    chooseOne=false;
	    break;
	  }
	}
	if (solveAgain) {
	  // need to solve problem again - signal this
	  variable_ = numberIntegers;
	  chooseOne=false;
	  break;
	}
      }
    }
    if (chooseOne) {
      // choose the one that makes most difference both ways
      double best = -1.0;
      double best2 = -1.0;
      for (i=0;i<STRONG_BRANCHING;i++) {
	int iInt = chosen[i];
	if (iInt>=0) {
	  model.messageHandler()->message(OSI_BAB_STRONG,model.messages())
	    <<i<<iInt<<downMovement[i]<<upMovement[i]
	    <<solutionValue[i]
	    <<OsiMessageEol;
	  bool better = false;
	  if (min(upMovement[i],downMovement[i])>best) {
	    // smaller is better
	    better=true;
	  } else if (min(upMovement[i],downMovement[i])>best-1.0e-5) {
	    if (max(upMovement[i],downMovement[i])>best2+1.0e-5) {
	      // smaller is about same, but larger is better
	      better=true;
	    }
	  }
	  if (better) {
	    best = min(upMovement[i],downMovement[i]);
	    best2 = max(upMovement[i],downMovement[i]);
	    variable_ = iInt;
	    double value = solutionValue[i];
	    value = max(value,(double) lower_[variable_]);
	    value = min(value,(double) upper_[variable_]);
	    value_=value;
	    if (upMovement[i]<=downMovement[i])
	      way_=1; // up
	    else
	      way_=-1; // down
	  }
	}
      }
    }
    // Delete the snapshot
    model.unmarkHotStart();
  }
#endif
}

OsiNode::OsiNode(const OsiNode & rhs) 
{  
  basis_=rhs.basis_;
  objectiveValue_=rhs.objectiveValue_;
  variable_=rhs.variable_;
  way_=rhs.way_;
  numberIntegers_=rhs.numberIntegers_;
  value_=rhs.value_;
  lower_=NULL;
  upper_=NULL;
  if (rhs.lower_!=NULL) {
    lower_ = new int [numberIntegers_];
    upper_ = new int [numberIntegers_];
    assert (upper_!=NULL);
    memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
    memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
  }
}

OsiNode &
OsiNode::operator=(const OsiNode & rhs)
{
  if (this != &rhs) {
    basis_=rhs.basis_;
    objectiveValue_=rhs.objectiveValue_;
    variable_=rhs.variable_;
    way_=rhs.way_;
    numberIntegers_=rhs.numberIntegers_;
    value_=rhs.value_;
    delete [] lower_;
    delete [] upper_;
    lower_=NULL;
    upper_=NULL;
    if (rhs.lower_!=NULL) {
      lower_ = new int [numberIntegers_];
      upper_ = new int [numberIntegers_];
      assert (upper_!=NULL);
      memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
      memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
    }
  }
  return *this;
}


OsiNode::~OsiNode ()
{
  delete [] lower_;
  delete [] upper_;
}

#include <vector>

// Vector of OsiNodes 
typedef std::vector<OsiNode>    OsiVectorNode;

// Invoke solver's built-in enumeration algorithm
void 
OsiClpSolverInterface::branchAndBound() {
  // solve LP
  initialSolve();
  if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
    // This is a really simple Branch and Bound code - mainly
    // to test strong branching
    // I should look at STL more to allow other than depth first
    int numberIntegers=0;
    int numberColumns = getNumCols();
    int iColumn;
    int i;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn))
	numberIntegers++;
    }
    if (!numberIntegers) {
      handler_->message(OSI_BAB_NOINT,messages_)
	<<OsiMessageEol;
      return;
    }
    int * which = new int[numberIntegers]; // which variables are integer
    numberIntegers=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn))
	which[numberIntegers++]=iColumn;
    }

    // empty tree
    OsiVectorNode branchingTree;

    // Add continuous to it;
    branchingTree.push_back(OsiNode(*this,numberIntegers,which));

    // For printing totals
    int numberIterations=0;
    int numberNodes =0;

    OsiNode bestNode;
    // while until nothing on stack
    while (branchingTree.size()) {
      // last node
      OsiNode node = branchingTree.back();
      branchingTree.pop_back();
      numberNodes++;
      if (node.variable_>=0) {
	// branch - do bounds
	for (i=0;i<numberIntegers;i++) {
	  iColumn=which[i];
	  setColBounds( iColumn,node.lower_[i],node.upper_[i]);
	}
	// move basis
	basis_=node.basis_;
	// do branching variable
	if (node.way_<0) {
	  setColUpper(which[node.variable_],floor(node.value_));
	  // now push back node if more to come
	  if (node.way_==-1) { 
	    node.way_=+2;	  // Swap direction
	    branchingTree.push_back(node);
	  }
	} else {
	  setColLower(which[node.variable_],ceil(node.value_));
	  // now push back node if more to come
	  if (node.way_==1) { 
	    node.way_=-2;	  // Swap direction
	    branchingTree.push_back(node);
	  }
	}
	// solve
	resolve();
	numberIterations += getIterationCount();
	if (!isIterationLimitReached()) {
	  OsiNode newNode(*this,numberIntegers,which);
	  // something extra may have been fixed by strong branching
	  // if so go round again
	  while (newNode.variable_==numberIntegers) {
	    resolve();
	    newNode = OsiNode(*this,numberIntegers,which);
	  }
	  if (newNode.objectiveValue_<1.0e100) {
	    // push on stack
	    branchingTree.push_back(newNode);
	  }
	} else {
	  // maximum iterations - exit
      handler_->message(OSI_BAB_MAXITS,messages_)
	<<OsiMessageEol;
	  break;
	}
      } else {
	// integer solution - save
	bestNode = node;
	// set cutoff (hard coded tolerance)
	setDblParam(OsiDualObjectiveLimit,bestNode.objectiveValue_-1.0e-5);
	handler_->message(OSI_BAB_SOLUTION,messages_)
	  <<bestNode.objectiveValue_<<numberIterations
	  <<numberNodes
	  <<OsiMessageEol;
      }
    }
    handler_->message(OSI_BAB_END,messages_)
	  <<numberIterations
	  <<numberNodes
	  <<OsiMessageEol;
    if (bestNode.numberIntegers_) {
      // we have a solution restore
      // do bounds
      for (i=0;i<numberIntegers;i++) {
	iColumn=which[i];
	setColBounds( iColumn,bestNode.lower_[i],bestNode.upper_[i]);
      }
      // move basis
      basis_=bestNode.basis_;
      resolve();
    }
    delete [] which;
  } else {
    handler_->message(OSI_BAB_INFEAS,messages_)
    <<OsiMessageEol;
    throw CoinError("The LP relaxation is infeasible or too expensive",
		    "branchAndBound", "OsiClpSolverInterface");
  }
}

//#############################################################################
// Problem information methods (original data)
//#############################################################################

//------------------------------------------------------------------
const char * OsiClpSolverInterface::getRowSense() const
{
  extractSenseRhsRange();
  return rowsense_;
}
//------------------------------------------------------------------
const double * OsiClpSolverInterface::getRightHandSide() const
{
  extractSenseRhsRange();
  return rhs_;
}
//------------------------------------------------------------------
const double * OsiClpSolverInterface::getRowRange() const
{
  extractSenseRhsRange();
  return rowrange_;
}
//------------------------------------------------------------------
// Return information on integrality
//------------------------------------------------------------------
bool OsiClpSolverInterface::isContinuous(int colNumber) const
{
  if ( integerInformation_==NULL ) return true;
  if ( integerInformation_[colNumber]==0 ) return true;
  return false;
}
//------------------------------------------------------------------

//------------------------------------------------------------------
// Row and column copies of the matrix ...
//------------------------------------------------------------------
const OsiPackedMatrix * OsiClpSolverInterface::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL ) {
    matrixByRow_ = new OsiPackedMatrix(); 
    matrixByRow_->reverseOrderedCopyOf(*modelPtr_->matrix());
    matrixByRow_->removeGaps();
#if 0
    OsiPackedMatrix back;
    std::cout<<"start check"<<std::endl;
    back.reverseOrderedCopyOf(*matrixByRow_);
    modelPtr_->matrix()->isEquivalent2(back);
    std::cout<<"stop check"<<std::endl;
#endif
  }
  return matrixByRow_;
}

const OsiPackedMatrix * OsiClpSolverInterface::getMatrixByCol() const
{
  return modelPtr_->matrix();
}

//------------------------------------------------------------------
std::vector<double*> OsiClpSolverInterface::getDualRays(int maxNumRays) const
{
  return std::vector<double*>(1, modelPtr_->infeasibilityRay());
}
//------------------------------------------------------------------
std::vector<double*> OsiClpSolverInterface::getPrimalRays(int maxNumRays) const
{
  return std::vector<double*>(1, modelPtr_->unboundedRay());
}
//------------------------------------------------------------------

//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setColSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  double * lower = modelPtr_->columnLower();
  double * upper = modelPtr_->columnUpper();
  while (indexFirst != indexLast) {
    const int iCol=*indexFirst++;
    lower[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
    upper[iCol]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowType(int i, char sense, double rightHandSide,
				  double range)
{
  // *TEST*
  double lower, upper;
  convertSenseToBound(sense, rightHandSide, range, lower, upper);
  setRowBounds(i, lower, upper);
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowSetBounds(const int* indexFirst,
					    const int* indexLast,
					    const double* boundList)
{
  double * lower = modelPtr_->rowLower();
  double * upper = modelPtr_->rowUpper();
  const int len = indexLast - indexFirst;
  while (indexFirst != indexLast) {
    const int iRow=*indexFirst++;
    lower[iRow]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(*boundList++, -OsiClpInfinity, OsiClpInfinity);
  }
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    indexFirst -= len;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      convertBoundToSense(lower[iRow], upper[iRow],
			  rowsense_[iRow], rhs_[iRow], rowrange_[iRow]);
    }
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setRowSetTypes(const int* indexFirst,
				      const int* indexLast,
				      const char* senseList,
				      const double* rhsList,
				      const double* rangeList)
{
  double * lower = modelPtr_->rowLower();
  double * upper = modelPtr_->rowUpper();
  const int len = indexLast - indexFirst;
  while (indexFirst != indexLast) {
    const int iRow= *indexFirst++;
    convertSenseToBound(*senseList++, *rhsList++, *rangeList++,
			lower[iRow], upper[iRow]);
  }
  if (rowsense_ != NULL) {
    assert ((rhs_ != NULL) && (rowrange_ != NULL));
    indexFirst -= len;
    senseList -= len;
    rhsList -= len;
    rangeList -= len;
    while (indexFirst != indexLast) {
      const int iRow=*indexFirst++;
      rowsense_[iRow] = *senseList++;
      rhs_[iRow] = *rhsList++;
      rowrange_[iRow] = *rangeList++;
    }
  }
}
//#############################################################################
void
OsiClpSolverInterface::setContinuous(int index)
{

  if (integerInformation_) {
    integerInformation_[index]=0;
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setInteger(int index)
{
  if (!integerInformation_) {
    integerInformation_ = new char[modelPtr_->numberColumns()];
    CoinFillN ( integerInformation_, modelPtr_->numberColumns(),(char) 0);
  }
  integerInformation_[index]=1;
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setContinuous(const int* indices, int len)
{
  if (integerInformation_) {
    int i;
    for (i=0; i<len;i++) {
      integerInformation_[i]=0;
    }
  }
}
//-----------------------------------------------------------------------------
void
OsiClpSolverInterface::setInteger(const int* indices, int len)
{
  if (!integerInformation_) {
    integerInformation_ = new char[modelPtr_->numberColumns()];
    CoinFillN ( integerInformation_, modelPtr_->numberColumns(),(char) 0);
  }
  int i;
  for (i=0; i<len;i++) {
    integerInformation_[indices[i]]=1;
  }
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setColSolution(const double * cs) 
{
  CoinDisjointCopyN(cs,modelPtr_->numberColumns(),
		    modelPtr_->primalColumnSolution());
}
//-----------------------------------------------------------------------------
void OsiClpSolverInterface::setRowPrice(const double * rs) 
{
  CoinDisjointCopyN(rs,modelPtr_->numberRows(),
		    modelPtr_->dualRowSolution());
}

//#############################################################################
// Problem modifying methods (matrix)
//#############################################################################
void 
OsiClpSolverInterface::addCol(const OsiPackedVectorBase& vec,
			      const double collb, const double colub,   
			      const double obj)
{
  int numberColumns = modelPtr_->numberColumns();
  modelPtr_->resize(modelPtr_->numberRows(),numberColumns+1);
  basis_.resize(modelPtr_->numberRows(),numberColumns+1);
  setColBounds(numberColumns,collb,colub);
  setObjCoeff(numberColumns,obj);
  modelPtr_->matrix()->appendCol(vec);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addCols(const int numcols,
			       const OsiPackedVectorBase * const * cols,
			       const double* collb, const double* colub,   
			       const double* obj)
{
  int numberColumns = modelPtr_->numberColumns();
  modelPtr_->resize(modelPtr_->numberRows(),numberColumns+numcols);
  basis_.resize(modelPtr_->numberRows(),numberColumns+numcols);
  double * lower = modelPtr_->columnLower()+numberColumns;
  double * upper = modelPtr_->columnUpper()+numberColumns;
  double * objective = modelPtr_->objective()+numberColumns;
  int iCol;
  for (iCol = 0; iCol < numcols; iCol++) {
    lower[iCol]= forceIntoRange(collb[iCol], -OsiClpInfinity, OsiClpInfinity);
    upper[iCol]= forceIntoRange(colub[iCol], -OsiClpInfinity, OsiClpInfinity);
    objective[iCol] = obj[iCol];
  }
  modelPtr_->matrix()->appendCols(numcols,cols);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::deleteCols(const int num, const int * columnIndices)
{
  modelPtr_->deleteColumns(num,columnIndices);
  basis_.deleteColumns(num,columnIndices);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const double rowlb, const double rowub)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  setRowBounds(numberRows,rowlb,rowub);
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRow(const OsiPackedVectorBase& vec,
			      const char rowsen, const double rowrhs,   
			      const double rowrng)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+1,modelPtr_->numberColumns());
  basis_.resize(numberRows+1,modelPtr_->numberColumns());
  double rowlb, rowub;
  convertSenseToBound(rowsen, rowrhs, rowrng, rowlb, rowub);
  setRowBounds(numberRows,rowlb,rowub);
  modelPtr_->matrix()->appendRow(vec);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const double* rowlb, const double* rowub)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    lower[iRow]= forceIntoRange(rowlb[iRow], -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(rowub[iRow], -OsiClpInfinity, OsiClpInfinity);
  }
  modelPtr_->matrix()->appendRows(numrows,rows);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::addRows(const int numrows,
			       const OsiPackedVectorBase * const * rows,
			       const char* rowsen, const double* rowrhs,   
			       const double* rowrng)
{
  int numberRows = modelPtr_->numberRows();
  modelPtr_->resize(numberRows+numrows,modelPtr_->numberColumns());
  basis_.resize(numberRows+numrows,modelPtr_->numberColumns());
  double * lower = modelPtr_->rowLower()+numberRows;
  double * upper = modelPtr_->rowUpper()+numberRows;
  int iRow;
  for (iRow = 0; iRow < numrows; iRow++) {
    double rowlb, rowub;
    convertSenseToBound(rowsen[iRow], rowrhs[iRow], rowrng[iRow], 
			rowlb, rowub);
    lower[iRow]= forceIntoRange(rowlb, -OsiClpInfinity, OsiClpInfinity);
    upper[iRow]= forceIntoRange(rowub, -OsiClpInfinity, OsiClpInfinity);
  }
  modelPtr_->matrix()->appendRows(numrows,rows);
  freeCachedResults();
}
//-----------------------------------------------------------------------------
void 
OsiClpSolverInterface::deleteRows(const int num, const int * rowIndices)
{
  modelPtr_->deleteRows(num,rowIndices);
  basis_.deleteRows(num,rowIndices);
  freeCachedResults();
}

//#############################################################################
// Methods to input a problem
//#############################################################################

void
OsiClpSolverInterface::loadProblem(const OsiPackedMatrix& matrix,
				   const double* collb, const double* colub,   
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  modelPtr_->loadProblem(matrix, collb, colub, obj, rowlb, rowub);
  freeCachedResults();

}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::assignProblem(OsiPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     double*& rowlb, double*& rowub)
{
   modelPtr_->loadProblem(*matrix, collb, colub, obj, rowlb, rowub);
   freeCachedResults();
   delete matrix;   matrix = NULL;
   delete[] collb;  collb = NULL;
   delete[] colub;  colub = NULL;
   delete[] obj;    obj = NULL;
   delete[] rowlb;  rowlb = NULL;
   delete[] rowub;  rowub = NULL;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const OsiPackedMatrix& matrix,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
   assert( rowsen != NULL );
   assert( rowrhs != NULL );
   int numrows = matrix.getNumRows();
   double * rowlb = new double[numrows];
   double * rowub = new double[numrows];
   for (int i = numrows-1; i >= 0; --i) {   
      convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
   }
   modelPtr_->loadProblem(matrix, collb, colub, obj, rowlb, rowub);
   freeCachedResults();
   delete [] rowlb;
   delete [] rowub;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::assignProblem(OsiPackedMatrix*& matrix,
				     double*& collb, double*& colub,
				     double*& obj,
				     char*& rowsen, double*& rowrhs,
				     double*& rowrng)
{
   loadProblem(*matrix, collb, colub, obj, rowsen, rowrhs, rowrng);
   delete matrix;   matrix = NULL;
   delete[] collb;  collb = NULL;
   delete[] colub;  colub = NULL;
   delete[] obj;    obj = NULL;
   delete[] rowsen; rowsen = NULL;
   delete[] rowrhs; rowrhs = NULL;
   delete[] rowrng; rowrng = NULL;
}

//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const double* rowlb, const double* rowub)
{
  modelPtr_->loadProblem(numcols, numrows, start,  index,
	    value, collb, colub, obj,
	    rowlb,  rowub);
  freeCachedResults();
}
//-----------------------------------------------------------------------------

void
OsiClpSolverInterface::loadProblem(const int numcols, const int numrows,
				   const int* start, const int* index,
				   const double* value,
				   const double* collb, const double* colub,
				   const double* obj,
				   const char* rowsen, const double* rowrhs,   
				   const double* rowrng)
{
   assert( rowsen != NULL );
   assert( rowrhs != NULL );
   double * rowlb = new double[numrows];
   double * rowub = new double[numrows];
   for (int i = numrows-1; i >= 0; --i) {   
      convertSenseToBound(rowsen[i],rowrhs[i],rowrng[i],rowlb[i],rowub[i]);
   }
   modelPtr_->loadProblem(numcols, numrows, start,  index,
	     value, collb, colub, obj,
	     rowlb,  rowub);
   freeCachedResults();
   delete[] rowlb;
   delete[] rowub;
}

//-----------------------------------------------------------------------------
// Write mps files
//-----------------------------------------------------------------------------

void OsiClpSolverInterface::writeMps(const char * filename,
				     const char * extension) const
{
  std::string f(filename);
  std::string e(extension);
  std::string fullname;
  if (e!="") {
    fullname = f + "." + e;
  } else {
    // no extension so no trailing period
    fullname = f;
  }
  throw CoinError("Sorry, the Clp Algorithm doesn't implement writeMps",
		  "writeMps", "OsiClpSolverInterface");

}

//#############################################################################
// CLP specific public interfaces
//#############################################################################

ClpModel * OsiClpSolverInterface::getModelPtr() const
{
  freeCachedResults();
  return modelPtr_;
}

//------------------------------------------------------------------- 

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::OsiClpSolverInterface ()
:
OsiSolverInterface(),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
basis_(),  
itlimOrig_(9999999),
lastAlgorithm_(0),
matrixByRow_(NULL),
integerInformation_(NULL)
{
  modelPtr_ = new ClpModel();
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * OsiClpSolverInterface::clone(bool CopyData) const
{
   if (CopyData) {
      return new OsiClpSolverInterface(*this);
   } else {
      return new OsiClpSolverInterface();
   }
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::OsiClpSolverInterface (
                  const OsiClpSolverInterface & rhs)
:
OsiSolverInterface(rhs),
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
ws_(NULL),
basis_(),
itlimOrig_(9999999),
lastAlgorithm_(0),
matrixByRow_(NULL),
integerInformation_(NULL)
{
  if ( rhs.modelPtr_  ) 
    modelPtr_ = new ClpModel(*rhs.modelPtr_);
  else
    modelPtr_ = new ClpModel();
  if ( rhs.ws_ ) 
    ws_ = new OsiWarmStartBasis(*rhs.ws_);
  basis_ = rhs.basis_;
  if (rhs.integerInformation_) {
    int numberColumns = modelPtr_->numberColumns();
    integerInformation_ = new char[numberColumns];
    memcpy(integerInformation_,rhs.integerInformation_,
	   numberColumns*sizeof(char));
  }
}


//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiClpSolverInterface::~OsiClpSolverInterface ()
{
  freeCachedResults();
  delete modelPtr_;
  delete ws_;
  delete [] integerInformation_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiClpSolverInterface &
OsiClpSolverInterface::operator=(const OsiClpSolverInterface& rhs)
{
  if (this != &rhs) {    
    OsiSolverInterface::operator=(rhs);
    freeCachedResults();
    delete modelPtr_;
    delete ws_;
    if ( rhs.modelPtr_  ) 
      modelPtr_ = new ClpModel(*rhs.modelPtr_);
    
    if ( rhs.ws_ ) 
      ws_ = new OsiWarmStartBasis(*rhs.ws_);
    basis_ = rhs.basis_;
  }
  return *this;
}

//#############################################################################
// Applying cuts
//#############################################################################

void OsiClpSolverInterface::applyRowCut( const OsiRowCut & rowCut )
{
  const OsiPackedVector & row=rowCut.row();
  addRow(row ,  rowCut.lb(),rowCut.ub());
}

//-----------------------------------------------------------------------------

void OsiClpSolverInterface::applyColCut( const OsiColCut & cc )
{
  double * lower = modelPtr_->columnLower();
  double * upper = modelPtr_->columnUpper();
  const OsiPackedVector & lbs = cc.lbs();
  const OsiPackedVector & ubs = cc.ubs();
  int i;

  for ( i=0; i<lbs.getNumElements(); i++ ) {
    int iCol = lbs.getIndices()[i];
    double value = lbs.getElements()[i];
    if ( value > lower[iCol] )
      lower[iCol]= value;
  }
  for ( i=0; i<ubs.getNumElements(); i++ ) {
    int iCol = ubs.getIndices()[i];
    double value = ubs.getElements()[i];
    if ( value < upper[iCol] )
      upper[iCol]= value;
  }
}
//#############################################################################
// Private methods
//#############################################################################


//------------------------------------------------------------------- 

void OsiClpSolverInterface::freeCachedResults() const
{  
  delete [] rowsense_;
  delete [] rhs_;
  delete [] rowrange_;
  delete matrixByRow_;
  delete ws_;
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
  matrixByRow_=NULL;
  ws_ = NULL;
}

//------------------------------------------------------------------
void OsiClpSolverInterface::extractSenseRhsRange() const
{
  if (rowsense_ == NULL) {
    // all three must be NULL
    assert ((rhs_ == NULL) && (rowrange_ == NULL));


    int nr=modelPtr_->numberRows();
    rowsense_ = new char[nr];
    rhs_ = new double[nr];
    rowrange_ = new double[nr];
    std::fill(rowrange_,rowrange_+nr,0.0);

    const double * lb = modelPtr_->rowLower();
    const double * ub = modelPtr_->rowUpper();

    int i;
    for ( i=0; i<nr; i++ ) {
      convertBoundToSense(lb[i], ub[i], rowsense_[i], rhs_[i], rowrange_[i]);
    }
  }
}
// Set language
void 
OsiClpSolverInterface::newLanguage(OsiMessages::Language language)
{
  modelPtr_->newLanguage(language);
  OsiSolverInterface::newLanguage(language);
}
//#############################################################################
