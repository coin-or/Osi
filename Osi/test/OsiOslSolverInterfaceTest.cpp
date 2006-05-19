// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiConfig.h"

#include <cassert>

#include "OsiOslSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinMessage.hpp"

// below needed for pathetic branch and bound code
#include <vector>
#include <map>
#ifndef _MSC_VER
  using std::max;
  using std::min;
#endif
  
// Added so build windows build with dsp files works,
// when not building with cplex.
#ifdef COIN_HAS_OSL

// Trivial class for Branch and Bound

class OsiNodeSimple  {
  
public:
    
  // Default Constructor 
  OsiNodeSimple ();

  // Constructor from current state (and list of integers)
  // Also chooses branching variable (if none set to -1)
  OsiNodeSimple (OsiSolverInterface &model,
	   int numberIntegers, int * integer);
  
  // Copy constructor 
  OsiNodeSimple ( const OsiNodeSimple &);
   
  // Assignment operator 
  OsiNodeSimple & operator=( const OsiNodeSimple& rhs);

  // Destructor 
  ~OsiNodeSimple ();
  
  // Public data
  // Basis (should use tree, but not as wasteful as bounds!)
  CoinWarmStartBasis basis_;
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


OsiNodeSimple::OsiNodeSimple() :
  basis_(CoinWarmStartBasis()),
  objectiveValue_(1.0e100),
  variable_(-100),
  way_(-1),
  numberIntegers_(0),
  value_(0.5),
  lower_(NULL),
  upper_(NULL)
{
}
OsiNodeSimple::OsiNodeSimple(OsiSolverInterface & model,
		 int numberIntegers, int * integer)
{
  const CoinWarmStartBasis* ws =
    dynamic_cast<const CoinWarmStartBasis*>(model.getWarmStart());

  assert (ws!=NULL); // make sure not volume
  basis_ = CoinWarmStartBasis(*ws);
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
	  std::cout<<"Strong branching on "
		   <<i<<""<<iInt<<" down "<<downMovement[i]
		   <<" up "<<upMovement[i]
		   <<" value "<<solutionValue[i]
		   <<std::endl;
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

OsiNodeSimple::OsiNodeSimple(const OsiNodeSimple & rhs) 
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

OsiNodeSimple &
OsiNodeSimple::operator=(const OsiNodeSimple & rhs)
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


OsiNodeSimple::~OsiNodeSimple ()
{
  delete [] lower_;
  delete [] upper_;
}

#include <vector>

// Vector of OsiNodeSimples 
typedef std::vector<OsiNodeSimple>    OsiVectorNode;

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif
class OsiOslMessageTest :
   public CoinMessageHandler {

public:
  virtual int print() ;
};

int
OsiOslMessageTest::print()
{
  static int numberOptimal=0,numberInfeasible=0;
  if (currentSource()=="Osl") {
    if (currentMessage().externalNumber()==1) { 
      numberOptimal++;
      if (numberOptimal%100==0)
	return CoinMessageHandler::print(); // print
    } else if (currentMessage().externalNumber()==3000) { 
      numberInfeasible++;
    }
  }  else if (currentSource()=="Osi") {
    if (currentMessage().externalNumber()==5) 
      std::cout<<"End of search trapped, "
	       <<numberOptimal<<" optimal Lp's, "
	       <<numberInfeasible
	       <<" infeasible (including strong branching)"<<std::endl;
    return CoinMessageHandler::print();
  }
  return 0;
}
#include <stdio.h>
// OSL function to trap messages
extern "C" void trapMessages(EKKModel * model, int  msgno, int nreal,
                            double * rvec, int nint, int * ivec,
                            int nchar, char * cvec)
{
  // break out message (assumes prefix there)
  char severity = cvec[nchar*128+8];
  char mainPart[153];
  strncpy(mainPart,cvec+nchar*128+9,152);
  mainPart[152]='\0';
  int length = strlen(mainPart);
  int i;
  for (i=length-1;i>=0;i--) {
    if (mainPart[i]!=' '&&mainPart[i]!='\t')
      break;
  }
  if (i>=0) {
    mainPart[i+1]='\0';
    CoinMessageHandler * handler =  (CoinMessageHandler *) ekk_userData(model);
    handler->message(msgno,"Osl",mainPart,severity);
    for (i=0;i<nreal;i++) 
      *handler<<rvec[i];
    for (i=0;i<nint;i++) 
      *handler<<ivec[i];
    for (i=0;i<nchar;i++) 
      *handler<<cvec+128*i;
    handler->finish();
  }
  return;
}
//--------------------------------------------------------------------------
// test EKKsolution methods.
void
OsiOslSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

  // Test default constructor
  {
    assert( OsiOslSolverInterface::getNumInstances()==0 );
    OsiOslSolverInterface m;
    assert( m.modelPtr_==NULL );
    assert( m.rowsense_==NULL );
    assert( m.rhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.matrixByColumn_==NULL );
    assert( OsiOslSolverInterface::getNumInstances()==1 );
    assert( m.getApplicationData() == NULL );
    int i=2346;
    m.setApplicationData(&i);
    assert( *((int *)(m.getApplicationData())) == i );
  }
  assert( OsiOslSolverInterface::getNumInstances()==0 );

  
  {    
    CoinRelFltEq eq;
    OsiOslSolverInterface m;
    assert( OsiOslSolverInterface::getNumInstances()==1 );
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    int ad = 13579;
    m.setApplicationData(&ad);
    assert( *((int *)(m.getApplicationData())) == ad );
    
    {
      OsiOslSolverInterface im;    
      assert( im.modelPtr_==NULL );
      
      assert( im.getNumCols() == 0 ); 
      
      assert( im.getModelPtr()!=NULL );
      assert( im.getMutableModelPtr()!=NULL );
      assert( im.getModelPtr() == im.getMutableModelPtr() );
    }
    
    // Test copy constructor and assignment operator
    {
      OsiOslSolverInterface lhs;
      {      
        assert( *((int *)(m.getApplicationData())) == ad );
        OsiOslSolverInterface im(m);        
        assert( *((int *)(im.getApplicationData())) == ad );

        OsiOslSolverInterface imC1(im);
        assert( imC1.getMutableModelPtr()!=im.getMutableModelPtr() );
        assert( imC1.getModelPtr()!=im.getModelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        assert( *((int *)(imC1.getApplicationData())) == ad ); 
        
        //im.setModelPtr(m);
        
        
        OsiOslSolverInterface imC2(im);
        assert( imC2.getMutableModelPtr()!=im.getMutableModelPtr() );
        assert( imC2.getModelPtr()!=im.getModelPtr() );
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        assert( *((int *)(imC2.getApplicationData())) == ad ); 
        
        assert( imC2.getMutableModelPtr()!=imC1.getMutableModelPtr() );
        assert( imC2.getModelPtr()!=imC1.getModelPtr() );
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope
      
      assert( lhs.getMutableModelPtr() != m.getMutableModelPtr() );
      assert( lhs.getModelPtr() != m.getModelPtr() );
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
      assert( *((int *)(lhs.getApplicationData())) == ad );
    }
    
    // Test clone
    {
      OsiOslSolverInterface oslSi(m);
      OsiSolverInterface * siPtr = &oslSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiOslSolverInterface * oslClone = dynamic_cast<OsiOslSolverInterface*>(siClone);
      assert( oslClone != NULL );
      assert( oslClone->getModelPtr() != oslSi.getModelPtr() );
      assert( oslClone->getModelPtr() != m.getModelPtr() );
      assert( oslClone->getNumRows() == oslSi.getNumRows() );
      assert( oslClone->getNumCols() == m.getNumCols() );
      
      assert( *((int *)(oslClone->getApplicationData())) == ad );
      // Test reset
      siClone->reset();
      assert( oslClone->rowsense_==NULL );
      assert( oslClone->rhs_==NULL );
      assert( oslClone->rowrange_==NULL );
      assert( oslClone->matrixByRow_==NULL );
      assert( oslClone->ws_==NULL);
      assert( oslClone->itlimOrig_==9999999);
      delete siClone;
    }
   
    // test infinity
    {
      OsiOslSolverInterface si;
      assert( eq(si.getInfinity(),OSL_INFINITY));
    }     

    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow
    {
      OsiOslSolverInterface lhs;
      {      
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rhs_==NULL );
        assert( m.matrixByRow_==NULL );
        
        OsiOslSolverInterface siC1(m);     
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
        assert( siC1.matrixByRow_==NULL );

        const char   * siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        
        const double * siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        
        const double * siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        
        const CoinPackedMatrix * siC1mbr = siC1.getMatrixByRow();
        assert( siC1mbr != NULL );
        
        const double * ev = siC1mbr->getElements();
        assert( eq(ev[0],   3.0) );
        assert( eq(ev[1],   1.0) );
        assert( eq(ev[2],  -2.0) );
        assert( eq(ev[3],  -1.0) );
        assert( eq(ev[4],  -1.0) );
        assert( eq(ev[5],   2.0) );
        assert( eq(ev[6],   1.1) );
        assert( eq(ev[7],   1.0) );
        assert( eq(ev[8],   1.0) );
        assert( eq(ev[9],   2.8) );
        assert( eq(ev[10], -1.2) );
        assert( eq(ev[11],  5.6) );
        assert( eq(ev[12],  1.0) );
        assert( eq(ev[13],  1.9) );
        
        const CoinBigIndex * mi = siC1mbr->getVectorStarts();
        assert( mi[0]==0 );
        assert( mi[1]==5 );
        assert( mi[2]==7 );
        assert( mi[3]==9 );
        assert( mi[4]==11 );
        assert( mi[5]==14 );
        
        const int * ei = siC1mbr->getIndices();
        assert( ei[0]  ==  0 );
        assert( ei[1]  ==  1 );
        assert( ei[2]  ==  3 );
        assert( ei[3]  ==  4 );
        assert( ei[4]  ==  7 );
        assert( ei[5]  ==  1 );
        assert( ei[6]  ==  2 );
        assert( ei[7]  ==  2 );
        assert( ei[8]  ==  5 );
        assert( ei[9]  ==  3 );
        assert( ei[10] ==  6 );
        assert( ei[11] ==  0 );
        assert( ei[12] ==  4 );
        assert( ei[13] ==  7 );    
        
        assert( siC1mbr->getMajorDim() == 5 ); 
        assert( siC1mbr->getNumElements() == 14 );
        

        assert( siC1rs  == siC1.getRowSense() );
        assert( siC1rhs == siC1.getRightHandSide() );
        assert( siC1rr  == siC1.getRowRange() );

        // Change OSL Model by adding free row
        OsiRowCut rc;
        rc.setLb(-DBL_MAX);
        rc.setUb( DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);
             
        // Since model was changed, test that cached
        // data is now freed.
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
        assert( siC1.matrixByRow_==NULL );
        
        siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        assert( siC1rs[5]=='N' );

        siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        assert( eq(siC1rhs[5],0.0 ) ); 

        siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        assert( eq(siC1rr[5],0.0) );
    
        lhs=siC1;
      }
      // Test that lhs has correct values even though siC1 has gone out of scope    
      assert( lhs.rowrange_==NULL );
      assert( lhs.rowsense_==NULL );
      assert( lhs.rhs_==NULL ); 
      assert( lhs.matrixByRow_==NULL ); 
      
      const char * lhsrs  = lhs.getRowSense();
      assert( lhsrs[0]=='G' );
      assert( lhsrs[1]=='L' );
      assert( lhsrs[2]=='E' );
      assert( lhsrs[3]=='R' );
      assert( lhsrs[4]=='R' );
      assert( lhsrs[5]=='N' );
      
      const double * lhsrhs = lhs.getRightHandSide();
      assert( eq(lhsrhs[0],2.5) );
      assert( eq(lhsrhs[1],2.1) );
      assert( eq(lhsrhs[2],4.0) );
      assert( eq(lhsrhs[3],5.0) );
      assert( eq(lhsrhs[4],15.) ); 
      assert( eq(lhsrhs[5],0.0) ); 
      
      const double *lhsrr  = lhs.getRowRange();
      assert( eq(lhsrr[0],0.0) );
      assert( eq(lhsrr[1],0.0) );
      assert( eq(lhsrr[2],0.0) );
      assert( eq(lhsrr[3],5.0-1.8) );
      assert( eq(lhsrr[4],15.0-3.0) );
      assert( eq(lhsrr[5],0.0) );      
      
      const CoinPackedMatrix * lhsmbr = lhs.getMatrixByRow();
      assert( lhsmbr != NULL );       
      const double * ev = lhsmbr->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const CoinBigIndex * mi = lhsmbr->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = lhsmbr->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    
      
      int md = lhsmbr->getMajorDim();
      assert(  md == 6 ); 
      assert( lhsmbr->getNumElements() == 14 );
    }
    
    assert(OsiOslSolverInterface::getNumInstances()==1);
  }
  assert(OsiOslSolverInterface::getNumInstances()==0);

  // Do common solverInterface testing 
  {
    OsiOslSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir, netlibDir );
  }

  // Do primitive branch and bound
  // And test message handling
  /* This could be moved down to OsiSolverInterface by
     putting in branchAndBound and taking off m.
  */
  int iPass;
  for (iPass=0;iPass<2;iPass++) {
    OsiOslSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
     // derived message handler (only used on second pass)
    OsiOslMessageTest messageHandler;
   if (iPass) {
      std::cout<<"Testing derived message handler"<<std::endl;
      m.passInMessageHandler(&messageHandler);
      ekk_registerMsguCallBack(m.getModelPtr(),trapMessages);
      // say all messages to be trapped
      ekk_mset(m.getModelPtr(),1,0,-1,2,9999,0);
      // pass handler to OSL
      ekk_setUserData(m.getModelPtr(),&messageHandler);
   }
    // solve LP
    m.initialSolve();
    // setColBounds prints every time - don't even get to message handler
    ekk_messagePrintOff(m.getModelPtr() ,317);
    ekk_messagePrintOff(m.getModelPtr() ,318);
    ekk_messagePrintOff(m.getModelPtr() ,3048);
    ekk_messagePrintOff(m.getModelPtr() ,85);
    ekk_messagePrintOff(m.getModelPtr() ,82);
    ekk_messagePrintOff(m.getModelPtr() ,38);

    if (m.isProvenOptimal()&&!m.isDualObjectiveLimitReached()) {
      // This is a really simple Branch and Bound code - mainly
      // to test strong branching
      // I should look at STL more to allow other than depth first
      int numberIntegers=0;
      int numberColumns = m.getNumCols();
      int iColumn;
      int i;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if( m.isInteger(iColumn))
	  numberIntegers++;
      }
      if (!numberIntegers) {
	std::cout<<"No integer variables"
	  <<std::endl;
	return;
      }
      int * which = new int[numberIntegers]; // which variables are integer
      numberIntegers=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
	if( m.isInteger(iColumn))
	  which[numberIntegers++]=iColumn;
      }
      
      // empty tree
      OsiVectorNode branchingTree;
      
      // Add continuous to it;
      branchingTree.push_back(OsiNodeSimple(m,numberIntegers,which));
      
      // For printing totals
      int numberIterations=0;
      int numberNodes =0;
      
      OsiNodeSimple bestNode;
      // while until nothing on stack
      while (branchingTree.size()) {
	// last node
	OsiNodeSimple node = branchingTree.back();
	branchingTree.pop_back();
	numberNodes++;
	if (node.variable_>=0) {
	  // branch - do bounds
	  for (i=0;i<numberIntegers;i++) {
	    iColumn=which[i];
	    m.setColBounds( iColumn,node.lower_[i],node.upper_[i]);
	  }
	  // move basis
	  m.setWarmStart(&node.basis_);
	  // do branching variable
	  if (node.way_<0) {
	    m.setColUpper(which[node.variable_],floor(node.value_));
	    // now push back node if more to come
	    if (node.way_==-1) { 
	      node.way_=+2;	  // Swap direction
	      branchingTree.push_back(node);
	    }
	  } else {
	    m.setColLower(which[node.variable_],ceil(node.value_));
	    // now push back node if more to come
	    if (node.way_==1) { 
	      node.way_=-2;	  // Swap direction
	      branchingTree.push_back(node);
	    }
	  }
	  // solve
	  m.resolve();
	  numberIterations += m.getIterationCount();
	  if (!m.isIterationLimitReached()) {
	    OsiNodeSimple newNode(m,numberIntegers,which);
	    // something extra may have been fixed by strong branching
	    // if so go round again
	    while (newNode.variable_==numberIntegers) {
	      m.resolve();
	      newNode = OsiNodeSimple(m,numberIntegers,which);
	    }
	    if (newNode.objectiveValue_<1.0e100) {
	      // push on stack
	      branchingTree.push_back(newNode);
	    }
	  } else {
	    // maximum iterations - exit
	    std::cout<<"Exiting on maximum iterations"
	      <<std::endl;
	  break;
	  }
	} else {
	  // integer solution - save
	  bestNode = node;
	  // set cutoff (hard coded tolerance)
	  m.setDblParam(OsiDualObjectiveLimit,bestNode.objectiveValue_-1.0e-5);
	  std::cout<<"Integer solution of "
		   <<bestNode.objectiveValue_
		   <<" found after "<<numberIterations
		   <<" iterations and "<<numberNodes<<" nodes"
		   <<std::endl;
	}
      }
      std::cout<<"Search took "
	       <<numberIterations
	       <<" iterations and "<<numberNodes<<" nodes"
	       <<std::endl;
      if (bestNode.numberIntegers_) {
	// we have a solution restore
	// do bounds
	for (i=0;i<numberIntegers;i++) {
	  iColumn=which[i];
	  m.setColBounds( iColumn,bestNode.lower_[i],bestNode.upper_[i]);
	}
	// move basis
	m.setWarmStart(&bestNode.basis_);
	m.resolve();
      }
      delete [] which;
    } else {
      std::cout<<"The LP relaxation is infeasible"
	       <<std::endl;
      throw CoinError("The LP relaxation is infeasible or too expensive",
		      "branchAndBound", "OsiClpSolverInterface");
    }
    ekk_messagePrintOn(m.getModelPtr() ,317);
    ekk_messagePrintOn(m.getModelPtr() ,318);
    ekk_messagePrintOn(m.getModelPtr() ,3048);
    ekk_messagePrintOn(m.getModelPtr() ,85);
    ekk_messagePrintOn(m.getModelPtr() ,82);
    ekk_messagePrintOn(m.getModelPtr() ,38);
    // Very important to normalize before going out of scope
    ekk_clearMsguCallBack(m.getModelPtr());
    ekk_setUserData(m.getModelPtr(),NULL);
   }

  assert(OsiOslSolverInterface::getNumInstances()==0);

}
#endif
