// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiOslSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiOsiMessage.hpp"

// below needed for pathetic branch and bound code
#include <vector>
#include <map>
using std::max;
using std::min;
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

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif
class OsiOslMessageTest :
   public OsiMessageHandler {

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
	return OsiMessageHandler::print(); // print
    } else if (currentMessage().externalNumber()==3000) { 
      numberInfeasible++;
    }
  }  else if (currentSource()=="Osi") {
    if (currentMessage().externalNumber()==5) 
      std::cout<<"End of search trapped, "
	       <<numberOptimal<<" optimal Lp's, "
	       <<numberInfeasible
	       <<" infeasible (including strong branching)"<<std::endl;
    return OsiMessageHandler::print();
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
    OsiMessageHandler * handler =  (OsiMessageHandler *) ekk_userData(model);
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
OsiOslSolverInterfaceUnitTest(const std::string & mpsDir)
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
    OsiRelFltEq eq;
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
      delete siClone;
    }
   
    // test infinity
    {
      OsiOslSolverInterface si;
      assert( eq(si.getInfinity(),OSL_INFINITY));
    }     
    
    // Test setting solution
    {
      OsiOslSolverInterface m1(m);
      int i;

      double * cs = new double[m1.getNumCols()];
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        cs[i] = i + .5;
      m1.setColSolution(cs);
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        assert(m1.getColSolution()[i] == i + .5);
      
      double * rs = new double[m1.getNumRows()];
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        rs[i] = i - .5;
      m1.setRowPrice(rs);
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        assert(m1.getRowPrice()[i] == i - .5);

      delete [] cs;
      delete [] rs;
    }
    
    
    // Test fraction Indices
    {
      OsiOslSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      //fim.setModelPtr(m);
      // exmip1.mps has 2 integer variables with index 2 & 3
      assert(  fim.isContinuous(0) );
      assert(  fim.isContinuous(1) );
      assert( !fim.isContinuous(2) );
      assert( !fim.isContinuous(3) );
      assert(  fim.isContinuous(4) );
      
      assert( !fim.isInteger(0) );
      assert( !fim.isInteger(1) );
      assert(  fim.isInteger(2) );
      assert(  fim.isInteger(3) );
      assert( !fim.isInteger(4) );
      
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert(  fim.isBinary(2) );
      assert(  fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert( !fim.isIntegerNonBinary(2) );
      assert( !fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );

      // Test fractionalIndices
      {
        double sol[]={2.9,3.0};
        ekk_copyColsol(fim.getModelPtr(),sol,2,4);
        OsiVectorInt fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        sol[0]=5 + .00001/2.;
        sol[1]=8 - .00001/2.;
        ekk_copyColsol(fim.getModelPtr(),sol,2,4);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        sol[0]=5 + .00001*2.;
        sol[1]=8 - .00001*2.;
        ekk_copyColsol(fim.getModelPtr(),sol,2,4);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );
      }


      
      // Change date so column 2 & 3 are integerNonBinary
      double ub[]={5.0,6.0};
      ekk_copyColupper(fim.getModelPtr(),ub,2,4);
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert( !fim.isBinary(2) );
      assert( !fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert(  fim.isIntegerNonBinary(2) );
      assert(  fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );
    }
    
    // Test apply cuts method
    {      
      OsiOslSolverInterface im(m);
      OsiCuts cuts;
      
      // Generate some cuts 
      {
        // Get number of rows and columns in model
        int nr=im.getNumRows();
        int nc=im.getNumCols();
        assert( nr == 5 );
        assert( nc == 8 );
        
        // Generate a valid row cut from thin air
        int c;
        {
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *el = new double[nc];
          for (c=0;c<nc;c++) el[c]=((double)c)*((double)c);
          
          OsiRowCut rc;
          rc.setRow(nc,inx,el);
          rc.setLb(-100.);
          rc.setUb(100.);
          rc.setEffectiveness(22);
          
          cuts.insert(rc);
          delete[]el;
          delete[]inx;
        }
        
        // Generate valid col cut from thin air
        {
          const double * oslColLB = im.getColLower();
          const double * oslColUB = im.getColUpper();
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *lb = new double[nc];
          double *ub = new double[nc];
          for (c=0;c<nc;c++) lb[c]=oslColLB[c]+0.001;
          for (c=0;c<nc;c++) ub[c]=oslColUB[c]-0.001;
          
          OsiColCut cc;
          cc.setLbs(nc,inx,lb);
          cc.setUbs(nc,inx,ub);
          
          cuts.insert(cc);
          delete [] ub;
          delete [] lb;
          delete [] inx;
        }
        
        {
          // Generate a row and column cut which are ineffective
          OsiRowCut * rcP= new OsiRowCut;
          rcP->setEffectiveness(-1.);
          cuts.insert(rcP);
          assert(rcP==NULL);
          
          OsiColCut * ccP= new OsiColCut;
          ccP->setEffectiveness(-12.);
          cuts.insert(ccP);
          assert(ccP==NULL);
        }
        {
          //Generate inconsistent Row cut
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          rc.setLb(3.);
          rc.setUb(4.);
          assert(!rc.consistent());
          cuts.insert(rc);
        }
        {
          //Generate inconsistent col cut
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          cc.setUbs(ne,inx,el);
          assert(!cc.consistent());
          cuts.insert(cc);
        }
        {
          // Generate row cut which is inconsistent for model m
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          assert(rc.consistent());
          assert(!rc.consistent(im));
          cuts.insert(rc);
        }
        {
          // Generate col cut which is inconsistent for model m
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={30};
          double el[ne]={2.0};
          cc.setLbs(ne,inx,el);
          assert(cc.consistent());
          assert(!cc.consistent(im));
          cuts.insert(cc);
        }
        {
          // Generate col cut which is infeasible
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={0};
          double el[ne]={2.0};
          cc.setUbs(ne,inx,el);
          cc.setEffectiveness(1000.);
          assert(cc.consistent());
          assert(cc.consistent(im));
          assert(cc.infeasible(im));
          cuts.insert(cc);
        }
      }
      assert(cuts.sizeRowCuts()==4);
      assert(cuts.sizeColCuts()==5);
      
      OsiSolverInterface::ApplyCutsReturnCode rc = im.applyCuts(cuts);
      assert( rc.getNumIneffective() == 2 );
      assert( rc.getNumApplied() == 2 );
      assert( rc.getNumInfeasible() == 1 );
      assert( rc.getNumInconsistentWrtIntegerModel() == 2 );
      assert( rc.getNumInconsistent() == 2 );
      assert( cuts.sizeCuts() == rc.getNumIneffective() +
        rc.getNumApplied() +
        rc.getNumInfeasible() +
        rc.getNumInconsistentWrtIntegerModel() +
        rc.getNumInconsistent() );
    }
    
    {    
      OsiOslSolverInterface oslSi(m);
      int nc = oslSi.getNumCols();
      int nr = oslSi.getNumRows();
      const double * cl = oslSi.getColLower();
      const double * cu = oslSi.getColUpper();
      const double * rl = oslSi.getRowLower();
      const double * ru = oslSi.getRowUpper();
      assert( nc == 8 );
      assert( nr == 5 );
      assert( eq(cl[0],2.5) );
      assert( eq(cl[1],0.0) );
      assert( eq(cu[1],4.1) );
      assert( eq(cu[2],1.0) );
      assert( eq(rl[0],2.5) );
      assert( eq(rl[4],3.0) );
      assert( eq(ru[1],2.1) );
      assert( eq(ru[4],15.0) );
      
      const double * cs = oslSi.getColSolution();
      assert( eq(cs[0],2.5) );
      assert( eq(cs[7],0.0) );
      
      assert( !eq(cl[3],1.2345) );
      oslSi.setColLower( 3, 1.2345 );
      assert( eq(oslSi.getColLower()[3],1.2345) );
      
      assert( !eq(cu[4],10.2345) );
      oslSi.setColUpper( 4, 10.2345 );
      assert( eq(oslSi.getColUpper()[4],10.2345) );

      assert( eq(oslSi.getObjValue(),0.0) );

      assert( eq( oslSi.getObjCoefficients()[0],  1.0) );
      assert( eq( oslSi.getObjCoefficients()[1],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[2],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[3],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[4],  2.0) );
      assert( eq( oslSi.getObjCoefficients()[5],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[6],  0.0) );
      assert( eq( oslSi.getObjCoefficients()[7], -1.0) );
    }
    
    // Test matrixByRow method
    { 
      const OsiOslSolverInterface si(m);
      const OsiPackedMatrix * smP = si.getMatrixByRow();
      // LL:      const OsiOslPackedMatrix * osmP = dynamic_cast<const OsiOslPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      OsiRelFltEq eq;
      const double * ev = smP->getElements();
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
      
      const int * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = smP->getIndices();
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
      
      assert( smP->getMajorDim() == 5 ); 
      assert( smP->getNumElements() == 14 );
      
    }
        // Test matrixByCol method
    {
  
      const OsiOslSolverInterface si(m);
      const OsiPackedMatrix * smP = si.getMatrixByCol();
      // LL:      const OsiOslPackedMatrix * osmP = dynamic_cast<const OsiOslPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      OsiRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   5.6) );
      assert( eq(ev[2],   1.0) );
      assert( eq(ev[3],   2.0) );
      assert( eq(ev[4],   1.1) );
      assert( eq(ev[5],   1.0) );
      assert( eq(ev[6],  -2.0) );
      assert( eq(ev[7],   2.8) );
      assert( eq(ev[8],  -1.0) );
      assert( eq(ev[9],   1.0) );
      assert( eq(ev[10],  1.0) );
      assert( eq(ev[11], -1.2) );
      assert( eq(ev[12], -1.0) );
      assert( eq(ev[13],  1.9) );
      
      const int * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==2 );
      assert( mi[2]==4 );
      assert( mi[3]==6 );
      assert( mi[4]==8 );
      assert( mi[5]==10 );
      assert( mi[6]==11 );
      assert( mi[7]==12 );
      assert( mi[8]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  4 );
      assert( ei[2]  ==  0 );
      assert( ei[3]  ==  1 );
      assert( ei[4]  ==  1 );
      assert( ei[5]  ==  2 );
      assert( ei[6]  ==  0 );
      assert( ei[7]  ==  3 );
      assert( ei[8]  ==  0 );
      assert( ei[9]  ==  4 );
      assert( ei[10] ==  2 );
      assert( ei[11] ==  3 );
      assert( ei[12] ==  0 );
      assert( ei[13] ==  4 );    
      
      assert( smP->getMajorDim() == 8 ); 
      assert( smP->getNumElements() == 14 );

      assert( smP->getSizeVectorStarts()==9 );
      assert( smP->getMinorDim() == 5 );
      
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
        
        const OsiPackedMatrix * siC1mbr = siC1.getMatrixByRow();
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
        
        const int * mi = siC1mbr->getVectorStarts();
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
      
      const OsiPackedMatrix * lhsmbr = lhs.getMatrixByRow();
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
      
      const int * mi = lhsmbr->getVectorStarts();
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
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir);
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
    OsiMessageHandler * handler = m.messageHandler();
    OsiMessages messages = m.messages() ;
    if (!iPass) {
      // switch on strong branching  message
      int oneMessage=7;
      // we want to change models copy of messages
      m.messagesPointer()->setDetailMessages(1,1,&oneMessage);
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
	handler->message(OSI_BAB_NOINT,messages)
	  <<OsiMessageEol;
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
      branchingTree.push_back(OsiNode(m,numberIntegers,which));
      
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
	    OsiNode newNode(m,numberIntegers,which);
	    // something extra may have been fixed by strong branching
	    // if so go round again
	    while (newNode.variable_==numberIntegers) {
	      m.resolve();
	      newNode = OsiNode(m,numberIntegers,which);
	    }
	    if (newNode.objectiveValue_<1.0e100) {
	      // push on stack
	      branchingTree.push_back(newNode);
	    }
	  } else {
	    // maximum iterations - exit
	    handler->message(OSI_BAB_MAXITS,messages)
	      <<OsiMessageEol;
	  break;
	  }
	} else {
	  // integer solution - save
	  bestNode = node;
	  // set cutoff (hard coded tolerance)
	  m.setDblParam(OsiDualObjectiveLimit,bestNode.objectiveValue_-1.0e-5);
	  handler->message(OSI_BAB_SOLUTION,messages)
	    <<bestNode.objectiveValue_<<numberIterations
	    <<numberNodes
	    <<OsiMessageEol;
	}
      }
      handler->message(OSI_BAB_END,messages)
	<<numberIterations
	<<numberNodes
	<<OsiMessageEol;
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
      handler->message(OSI_BAB_INFEAS,messages)
	<<OsiMessageEol;
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
