// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
  
#include <cassert>

#include "OsiClpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinMessage.hpp"
#include "ClpMessage.hpp"
#include "ClpFactorization.hpp"

//#############################################################################

#ifdef NDEBUG
#undef NDEBUG
#endif
class OsiClpMessageTest :
   public CoinMessageHandler {

public:
  virtual int print() ;
  OsiClpMessageTest();
};

OsiClpMessageTest::OsiClpMessageTest() : CoinMessageHandler()
{
}
int
OsiClpMessageTest::print()
{
  if (currentMessage().externalNumber()==0&&currentSource()=="Clp") 
    std::cout<<"This is not actually an advertisement by Dash Associates - just my feeble attempt to test message handling and language - JJHF"<<std::endl;
  else if (currentMessage().externalNumber()==5&&currentSource()=="Osi") 
    std::cout<<"End of search trapped"<<std::endl;
  return CoinMessageHandler::print();
}

//--------------------------------------------------------------------------
// test EKKsolution methods.
void
OsiClpSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{
  
  // Test default constructor
  {
    OsiClpSolverInterface m;
    assert( m.rowsense_==NULL );
    assert( m.rhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.ws_==NULL);
    assert( m.itlimOrig_==9999999);
    assert( m.lastAlgorithm_==0);
    assert( m.integerInformation_==NULL);
  }
  
  
  {    
    CoinRelFltEq eq;
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    
    {
      OsiClpSolverInterface im;    
      
      assert( im.getNumCols() == 0 ); 
      
      assert( im.getModelPtr()!=NULL );
    }
    
    // Test copy constructor and assignment operator
    {
      OsiClpSolverInterface lhs;
      {      
        OsiClpSolverInterface im(m);        
	
        OsiClpSolverInterface imC1(im);
        assert( imC1.getModelPtr()!=im.getModelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        
        OsiClpSolverInterface imC2(im);
        assert( imC2.getModelPtr()!=im.getModelPtr() );
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
	
        assert( imC2.getModelPtr()!=imC1.getModelPtr() );
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope
      
      assert( lhs.getModelPtr() != m.getModelPtr() );
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
    }
    // Test clone
    {
      OsiClpSolverInterface oslSi(m);
      OsiSolverInterface * siPtr = &oslSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiClpSolverInterface * oslClone = dynamic_cast<OsiClpSolverInterface*>(siClone);
      assert( oslClone != NULL );
      assert( oslClone->getModelPtr() != oslSi.getModelPtr() );
      assert( oslClone->getNumRows() == oslSi.getNumRows() );
      assert( oslClone->getNumCols() == m.getNumCols() );
      
      delete siClone;
    }
  
    // test infinity
    {
      OsiClpSolverInterface si;
      assert( eq(si.getInfinity(),OsiClpInfinity));
    }     
    
    // Test setting solution
    {
      OsiClpSolverInterface m1(m);
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
      OsiClpSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
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
	memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        OsiVectorInt fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        sol[0]=5 + .00001/2.;
        sol[1]=8 - .00001/2.;
	memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        sol[0]=5 + .00001*2.;
        sol[1]=8 - .00001*2.;
	memcpy(fim.modelPtr_->primalColumnSolution()+2,sol,2*sizeof(double));
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );
      }


      
      // Change date so column 2 & 3 are integerNonBinary
      double ub[]={5.0,6.0};
      memcpy(fim.modelPtr_->columnUpper()+2,ub,2*sizeof(double));
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
      OsiClpSolverInterface im(m);
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
      OsiClpSolverInterface oslSi(m);
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

      double objValue = oslSi.getObjValue();
      assert( eq(objValue,3.5) );

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
      const OsiClpSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByRow();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
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
    // Test adding several cuts
    {
      OsiClpSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      // exmip1.mps has 2 integer variables with index 2 & 3
      fim.initialSolve();
      OsiRowCut cuts[3];
      
      
      // Generate one ineffective cut plus two trivial cuts
      int c;
      int nc = fim.getNumCols();
      int *inx = new int[nc];
      for (c=0;c<nc;c++) inx[c]=c;
      double *el = new double[nc];
      for (c=0;c<nc;c++) el[c]=1.0e-50+((double)c)*((double)c);
      
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      el[4]=0.0; // to get inf later
      
      for (c=2;c<4;c++) {
	el[0]=1.0;
	inx[0]=c;
	cuts[c-1].setRow(1,inx,el);
	cuts[c-1].setLb(1.);
	cuts[c-1].setUb(100.);
	cuts[c-1].setEffectiveness(c);
      }
      fim.writeMps("x1.mps");
      fim.applyRowCuts(3,cuts);
      fim.writeMps("x2.mps");
      // resolve - should get message about zero elements
      fim.resolve();
      fim.writeMps("x3.mps");
      // check integer solution
      const double * cs = fim.getColSolution();
      CoinRelFltEq eq;
      assert( eq(cs[2],   1.0) );
      assert( eq(cs[3],   1.0) );
      // check will find invalid matrix
      el[0]=1.0/el[4];
      inx[0]=0;
      cuts[0].setRow(nc,inx,el);
      cuts[0].setLb(-100.);
      cuts[0].setUb(500.);
      cuts[0].setEffectiveness(22);
      fim.applyRowCut(cuts[0]);
      // resolve - should get message about zero elements
      fim.resolve();
      assert (fim.isAbandoned());
      delete[]el;
      delete[]inx;
    }
        // Test matrixByCol method
    {
  
      const OsiClpSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByCol();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
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
      OsiClpSolverInterface lhs;
      {      
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rhs_==NULL );
        assert( m.matrixByRow_==NULL );
        
        OsiClpSolverInterface siC1(m);     
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
    
  }

  // Solve an integer problem 
#if 0
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    std::cout<<"Testing languages and derived message handler"<<std::endl;
    m.setLanguage(CoinMessages::uk_en);
    // derived message handler
    OsiClpMessageTest messageHandler;
    // Vanilla one but to FILE ptr

    CoinMessageHandler defaultHandler(stderr);
    
    /* This is a bit of a mess - or it could be a feature.  Osi and
       Clp can be using different message handlers.  On reflection I
       think this makes sense */
    // First just trap Clp messages and see when optimal
    m.getModelPtr()->passInMessageHandler(&messageHandler);
    
    m.initialSolve();
    m.setLanguage(CoinMessages::us_en);
    // now put Clp messages to stderr and trap end of search
    m.getModelPtr()->passInMessageHandler(&defaultHandler);

    m.passInMessageHandler(&messageHandler);
    // This also tests slow strong branching
    m.branchAndBound();
    m.resolve();
    //m.getObjValue();
  }
#endif

  // Solve an lp by hand
  {    
    OsiClpSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    m.setObjSense(-1.0);
    m.getModelPtr()->messageHandler()->setLogLevel(4);
    m.initialSolve();
    m.getModelPtr()->factorization()->maximumPivots(5);
    m.setObjSense(1.0);
    // enable special mode
    m.enableSimplexInterface(true);
    // we happen to know that variables are 0-1 and rows are L
    int numberIterations=0;
    int numberColumns = m.getNumCols();
    int numberRows = m.getNumRows();
    double * fakeCost = new double[numberColumns];
    double * duals = new double [numberRows];
    double * djs = new double [numberColumns];
    const double * solution = m.getColSolution();
    memcpy(fakeCost,m.getObjCoefficients(),numberColumns*sizeof(double));
    while (1) {
      const double * dj;
      const double * dual;
      if ((numberIterations&1)==0) {
	// use given ones
	dj = m.getReducedCost();
	dual = m.getRowPrice();
      } else {
	// create
	dj = djs;
	dual = duals;
	m.getReducedGradient(djs,duals,fakeCost);
      }
      int i;
      int colIn=9999;
      int direction=1;
      double best=1.0e-6;
      // find most negative reduced cost
      // Should check basic - but should be okay on this problem
      for (i=0;i<numberRows;i++) {
	double value=dual[i];
	if (value>best) {
	  direction=-1;
	  best=value;
	  colIn=-i-1;
	}
      }
      for (i=0;i<numberColumns;i++) {
	double value=dj[i];
	if (value<-best&&solution[i]<1.0e-6) {
	  direction=1;
	  best=-value;
	  colIn=i;
	} else if (value>best&&solution[i]>1.0-1.0e-6) {
	  direction=-1;
	  best=value;
	  colIn=i;
	}
      }
      if (colIn==9999)
	break; // should be optimal
      int colOut;
      int outStatus;
      double theta;
      assert(!m.primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL));
      printf("out %d, direction %d theta %g\n",
	     colOut,outStatus,theta);
      numberIterations++;
    }
    delete [] fakeCost;
    delete [] duals;
    delete [] djs;
    // exit special mode
    m.disableSimplexInterface();
    m.getModelPtr()->messageHandler()->setLogLevel(4);
    m.resolve();
    assert (!m.getIterationCount());
    m.setObjSense(-1.0);
    m.initialSolve();
  }
  // Do common solverInterface testing 
  {
    OsiClpSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
}
