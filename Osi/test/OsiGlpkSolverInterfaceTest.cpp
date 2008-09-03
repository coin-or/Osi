// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// Copyright (C) 2004  University of Pittsburgh
//   University of Pittsburgh coding done by Brady Hunsaker

// OsiGlpkSolverInterfaceTest.cpp adapted from OsiClpSolverInterfaceTest.cpp 
//  on 2004/10/16
// Check for ??? to see tests that aren't working correctly.
// Also note that OsiPresolve doesn't appear to work with OsiGlpk.  This
// needs to be examined.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
  
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "OsiConfig.h"

#include <cassert>
//#include <cstdlib>
//#include <cstdio>
//#include <iostream>

#include "OsiGlpkSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinMessage.hpp"
#include "CoinFinite.hpp"

// Added so windows build with dsp files works,
// when not building with glpk.
#ifdef COIN_HAS_GLPK

//#############################################################################

/*
  Define helper routines in the file-local namespace. Makes it easier to see
  the main flow of tests.
*/


//--------------------------------------------------------------------------
int
OsiGlpkSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{
	int errCnt = 0;
  // Do common solverInterface testing 
  {
    OsiGlpkSolverInterface m;
    errCnt += OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
  
  // Test default constructor
  {
    OsiGlpkSolverInterface m;
    assert( m.obj_==NULL );
    assert( m.collower_==NULL );
    assert( m.colupper_==NULL );
    assert( m.ctype_==NULL );
    assert( m.rowsense_==NULL );
    assert( m.rhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.rowlower_==NULL );
    assert( m.rowupper_==NULL );
    assert( m.colsol_==NULL );
    assert( m.rowsol_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.matrixByCol_==NULL );
    assert( m.getApplicationData() == NULL );
    int i=2346;
    m.setApplicationData(&i);
    assert( *((int *)(m.getApplicationData())) == i );
  }
  
  
  {    
    CoinRelFltEq eq;
    OsiGlpkSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    
    {
      OsiGlpkSolverInterface im;    
      
      assert( im.getNumCols() == 0 ); 
      
      assert( im.getModelPtr()!=NULL );
      // Test reset
      im.reset();
      assert( im.obj_==NULL );
      assert( im.collower_==NULL );
      assert( im.colupper_==NULL );
      assert( im.ctype_==NULL );
      assert( im.rowsense_==NULL );
      assert( im.rhs_==NULL );
      assert( im.rowrange_==NULL );
      assert( im.rowlower_==NULL );
      assert( im.rowupper_==NULL );
      assert( im.colsol_==NULL );
      assert( im.rowsol_==NULL );
      assert( im.matrixByRow_==NULL );
      assert( im.matrixByCol_==NULL );
      assert( im.getApplicationData() == NULL );
    }
    
    // Test copy constructor and assignment operator
    {
      OsiGlpkSolverInterface lhs;
      {      
        OsiGlpkSolverInterface im(m);        
	
        OsiGlpkSolverInterface imC1(im);
        assert( imC1.getModelPtr()!=im.getModelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        
        OsiGlpkSolverInterface imC2(im);
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
      OsiGlpkSolverInterface oslSi(m);
      OsiSolverInterface * siPtr = &oslSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiGlpkSolverInterface * oslClone = dynamic_cast<OsiGlpkSolverInterface*>(siClone);
      assert( oslClone != NULL );
      assert( oslClone->getModelPtr() != oslSi.getModelPtr() );
      assert( oslClone->getNumRows() == oslSi.getNumRows() );
      assert( oslClone->getNumCols() == m.getNumCols() );
      
      delete siClone;
    }
  
    // test infinity
    {
      OsiGlpkSolverInterface si;
      assert( eq(si.getInfinity(),COIN_DBL_MAX) );
    }     
    
    // Test fraction Indices
    {
      OsiGlpkSolverInterface fim;
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
	// Set a solution vector
	double * sol = new double[fim.getNumCols()];
	for ( int i = 0;  i < fim.getNumCols();  sol[i++] = 0.0 );
	sol[2] = 2.9;
	sol[3] = 3.0;
	fim.setColSolution(sol);

        OsiVectorInt fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        sol[2]=5 + .00001/2.;
        sol[3]=8 - .00001/2.;
        fim.setColSolution(sol);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        sol[2]=5 + .00001*2.;
        sol[3]=8 - .00001*2.;
        fim.setColSolution(sol);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );

	delete [] sol;
      }


      
      // Change data so column 2 & 3 are integerNonBinary
      fim.setColUpper(2, 5);
      fim.setColUpper(3, 6.0);
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
#if 0
    // ??? These index error 'throw's aren't in OsiGlpk 

    // Test some catches
    {
      OsiGlpkSolverInterface solver;
      try {
	solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
	std::cout<<"Correct throw"<<std::endl;
      }
      std::string fn = mpsDir+"exmip1";
      solver.readMps(fn.c_str(),"mps");
      try {
	solver.setObjCoeff(0,0.0);
      }
      catch (CoinError e) {
	std::cout<<"** Incorrect throw"<<std::endl;
	abort();
      }
      try {
	int index[]={0,20};
	double value[]={0.0,0.0,0.0,0.0};
	solver.setColSetBounds(index,index+2,value);
      }
      catch (CoinError e) {
	std::cout<<"Correct throw"<<std::endl;
      }
    }
#endif
    // Test apply cuts method
    {      
      OsiGlpkSolverInterface im(m);
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
      OsiGlpkSolverInterface oslSi(m);
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
      
      // NO - cu has been refreshed assert( !eq(cu[4],10.2345) );
      assert( !eq(oslSi.getColUpper()[4],10.2345) );
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
      const OsiGlpkSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByRow();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
      const double * ev = smP->getElements();

      // GLPK returns each row in reverse order. This is consistent with 
      // the sparse matrix format but is not what most solvers do.  That's
      // why this section is different.
#if 1
      assert( eq(ev[4],   3.0) );
      assert( eq(ev[3],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[1],  -1.0) );
      assert( eq(ev[0],  -1.0) );
      assert( eq(ev[6],   2.0) );
      assert( eq(ev[5],   1.1) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[10],   2.8) );
      assert( eq(ev[9], -1.2) );
      assert( eq(ev[13],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[11],  1.9) );
#else  // this is how it would "normally" be
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
#endif
      
      const CoinBigIndex * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = smP->getIndices();
#if 1  // GLPK's way
      assert( ei[4]  ==  0 );
      assert( ei[3]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[1]  ==  4 );
      assert( ei[0]  ==  7 );
      assert( ei[6]  ==  1 );
      assert( ei[5]  ==  2 );
      assert( ei[8]  ==  2 );
      assert( ei[7]  ==  5 );
      assert( ei[10]  ==  3 );
      assert( ei[9] ==  6 );
      assert( ei[13] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[11] ==  7 );    
#else  // the "normal" way
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
#endif
      
      assert( smP->getMajorDim() == 5 ); 
      assert( smP->getNumElements() == 14 );
      
    }
    // Test adding several cuts
    {
      OsiGlpkSolverInterface fim;
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
#if 0  // ??? Not working for some reason.
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
#endif
      delete[]el;
      delete[]inx;
    }
        // Test matrixByCol method
    {
  
      const OsiGlpkSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByCol();
      // LL:      const OsiClpPackedMatrix * osmP = dynamic_cast<const OsiClpPackedMatrix*>(smP);
      // LL: assert( osmP!=NULL );
      
      CoinRelFltEq eq;
      const double * ev = smP->getElements();
      // Unlike row-ordered matrices, GLPK does column-ordered the "normal" way
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
      
      const CoinBigIndex * mi = smP->getVectorStarts();
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
      OsiGlpkSolverInterface lhs;
      {      
#if 0  // ??? these won't work because the copy constructor changes the
	// values in m 
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rhs_==NULL );
        assert( m.matrixByRow_==NULL );
#endif
        
        OsiGlpkSolverInterface siC1(m);     
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
#if 1 // GLPK's way
        assert( eq(ev[4],   3.0) );
        assert( eq(ev[3],   1.0) );
        assert( eq(ev[2],  -2.0) );
        assert( eq(ev[1],  -1.0) );
        assert( eq(ev[0],  -1.0) );
        assert( eq(ev[6],   2.0) );
        assert( eq(ev[5],   1.1) );
        assert( eq(ev[8],   1.0) );
        assert( eq(ev[7],   1.0) );
        assert( eq(ev[10],   2.8) );
        assert( eq(ev[9], -1.2) );
        assert( eq(ev[13],  5.6) );
        assert( eq(ev[12],  1.0) );
        assert( eq(ev[11],  1.9) );
#else // the "normal" way
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
#endif
        
        const CoinBigIndex * mi = siC1mbr->getVectorStarts();
        assert( mi[0]==0 );
        assert( mi[1]==5 );
        assert( mi[2]==7 );
        assert( mi[3]==9 );
        assert( mi[4]==11 );
        assert( mi[5]==14 );
        
        const int * ei = siC1mbr->getIndices();
#if 1  // GLPK's way
        assert( ei[4]  ==  0 );
        assert( ei[3]  ==  1 );
        assert( ei[2]  ==  3 );
        assert( ei[1]  ==  4 );
        assert( ei[0]  ==  7 );
        assert( ei[6]  ==  1 );
        assert( ei[5]  ==  2 );
        assert( ei[8]  ==  2 );
        assert( ei[7]  ==  5 );
        assert( ei[10]  ==  3 );
        assert( ei[9] ==  6 );
        assert( ei[13] ==  0 );
        assert( ei[12] ==  4 );
        assert( ei[11] ==  7 );    
#else  // the "normal" way
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
#endif
        
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
#if 1  // GLPK's way
      assert( eq(ev[4],   3.0) );
      assert( eq(ev[3],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[1],  -1.0) );
      assert( eq(ev[0],  -1.0) );
      assert( eq(ev[6],   2.0) );
      assert( eq(ev[5],   1.1) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[10],   2.8) );
      assert( eq(ev[9], -1.2) );
      assert( eq(ev[13],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[11],  1.9) );
#else // the "normal" way
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
#endif
      
      const CoinBigIndex * mi = lhsmbr->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = lhsmbr->getIndices();
#if 1  // GLPK's way
      assert( ei[4]  ==  0 );
      assert( ei[3]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[1]  ==  4 );
      assert( ei[0]  ==  7 );
      assert( ei[6]  ==  1 );
      assert( ei[5]  ==  2 );
      assert( ei[8]  ==  2 );
      assert( ei[7]  ==  5 );
      assert( ei[10]  ==  3 );
      assert( ei[9] ==  6 );
      assert( ei[13] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[11] ==  7 );    
#else  // the "normal" way
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
#endif
      
      int md = lhsmbr->getMajorDim();
      assert(  md == 6 ); 
      assert( lhsmbr->getNumElements() == 14 );
    }
    
  }

  // Test add/delete columns
  {    
    OsiGlpkSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    double inf = m.getInfinity();

    CoinPackedVector c0;
    c0.insert(0, 4);
    c0.insert(1, 1);
    m.addCol(c0, 0, inf, 3);
    m.initialSolve();
    double objValue = m.getObjValue();
    CoinRelFltEq eq(1.0e-2);
    assert( eq(objValue,2520.57) );
    // Try deleting first column that's nonbasic at lower bound (0).
    int * d = new int[1];
    CoinWarmStartBasis *cwsb =
	dynamic_cast<CoinWarmStartBasis *>(m.getWarmStart()) ;
    assert(cwsb) ;
    CoinWarmStartBasis::Status stati ;
    int iCol ;
    for (iCol = 0 ;  iCol < cwsb->getNumStructural() ; iCol++)
    { stati = cwsb->getStructStatus(iCol) ;
      if (stati == CoinWarmStartBasis::atLowerBound) break ; }
    d[0]=iCol;
    m.deleteCols(1,d);
    delete [] d;
    delete cwsb;
    d=NULL;
    m.resolve();
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );
    // Try deleting column we added. If basic, go to initialSolve as deleting
    // basic variable trashes basis required for warm start.
    iCol = m.getNumCols()-1;
    cwsb = dynamic_cast<CoinWarmStartBasis *>(m.getWarmStart()) ;
    stati =  cwsb->getStructStatus(iCol) ;
    delete cwsb;
    m.deleteCols(1,&iCol);
    if (stati == CoinWarmStartBasis::basic)
    { m.initialSolve() ; }
    else
    { m.resolve(); }
    objValue = m.getObjValue();
    assert( eq(objValue,2520.57) );

  }

#if 0
  // ??? Simplex routines not adapted to OsiGlpk yet

  // Solve an lp by hand
  {    
    OsiGlpkSolverInterface m;
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
      int returnCode=m.primalPivotResult(colIn,direction,colOut,outStatus,theta,NULL);
      assert (!returnCode);
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
  // Solve an lp when interface is on
  {    
    OsiGlpkSolverInterface m;
    std::string fn = mpsDir+"p0033";
    m.readMps(fn.c_str(),"mps");
    // enable special mode
    m.setHintParam(OsiDoScale,false,OsiHintDo);
    m.setHintParam(OsiDoPresolveInInitial,false,OsiHintDo);
    m.setHintParam(OsiDoDualInInitial,false,OsiHintDo);
    m.setHintParam(OsiDoPresolveInResolve,false,OsiHintDo);
    m.setHintParam(OsiDoDualInResolve,false,OsiHintDo);
    m.enableSimplexInterface(true);
    m.initialSolve();
  }
  // Check tableau stuff when simplex interface is on
  {    
    OsiGlpkSolverInterface m;
    /* 
       Wolsey : Page 130
       max 4x1 -  x2
       7x1 - 2x2    <= 14
       x2    <= 3
       2x1 - 2x2    <= 3
       x1 in Z+, x2 >= 0
    */
    
    double inf_ = m.getInfinity();
    int n_cols = 2;
    int n_rows = 3;
    
    double obj[2] = {-4.0, 1.0};
    double collb[2] = {0.0, 0.0};
    double colub[2] = {inf_, inf_};
    double rowlb[3] = {-inf_, -inf_, -inf_};
    double rowub[3] = {14.0, 3.0, 3.0};
    
    int rowIndices[5] =  {0,     2,    0,    1,    2};
    int colIndices[5] =  {0,     0,    1,    1,    1};
    double elements[5] = {7.0, 2.0, -2.0,  1.0, -2.0};
    CoinPackedMatrix M(true, rowIndices, colIndices, elements, 5);
    
    m.loadProblem(M, collb, colub, obj, rowlb, rowub);
    m.enableSimplexInterface(true);
    
    m.initialSolve();
    
    //check that the tableau matches wolsey (B-1 A)
    // slacks in second part of binvA
    double * binvA = (double*) malloc((n_cols+n_rows) * sizeof(double));
    
    printf("B-1 A");
    for(int i = 0; i < n_rows; i++){
      m.getBInvARow(i, binvA,binvA+n_cols);
      printf("\nrow: %d -> ",i);
      for(int j=0; j < n_cols+n_rows; j++){
	printf("%g, ", binvA[j]);
      }
    }
    printf("\n");
    m.disableSimplexInterface();
    free(binvA);
  }
#endif 

/*
  Read in exmip1 and solve it with verbose output setting.
*/
  { OsiGlpkSolverInterface osi ;
    std::cout << "Boosting verbosity.\n" ;
    osi.setHintParam(OsiDoReducePrint,false,OsiForceDo) ;

    std::string exmpsfile = mpsDir+"exmip1" ;
    std::string probname ;
    std::cout << "Reading mps file \"" << exmpsfile << "\"\n" ;
    osi.readMps(exmpsfile.c_str(), "mps") ;
    assert(osi.getStrParam(OsiProbName,probname)) ;
    std::cout << "Solving " << probname << " ... \n" ;
    osi.initialSolve() ;
    double val = osi.getObjValue() ;
    std::cout << "And the answer is " << val << ".\n" ;
    assert(fabs(val - 3.23) < 0.01) ;
  }

  return errCnt;
}
#endif
