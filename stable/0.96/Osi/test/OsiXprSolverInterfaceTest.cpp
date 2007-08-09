// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "OsiConfig.h"

#include <cassert>

#ifndef __ANSIC_
#  define __ANSIC_
#  include <xpresso.h>
#  undef  __ANSIC_
#else
#  include <xpresso.h>
#endif

#include "OsiXprSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "CoinPackedMatrix.hpp"

//-----------------------------------------------------------------------
// Test XPRESS-MP solution methods.
void
OsiXprSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

#if 0
  // Test to at least see if licence managment is working
  {
    int iret = initlz(NULL, 0);
    if ( iret != 0 ) getipv(N_ERRNO, &iret);
    assert(iret == 0);
  }
#endif


  // Test default constructor
  {
    assert( OsiXprSolverInterface::getNumInstances()==0 );
    OsiXprSolverInterface m;
    //    assert( m.xprSaved_ == false );
    //    assert( m.xprMatrixId_ = -1 );
    assert( m.xprProbname_ == "" );
    assert( m.matrixByRow_ == NULL );
    assert( m.colupper_ == NULL );
    assert( m.collower_ == NULL );
    assert( m.rowupper_ == NULL );
    assert( m.rowlower_ == NULL );
    assert( m.rowsense_ == NULL );
    assert( m.rhs_      == NULL );
    assert( m.rowrange_ == NULL );
    assert( m.colsol_   == NULL );
    assert( m.rowprice_   == NULL );
    assert( m.ivarind_  == NULL );
    assert( m.ivartype_ == NULL );
    assert( m.vartype_  == NULL );
    assert( OsiXprSolverInterface::getNumInstances() == 1 );
    //    assert( OsiXprSolverInterface::xprCurrentProblem_ == NULL );
    assert( m.getApplicationData() == NULL );
    int i = 2346;
    m.setApplicationData(&i);
    assert( *((int *)(m.getApplicationData())) == i );
  }
  assert( OsiXprSolverInterface::getNumInstances() == 0 );

  {
    CoinRelFltEq eq;
    OsiXprSolverInterface m;
    assert( OsiXprSolverInterface::getNumInstances() == 1 );
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str());
    //    assert( OsiXprSolverInterface::xprCurrentProblem_ == &m );
    // This assert fails on windows because fn is mixed case and xprProbname_is uppercase.
    //assert( m.xprProbname_ == fn );
    int ad = 13579;
    m.setApplicationData(&ad);
    assert( *((int *)(m.getApplicationData())) == ad );

    {
      OsiXprSolverInterface im;    
      //      assert( im.modelPtr_==NULL );
      
      assert( im.getNumCols() == 0 ); 
      
      //      assert( im.modelPtr()!=NULL );
      //      assert( im.mutableModelPtr()!=NULL );
      //      assert( im.modelPtr() == im.mutableModelPtr() );
    }

    // Test copy constructor and assignment operator
    {
      OsiXprSolverInterface lhs;
      {      
        assert( *((int *)(m.getApplicationData())) == ad );
        OsiXprSolverInterface im(m);        
        assert( *((int *)(im.getApplicationData())) == ad );

        OsiXprSolverInterface imC1(im);
	//        assert( imC1.mutableModelPtr()!=im.mutableModelPtr() );
	//        assert( imC1.modelPtr()!=im.modelPtr() );
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        assert( *((int *)(imC1.getApplicationData())) == ad ); 
        
        //im.setModelPtr(m);
        
        
        OsiXprSolverInterface imC2(im);
	//        assert( imC2.mutableModelPtr()!=im.mutableModelPtr() );
	//        assert( imC2.modelPtr()!=im.modelPtr() );
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        assert( *((int *)(imC2.getApplicationData())) == ad ); 
        
	//        assert( imC2.mutableModelPtr()!=imC1.mutableModelPtr() );
	//        assert( imC2.modelPtr()!=imC1.modelPtr() );
        
        lhs=imC2;
      }

      // Test that lhs has correct values even though rhs has gone out of scope
      //      assert( lhs.mutableModelPtr() != m.mutableModelPtr() );
      //      assert( lhs.modelPtr() != m.modelPtr() );
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
      assert( *((int *)(lhs.getApplicationData())) == ad );
    }

    // Test clone
    {
      OsiXprSolverInterface xprSi(m);
      OsiSolverInterface * siPtr = &xprSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiXprSolverInterface * xprClone = dynamic_cast<OsiXprSolverInterface*>(siClone);
      assert( xprClone != NULL );
      //      assert( xprClone->modelPtr() != xprSi.modelPtr() );
      //      assert( xprClone->modelPtr() != m.modelPtr() );
      assert( xprClone->getNumRows() == xprSi.getNumRows() );
      assert( xprClone->getNumCols() == m.getNumCols() );
      
      assert( *((int *)(xprClone->getApplicationData())) == ad );
      delete siClone;
    }

    // Test infinity
    {
      OsiXprSolverInterface si;
      assert( eq(si.getInfinity(), DPLINF) );
    }

    // Test setting solution
    {
      OsiXprSolverInterface m1(m);
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
      OsiXprSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str());
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

      // XPRESS-MP incorrectly treats unbounded integer variables as
      // general integers instead of binary (as in MPSX standard)
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
	double * cs = new double[fim.getNumCols()];
	for ( int i = 0;  i < fim.getNumCols();  cs[i++] = 0.0 );
	cs[2] = 2.9;
	cs[3] = 3.0;
	fim.setColSolution(cs);

        OsiVectorInt fi = fim.getFractionalIndices();
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        cs[2] = 5 + .00001/2.;
        cs[3] = 8 - .00001/2.;
	fim.setColSolution(cs);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        cs[2] = 5 + .00001*2.;
        cs[3] = 8 - .00001*2.;
	fim.setColSolution(cs);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );

	delete [] cs;
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
    
    // Test apply cut method
    {      
      OsiXprSolverInterface im(m);
      OsiCuts cuts;
      
      // Generate some cuts 
      //cg.generateCuts(im,cuts);
      {
        // Get number of rows and columns in model
        int nr=im.getNumRows();
        int nc=im.getNumCols();
	assert ( nr == 5 );
	assert ( nc == 8 );
        
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
          const double * xprColLB = im.getColLower();
          const double * xprColUB = im.getColUpper();
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *lb = new double[nc];
          double *ub = new double[nc];
          for (c=0;c<nc;c++) lb[c]=xprColLB[c]+0.001;
          for (c=0;c<nc;c++) ub[c]=xprColUB[c]-0.001;
          
          OsiColCut cc;
          cc.setLbs(nc,inx,lb);
          cc.setUbs(nc,inx,ub);
          
          cuts.insert(cc);
          delete [] ub;
          delete [] lb;
          delete [] inx;
        }
        
        {
          // Generate a row and column cut which have are ineffective
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
      OsiXprSolverInterface xprSi(m);
      int nc = xprSi.getNumCols();
      int nr = xprSi.getNumRows();
      assert( nc == 8 );
      assert( nr == 5 );
      assert( eq(xprSi.getColLower()[0],2.5) );
      assert( eq(xprSi.getColLower()[1],0.0) );
      assert( eq(xprSi.getColUpper()[1],4.1) );
      assert( eq(xprSi.getRowLower()[0],2.5) );
      assert( eq(xprSi.getRowLower()[4],3.0) );
      assert( eq(xprSi.getRowUpper()[1],2.1) );
      assert( eq(xprSi.getRowUpper()[4],15.0) );
      
      //      const double * cs = xprSi.getColSolution();
      //      assert( eq(cs[0],2.5) );
      //      assert( eq(cs[7],0.0) );
      
      assert( !eq(xprSi.getColLower()[3],1.2345) );
      xprSi.setColLower( 3, 1.2345 );
      assert( eq(xprSi.getColLower()[3],1.2345) );
      
      assert( !eq(xprSi.getColUpper()[4],10.2345) );
      xprSi.setColUpper( 4, 10.2345 );
      assert( eq(xprSi.getColUpper()[4],10.2345) );

      //assert( eq(xprSi.getObjValue(),0.0) );

      assert( eq( xprSi.getObjCoefficients()[0],  1.0) );
      assert( eq( xprSi.getObjCoefficients()[1],  0.0) );
      assert( eq( xprSi.getObjCoefficients()[2],  0.0) );
      assert( eq( xprSi.getObjCoefficients()[3],  0.0) );
      assert( eq( xprSi.getObjCoefficients()[4],  2.0) );
      assert( eq( xprSi.getObjCoefficients()[5],  0.0) );
      assert( eq( xprSi.getObjCoefficients()[6],  0.0) );
      assert( eq( xprSi.getObjCoefficients()[7], -1.0) );
    }
    
    // Test matrixByRow method
    { 
      const OsiXprSolverInterface si(m);
      const CoinPackedMatrix * smP = si.getMatrixByRow();
      
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
      assert( smP->getMinorDim() == 8 ); 
      assert( smP->getNumElements() == 14 ); 
      assert( smP->getSizeVectorStarts()==6 );

    }
    
    
    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow
    {
      OsiXprSolverInterface lhs;
      {      
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rhs_==NULL );
        
        OsiXprSolverInterface siC1(m);     
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );

        const char   * siC1rs  = siC1.getRowSense();
        assert( siC1rs[0] == 'G' );
        assert( siC1rs[1] == 'L' );
        assert( siC1rs[2] == 'E' );
        assert( siC1rs[3] == 'R' );
	assert( siC1rs[4] == 'R' );

        const double * siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0], 2.5) );
        assert( eq(siC1rhs[1], 2.1) );
        assert( eq(siC1rhs[2], 4.0) );
        assert( eq(siC1rhs[3], 5.0) ); 
	assert( eq(siC1rhs[4], 15.0) );
	
        const double * siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0], 0.0) );
        assert( eq(siC1rr[1], 0.0) );
        assert( eq(siC1rr[2], 0.0) );
        assert( eq(siC1rr[3], 5.0 - 1.8) );
        assert( eq(siC1rr[4], 15.0 - 3.0) );
	
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
	assert( mi[2]==7  );
	assert( mi[3]==9  );
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
	assert( ei[7 ] ==  2 );
	assert( ei[8 ] ==  5 );
	assert( ei[9 ] ==  3 );
	assert( ei[10] ==  6 );  
	assert( ei[11] ==  0 );  
	assert( ei[12] ==  4 );  
	assert( ei[13] ==  7 );  
    
	assert( siC1mbr->getMajorDim() == 5 ); 
	assert( siC1mbr->getMinorDim() == 8 ); 
	assert( siC1mbr->getNumElements() == 14 ); 
	assert( siC1mbr->getSizeVectorStarts()==6 );

        assert( siC1rs  == siC1.getRowSense() );
        assert( siC1rhs == siC1.getRightHandSide() );
        assert( siC1rr  == siC1.getRowRange() );

        // Change XPRESS Model by adding free row
        OsiRowCut rc;
        rc.setLb(-DBL_MAX);
        rc.setUb(DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);
             
        // Since model was changed, test that cached
        // data is now freed.
        assert( siC1.rowrange_     == NULL );
        assert( siC1.rowsense_     == NULL );
        assert( siC1.rhs_          == NULL );
	assert( siC1.matrixByRow_ == NULL );
        
        siC1rs  = siC1.getRowSense();
        assert( siC1rs[0] == 'G' );
        assert( siC1rs[1] == 'L' );
        assert( siC1rs[2] == 'E' );
        assert( siC1rs[3] == 'R' );
        assert( siC1rs[4] == 'R' );
        assert( siC1rs[5] == 'N' );
	
        siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.0) ); 
        assert( eq(siC1rhs[5],0.0) ); 
	
        siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0], 0.0) );
        assert( eq(siC1rr[1], 0.0) );
        assert( eq(siC1rr[2], 0.0) );
        assert( eq(siC1rr[3], 5.0 - 1.8) );
        assert( eq(siC1rr[4], 15.0 - 3.0) );
        assert( eq(siC1rr[5], 0.0) );
        lhs=siC1;
      }
      // Test that lhs has correct values even though siC1 has gone out of scope    
      assert( lhs.rowrange_    == NULL );
      assert( lhs.rowsense_    == NULL );
      assert( lhs.rhs_         == NULL ); 
      assert( lhs.matrixByRow_ == NULL );
      
      const char * lhsrs  = lhs.getRowSense();
      assert( lhsrs[0] == 'G' );
      assert( lhsrs[1] == 'L' );
      assert( lhsrs[2] == 'E' );
      assert( lhsrs[3] == 'R' );
      assert( lhsrs[4] == 'R' );
      assert( lhsrs[5] == 'N' );

      
      const double * lhsrhs = lhs.getRightHandSide();
      assert( eq(lhsrhs[0], 2.5) );
      assert( eq(lhsrhs[1], 2.1) );
      assert( eq(lhsrhs[2], 4.0) );
      assert( eq(lhsrhs[3], 5.0) );
      assert( eq(lhsrhs[4], 15.0) ); 
      assert( eq(lhsrhs[5], 0.0) ); 
      
      const double *lhsrr  = lhs.getRowRange();
      assert( eq(lhsrr[0], 0.0) );
      assert( eq(lhsrr[1], 0.0) );
      assert( eq(lhsrr[2], 0.0) );
      assert( eq(lhsrr[3], 5.0 - 1.8) );
      assert( eq(lhsrr[4], 15.0 - 3.0) );
      assert( eq(lhsrr[5], 0.0) );
      
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
      assert( mi[3]==9  );
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
      
      assert( lhsmbr->getMajorDim() == 6 ); 
      assert( lhsmbr->getMinorDim() == 8 ); 
      assert( lhsmbr->getNumElements() == 14 ); 
      assert( lhsmbr->getSizeVectorStarts()==7 );
    }

    //--------------
    
    // Test load problem
    {
      OsiXprSolverInterface base(m);
      base.initialSolve();
      assert(m.getNumRows() == base.getNumRows());

      OsiXprSolverInterface si1,si2,si3,si4;

      si1.loadProblem(
        *base.getMatrixByCol(),
        base.getColLower(),base.getColUpper(),base.getObjCoefficients(),
        base.getRowSense(),base.getRightHandSide(),base.getRowRange());
      si1.initialSolve();
      assert(eq(base.getObjValue(), si1.getObjValue()));
      assert(m.getNumRows() == si1.getNumRows());

      si2.loadProblem(
        *base.getMatrixByRow(),
        base.getColLower(),base.getColUpper(),base.getObjCoefficients(),
        base.getRowSense(),base.getRightHandSide(),base.getRowRange());
      si2.initialSolve();
      assert(eq(base.getObjValue(), si2.getObjValue()));
      assert(m.getNumRows() == si2.getNumRows());

      si3.loadProblem(
        *base.getMatrixByCol(),
        base.getColLower(),base.getColUpper(),base.getObjCoefficients(),
        base.getRowLower(),base.getRowUpper() );
      si3.initialSolve();
      assert(eq(base.getObjValue(), si3.getObjValue()));
      assert(m.getNumRows() == si3.getNumRows());

      si4.loadProblem(
        *base.getMatrixByCol(),
        base.getColLower(),base.getColUpper(),base.getObjCoefficients(),
        base.getRowLower(),base.getRowUpper() );
      si4.initialSolve();
      assert(eq(base.getObjValue(), si4.getObjValue()));
      assert(m.getNumRows() == si4.getNumRows());

      base.initialSolve();
      si1.initialSolve();
      si2.initialSolve();
      si3.initialSolve();
      si4.initialSolve();

      // Create an indices vector
      assert(base.getNumCols()<10);
      assert(base.getNumRows()<10);
      int indices[10];
      int i;
      for (i=0; i<10; i++) indices[i]=i;

      // Test collower
      CoinPackedVector basePv,pv;
      basePv.setVector(base.getNumCols(),indices,base.getColLower());
      pv.setVector( si1.getNumCols(),indices, si1.getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2.getNumCols(),indices, si2.getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3.getNumCols(),indices, si3.getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4.getNumCols(),indices, si4.getColLower());
      assert(basePv.isEquivalent(pv));
      
      // Test colupper
      basePv.setVector(base.getNumCols(),indices,base.getColUpper());
      pv.setVector( si1.getNumCols(),indices, si1.getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2.getNumCols(),indices, si2.getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3.getNumCols(),indices, si3.getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4.getNumCols(),indices, si4.getColUpper());
      assert(basePv.isEquivalent(pv));
     
      // Test getObjCoefficients
      basePv.setVector(base.getNumCols(),indices,base.getObjCoefficients());
      pv.setVector( si1.getNumCols(),indices, si1.getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2.getNumCols(),indices, si2.getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3.getNumCols(),indices, si3.getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4.getNumCols(),indices, si4.getObjCoefficients());
      assert(basePv.isEquivalent(pv));

      // Test rowlower
      basePv.setVector(base.getNumRows(),indices,base.getRowLower());
      pv.setVector( si1.getNumRows(),indices, si1.getRowLower());
      assert( eq(base.getRowLower()[3],si1.getRowLower()[3]) );
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2.getNumRows(),indices, si2.getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3.getNumRows(),indices, si3.getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4.getNumRows(),indices, si4.getRowLower());
      assert(basePv.isEquivalent(pv));
   
      // Test rowupper
      basePv.setVector(base.getNumRows(),indices,base.getRowUpper());
      pv.setVector( si1.getNumRows(),indices, si1.getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2.getNumRows(),indices, si2.getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3.getNumRows(),indices, si3.getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4.getNumRows(),indices, si4.getRowUpper());
      assert(basePv.isEquivalent(pv));

      // Test Constraint Matrix
      assert( base.getMatrixByCol()->isEquivalent(*si1.getMatrixByCol()) );
      assert( base.getMatrixByRow()->isEquivalent(*si1.getMatrixByRow()) );
      assert( base.getMatrixByCol()->isEquivalent(*si2.getMatrixByCol()) );
      assert( base.getMatrixByRow()->isEquivalent(*si2.getMatrixByRow()) );
      assert( base.getMatrixByCol()->isEquivalent(*si3.getMatrixByCol()) );
      assert( base.getMatrixByRow()->isEquivalent(*si3.getMatrixByRow()) );
      assert( base.getMatrixByCol()->isEquivalent(*si4.getMatrixByCol()) );
      assert( base.getMatrixByRow()->isEquivalent(*si4.getMatrixByRow()) );

      // Test Objective Value
      assert( eq(base.getObjValue(),si1.getObjValue()) );
      assert( eq(base.getObjValue(),si2.getObjValue()) );
      assert( eq(base.getObjValue(),si3.getObjValue()) );
      assert( eq(base.getObjValue(),si4.getObjValue()) );
    }
    //--------------

    assert(OsiXprSolverInterface::getNumInstances()==1);
  }
  assert(OsiXprSolverInterface::getNumInstances()==0);

  
  // Do common solverInterface testing by calling the
  // base class testing method.
  {
    OsiXprSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir,netlibDir);
  }
}
