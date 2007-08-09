// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#ifdef NDEBUG
# undef NDEBUG
#endif

#include "OsiConfig.h"

#include <cassert>

#include "OsiFmpSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#include <vector>
#include <map>
#ifndef _MSC_VER
  using std::max;
  using std::min;
#endif


// Added so build windows build with dsp files works,
// when not building with cplex.
#ifdef COIN_HAS_FMP

//#############################################################################


void
OsiFmpSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir)
{

  // Test default constructor
  {
    OsiFmpSolverInterface m;
    assert( m.rowsense_==NULL );
    assert( m.rowrhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.matrixByColumn_==NULL );
  }
  
  
  {    
    CoinRelFltEq eq;
    OsiFmpSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    int ad = 13579;
    m.setApplicationData(&ad);
    assert( *((int *)(m.getApplicationData())) == ad );
    
    {
      OsiFmpSolverInterface im;    
      assert( im.getNumCols() == 0 ); 
    
    }
    
    // Test copy constructor and assignment operator
    {
      OsiFmpSolverInterface lhs;
      {      
        assert( *((int *)(m.getApplicationData())) == ad );
        OsiFmpSolverInterface im(m);        
        assert( *((int *)(im.getApplicationData())) == ad );

        OsiFmpSolverInterface imC1(im);

        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        assert( *((int *)(imC1.getApplicationData())) == ad ); 
        
        OsiFmpSolverInterface imC2(im);
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        assert( *((int *)(imC2.getApplicationData())) == ad ); 
        
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope  
      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
      assert( *((int *)(lhs.getApplicationData())) == ad );
    }
    
    // Test clone
    {
      OsiFmpSolverInterface fmpSi(m);
      OsiSolverInterface * siPtr = &fmpSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiFmpSolverInterface * fmpClone = dynamic_cast<OsiFmpSolverInterface*>(siClone);
      assert( fmpClone != NULL );
      assert( fmpClone->getNumRows() == fmpSi.getNumRows() );
      assert( fmpClone->getNumCols() == m.getNumCols() );
      
      assert( *((int *)(fmpClone->getApplicationData())) == ad );
      // Test reset
      siClone->reset();
      assert( fmpClone->rowsense_==NULL );
      assert( fmpClone->rowrhs_==NULL );
      assert( fmpClone->rowrange_==NULL );
      assert( fmpClone->matrixByRow_==NULL );
      assert( fmpClone->ws_==NULL);
      assert( fmpClone->itlimOrig_==9999999);
      delete siClone;
    }
   
    // test infinity
    {
      OsiFmpSolverInterface si;
      assert( eq(si.getInfinity(),FMP_INFINITY));
    }     

    //--------------
    // Test rowsense, rhs, rowrange, matrixByRow
    {
      OsiFmpSolverInterface lhs;
      {      
        assert( m.rowrange_==NULL );
        assert( m.rowsense_==NULL );
        assert( m.rowrhs_==NULL );
        assert( m.matrixByRow_==NULL );
        
        OsiFmpSolverInterface siC1(m);     
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rowrhs_==NULL );
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

        // Change the model by adding free row
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
        assert( siC1.rowrhs_==NULL );
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
      assert( lhs.rowrhs_==NULL ); 
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
    
  }

	// Test iterations count
	{
		int k;
		OsiFmpSolverInterface m;
    std::string fn = netlibDir+"25fv47";
    m.readMps(fn.c_str(),"mps");
    m.initialSolve();
		k=m.getIterationCount();
    assert(k!=0);
	}

	// TEST issue with 80bau3b
	{
		int k;
		OsiFmpSolverInterface m;
    std::string fn = netlibDir+"80bau3b";
    m.readMps(fn.c_str(),"mps");
		m.writeMps("80bau3bpat.mps");
    m.initialSolve();
		k=m.getIterationCount();
    assert(k!=0);
	}

  // Do common solverInterface testing 
  {
    OsiFmpSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir, netlibDir);
  }
}
#endif
