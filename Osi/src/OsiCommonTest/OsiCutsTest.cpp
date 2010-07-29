// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CoinPragma.hpp"

#include "OsiUnitTests.hpp"

#include "OsiCuts.hpp"

//--------------------------------------------------------------------------
void
OsiCutsUnitTest()
{
  CoinRelFltEq eq;
  // Test default constructor
  {
    OsiCuts r;
    assert( r.colCutPtrs_.empty() );
    assert( r.rowCutPtrs_.empty() );
    assert( r.sizeColCuts()==0 );
    assert( r.sizeRowCuts()==0 );
    assert( r.sizeCuts()==0 );
    assert( r.mostEffectiveCutPtr()==NULL );
  }
  
  // Create some cuts
  OsiRowCut rcv[5];
  OsiColCut ccv[5];
  OsiCuts cuts;
  int i;
  for ( i=0; i<5; i++ ) {
    rcv[i].setEffectiveness(100.+i);
    ccv[i].setEffectiveness(200.+i);
    cuts.insert(rcv[i]);
    cuts.insert(ccv[i]);
  }
  
  OsiCuts rhs;
  {
    OsiCuts cs;
    
    // test inserting & accessing cut
    for ( i=0; i<5; i++ ) {
      assert( cs.sizeRowCuts()==i);
      assert( cs.sizeCuts()==2*i);
      cs.insert(rcv[i]);
      assert( cs.sizeRowCuts()==i+1);      
      assert( cs.sizeCuts()==2*i+1);
      assert( cs.rowCut(i)==rcv[i] );
#if 0
      const OsiCut * cp = cs.cutPtr(2*i);
      assert( cs.rowCut(i).effectiveness() == cp->effectiveness() );
      const OsiRowCut * rcP = dynamic_cast<const OsiRowCut*>( cp );
      assert( rcP != NULL );
      assert( *rcP == rcv[i] );
      const OsiColCut * ccP = dynamic_cast<const OsiColCut*>( cs.cutPtr(2*i) );
      assert( ccP == NULL );
#endif
      
      assert( cs.sizeColCuts()==i);
      assert( cs.sizeCuts()==2*i+1);
      cs.insert(ccv[i]);
      assert( cs.sizeColCuts()==i+1);
      assert( cs.sizeCuts()==2*i+2);
      assert( cs.colCut(i)==ccv[i] );
#if 0
      rcP = dynamic_cast<const OsiRowCut*>( cs.cutPtr(2*i+1) );
      assert( rcP == NULL );
      ccP = dynamic_cast<const OsiColCut*>( cs.cutPtr(2*i+1) );
      assert( ccP != NULL );
      assert( *ccP == ccv[i] );
#endif
      assert( eq(cs.mostEffectiveCutPtr()->effectiveness(),200.0+i) );

    }  
    
    // Test inserting collection of OsiCuts
    {
      OsiCuts cs;
      cs.insert(cuts);
      assert(cs.sizeColCuts()==5);
      assert(cs.sizeRowCuts()==5);
      assert(cs.sizeCuts()==10);
    }
/*
  Test handling of cut pointers. Create a vector of cut pointers, then add
  them to the collection. Dump the row cuts, then add again. The cut
  objects should be deleted when the OsiCuts object is deleted at the end of
  the block.  This test should be monitored with some flavour of runtime
  access checking to be sure that cuts are not freed twice or not freed at
  all.
*/
    {
      OsiCuts cs;
      OsiRowCut *rcPVec[5] ;
      OsiColCut *ccPVec[5] ;
      for (i = 0 ; i < 5 ; i++) {
	rcPVec[i] = rcv[i].clone() ;
	ccPVec[i] = ccv[i].clone() ;
      }
      // test inserting cut ptr & accessing cut
      for ( i=0; i<5; i++ ) {
        assert( cs.sizeRowCuts()==i);
        OsiRowCut * rcP = rcPVec[i] ;
        assert( rcP != NULL );
        cs.insert(rcP);
        assert( rcP == NULL );
        assert( cs.sizeRowCuts()==i+1);
        assert( cs.rowCut(i)==rcv[i] );
        
        
        OsiColCut * ccP = ccPVec[i] ;
        assert( ccP != NULL );
        assert( cs.sizeColCuts()==i);
        cs.insert(ccP);
        assert( ccP == NULL );
        assert( cs.sizeColCuts()==i+1);
        assert( cs.colCut(i)==ccv[i] );

        assert( eq(cs.mostEffectiveCutPtr()->effectiveness(),200.0+i) );
      }
      cs.dumpCuts() ;		// row cuts only
      for ( i=0; i<5; i++ ) {
        assert( cs.sizeRowCuts()==i);
        OsiRowCut * rcP = rcPVec[i] ;
        assert( rcP != NULL );
        cs.insert(rcP);
        assert( rcP == NULL );
        assert( cs.sizeRowCuts()==i+1);
        assert( cs.rowCut(i)==rcv[i] );
        
        assert( cs.sizeColCuts()==5);
      }
    }
    
    // Test copy constructor
    OsiCuts csC(cs);
    assert( csC.sizeRowCuts()==5 );
    assert( csC.sizeColCuts()==5 );
    for ( i=0; i<5; i++ ) {
      assert( csC.rowCut(i) == rcv[i] );
      assert( csC.colCut(i) == ccv[i] );
      assert( csC.rowCut(i) == cs.rowCut(i) );
      assert( csC.colCut(i) == cs.colCut(i) );
    }
    assert( eq(csC.mostEffectiveCutPtr()->effectiveness(),204.0) );
    
    rhs=cs;
  }
  
  // Test results of assigmnet operation
  for ( i=0; i<5; i++ ) {
    assert( rhs.rowCut(i) == rcv[i] );
    assert( rhs.colCut(i) == ccv[i] );
  }
  assert( eq(rhs.mostEffectiveCutPtr()->effectiveness(),204.0) );

  // Test removing cuts
  {
    OsiCuts t(rhs);  
    assert( t.sizeRowCuts()==5 );
    assert( t.sizeColCuts()==5 );
    assert( eq(rhs.mostEffectiveCutPtr()->effectiveness(),204.0) );
    assert( eq(t.mostEffectiveCutPtr()->effectiveness(),204.0) );

    t.eraseRowCut(3);     
    assert( t.sizeRowCuts()==4 );
    assert( t.sizeColCuts()==5 );
    for ( i=0; i<5; i++ ) {
      assert( t.colCut(i) == ccv[i] );
    }
    assert( t.rowCut(0)==rcv[0] );
    assert( t.rowCut(1)==rcv[1] );
    assert( t.rowCut(2)==rcv[2] );
    assert( t.rowCut(3)==rcv[4] );
    assert( eq(t.mostEffectiveCutPtr()->effectiveness(),204.0) );

    t.eraseColCut(4);     
    assert( t.sizeRowCuts()==4 );
    assert( t.sizeColCuts()==4 );
    for ( i=0; i<4; i++ ) {
      assert( t.colCut(i) == ccv[i] );
    }
    assert( t.rowCut(0)==rcv[0] );
    assert( t.rowCut(1)==rcv[1] );
    assert( t.rowCut(2)==rcv[2] );
    assert( t.rowCut(3)==rcv[4] );
    assert( eq(t.mostEffectiveCutPtr()->effectiveness(),203.0) );
  }
  

  
  // sorting cuts
  {
    OsiCuts t(rhs);
    assert( t.sizeRowCuts()==5 );
    assert( t.sizeColCuts()==5 );
    t.rowCut(0).setEffectiveness(9.);
    t.rowCut(1).setEffectiveness(1.);
    t.rowCut(2).setEffectiveness(7.);
    t.rowCut(3).setEffectiveness(3.);
    t.rowCut(4).setEffectiveness(5.);
    t.colCut(0).setEffectiveness(2.);
    t.colCut(1).setEffectiveness(8.);
    t.colCut(2).setEffectiveness(4.);
    t.colCut(3).setEffectiveness(6.);
    t.colCut(4).setEffectiveness(.5);
    double totEff=1.+2.+3.+4.+5.+6.+7.+8.+9.+0.5;

    {
      // Test iterator over all cuts
      double sumEff=0.;
      for ( OsiCuts::iterator it=t.begin(); it!=t.end(); ++it ) {
        double eff=(*it)->effectiveness();
        sumEff+= eff;
      }
      assert( sumEff == totEff );
    }

    t.sort();
    for ( i=1; i<5; i++ ) assert( t.colCut(i-1)>t.colCut(i) );
    for ( i=1; i<5; i++ ) assert( t.rowCut(i-1)>t.rowCut(i) );

    {
      // Test iterator over all cuts
      double sumEff=0.;
      for ( OsiCuts::iterator it=t.begin(); it!=t.end(); ++it ) {
        sumEff+= (*it)->effectiveness();
      }
      assert( sumEff == totEff );
    }

    {
      OsiCuts::iterator it=t.begin();
      OsiCut * cm1 = *it;
      ++it;
      for( ; it!=t.end(); it++ ) {
        OsiCut * c = *it;
        assert( (*cm1)>(*c) );
        cm1 = c;
      }
    }
  }
}
