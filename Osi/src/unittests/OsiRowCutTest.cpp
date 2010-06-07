// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
\
#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <cfloat>

#include "CoinPragma.hpp"

#include "OsiUnitTests.hpp"

#include "OsiRowCut.hpp"
#include "CoinFloatEqual.hpp"

//--------------------------------------------------------------------------
void
OsiRowCutUnitTest(const OsiSolverInterface * baseSiP,
		  const std::string & mpsDir )
{

  CoinRelFltEq eq;

  // Test default constructor
  {
    OsiRowCut r;
    assert( r.row_.getIndices()==NULL );
    assert( r.row_.getElements()==NULL );
    assert( r.lb_==-/*std::numeric_limits<double>::max()*/DBL_MAX );
    assert( r.ub_== /*std::numeric_limits<double>::max()*/DBL_MAX );
  }

  // Test set and get methods
  const int ne = 4;
  int inx[ne] = { 1, 3, 4, 7 };
  double el[ne] = { 1.2, 3.4, 5.6, 7.8 };
  {
    OsiRowCut r;    
    assert( r.row().getNumElements()==0 );
    assert( r.effectiveness()==0. );
    //assert( r.timesUsed()==0 );
    //assert( r.timesTested()==0 );

    // Test setting getting bounds
    r.setLb(65.432);
    r.setUb(123.45);
    assert( r.lb()==65.432 );
    assert( r.ub()==123.45 );

    // Test setting/getting of effectiveness,timesUsed,timesTested
    r.setEffectiveness(45.);
    assert( r.effectiveness()==45. );
#if 0
    r.setTimesUsed(11);
    assert( r.timesUsed()==11 );
    r.incrementTimesUsed();
    assert( r.timesUsed()==12 );
    r.setTimesTested(111);
    assert( r.timesTested()==111 );
    r.incrementTimesTested();
    assert( r.timesTested()==112 );
#endif
    
    // Test setting/getting elements with int* & float* vectors
    r.setRow( ne, inx, el );
    assert( r.row().getNumElements()==ne );
    for ( int i=0; i<ne; i++ ) {
      assert( r.row().getIndices()[i] == inx[i] );
      assert( r.row().getElements()[i] == el[i] );
    }
  } 

  // Repeat test with ownership constructor
  {
    int *inxo = new int[ne];
    double *elo = new double[ne];
    double lb = 65.432;
    double ub = 123.45;

    int i;
    for ( i = 0 ; i < ne ; i++) inxo[i] = inx[i] ;
    for ( i = 0 ; i < ne ; i++) elo[i] = el[i] ;
    OsiRowCut r(lb,ub,ne,ne,inxo,elo);    

    assert( r.row().getNumElements()==ne );
    assert( r.effectiveness()==0. );
    //assert( r.timesUsed()==0 );
    //assert( r.timesTested()==0 );

    // Test getting bounds
    assert( r.lb()==lb );
    assert( r.ub()==ub );

    // Test setting/getting of effectiveness,timesUsed,timesTested
    r.setEffectiveness(45.);
    assert( r.effectiveness()==45. );
#if 0
    r.setTimesUsed(11);
    assert( r.timesUsed()==11 );
    r.incrementTimesUsed();
    assert( r.timesUsed()==12 );
    r.setTimesTested(111);
    assert( r.timesTested()==111 );
    r.incrementTimesTested();
    assert( r.timesTested()==112 );
#endif
    
    // Test getting elements with int* & float* vectors
    assert( r.row().getNumElements()==ne );
    for ( i=0; i<ne; i++ ) {
      assert( r.row().getIndices()[i] == inx[i] );
      assert( r.row().getElements()[i] == el[i] );
    }
  } 

  // Test sense, rhs, range
  {
    {
      OsiRowCut r;
      assert( r.sense() == 'N' );
      assert( r.rhs() == 0.0 );
      assert( r.range() == 0.0 );
    }
    {
      OsiRowCut r;
      r.setLb(65.432);
      assert( r.sense() == 'G' );
      assert( r.rhs() == 65.432 );
      assert( r.range() == 0.0 );
    }
    {
      OsiRowCut r;
      r.setLb(65.432);
      r.setUb(65.432);
      assert( r.sense() == 'E' );
      assert( r.rhs() == 65.432 );
      assert( r.range() == 0.0 );
    }
    {
      OsiRowCut r;
      r.setUb(123.45);
      assert( r.sense() == 'L' );
      assert( r.rhs() == 123.45 );
      assert( r.range() == 0.0 );
    }
    {
      OsiRowCut r;
      r.setLb(65.432);
      r.setUb(123.45);
      assert( r.sense() == 'R' );
      assert( r.rhs() == 123.45 );
      assert( eq(r.range(),123.45 - 65.432) );
    }
  }
  
  // Test copy constructor and assignment operator
  {
    OsiRowCut rhs;
    {
      OsiRowCut r;
      OsiRowCut rC1(r);
      assert( rC1.row().getNumElements()==r.row().getNumElements() );
      assert( rC1.row().getIndices()==r.row().getIndices() );
      assert( rC1.row().getElements()==r.row().getElements() );
      assert( rC1.lb()==r.lb() );
      assert( rC1.ub()==r.ub() );
      
      r.setLb(65.432);
      r.setUb(123.45);
      r.setRow( ne, inx, el );
      r.setEffectiveness(123.);
      
      assert( rC1.row().getNumElements()!=r.row().getNumElements() );
      assert( rC1.row().getIndices()!=r.row().getIndices() );
      assert( rC1.row().getElements()!=r.row().getElements() );
      assert( rC1.lb()!=r.lb() );
      assert( rC1.ub()!=r.ub() );
      assert( rC1.effectiveness()!=r.effectiveness() );
      
      OsiRowCut rC2(r);
      assert( rC2.row().getNumElements()==r.row().getNumElements() );
      for( int i=0; i<r.row().getNumElements(); i++ ){
        assert( rC2.row().getIndices()[i]==r.row().getIndices()[i] );
        assert( rC2.row().getElements()[i]==r.row().getElements()[i] );
      }
      assert( rC2.lb()==r.lb() );
      assert( rC2.ub()==r.ub() );       
      assert( rC2.effectiveness()==r.effectiveness() );
      
      rhs=rC2;
    }
    // Test that rhs has correct values even though lhs has gone out of scope
    assert( rhs.row().getNumElements()==ne );
    for ( int i=0; i<ne; i++ ) {
      assert( rhs.row().getIndices()[i] == inx[i] );
      assert( rhs.row().getElements()[i] == el[i] );
    }  
    assert( rhs.effectiveness()==123. );
    assert( rhs.lb()==65.432 ); 
    assert( rhs.ub()==123.45 ); 
  }

  // Test setting row with packed vector
  {
    CoinPackedVector r;
    r.setVector(ne,inx,el);

    OsiRowCut rc;
    assert( rc.row()!=r );
    rc.setRow(r);
    assert( rc.row()==r );
  }

  // Test operator==
  {
    CoinPackedVector r;
    r.setVector(ne,inx,el);

    OsiRowCut rc;
    rc.setRow(r);
    rc.setEffectiveness(2.);
    rc.setLb(3.0);
    rc.setUb(4.0);

    {
      OsiRowCut c(rc);
      assert(   c == rc  );
      assert( !(c != rc) );
      assert( !(c  < rc)  );
      assert( !(rc < c )  );
      assert( !(c  > rc)  );
      assert( !(rc > c )  );
    }
    {
      OsiRowCut c(rc);      
      const int ne1 = 4;
      int inx1[ne] = { 1, 3, 4, 7 };
      double el1[ne] = { 1.2, 3.4, 5.6, 7.9 };
      c.setRow(ne1,inx1,el1);
      assert( !(c == rc) );
      assert(   c != rc  );
      assert( !(rc < c)  );
    }
    {
      OsiRowCut c(rc); 
      c.setEffectiveness(3.);
      assert( !(c == rc) );
      assert(   c != rc  );
      assert(   rc < c   );
      assert( !(c < rc)  );
      assert( !(rc > c)  );
      assert(   c > rc  );
    }
    {
      OsiRowCut c(rc); 
      c.setLb(4.0);
      assert( !(c == rc) );
      assert(   c != rc  );
    }
    {
      OsiRowCut c(rc); 
      c.setUb(5.0);
      assert( !(c == rc) );
      assert(   c != rc  );
    }
  }
#ifndef COIN_NOTEST_DUPLICATE
  {
    // Test consistent
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");

    OsiRowCut c;      
    const int ne = 3;
    int inx[ne] = { -1, 5, 4 };
    double el[ne] = { 1., 1., 1. };
    c.setRow(ne,inx,el);
    assert( !c.consistent() );;
    
    inx[0]=5;
#if 0
    c.setRow(ne,inx,el);
    assert( !c.consistent() );
#else
    bool errorThrown = false;
    try {
       c.setRow(ne,inx,el);
    }
    catch (CoinError e) {
       errorThrown = true;
    }
    assert(errorThrown);
#endif
    
    inx[0]=3;
    c.setRow(ne,inx,el);
    assert( c.consistent() );
    
    c.setLb(5.);
    c.setUb(5.);
    assert( c.consistent() );

    c.setLb(5.5);
    assert( c.consistent() );
    assert( c.infeasible(*imP) );
    delete imP;
  }
#endif
  {
    // Test consistent(IntegerModel) method.
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    
    OsiRowCut cut;    
    const int ne = 3;
    int inx[ne] = { 3, 5, 4 };
    double el[ne] = { 1., 1., 1. };
    cut.setRow(ne,inx,el);
    assert( cut.consistent() );
    
    inx[0]=7;
    cut.setRow(ne,inx,el);
    assert( cut.consistent(*imP) );
    
    inx[0]=8;
    cut.setRow(ne,inx,el);
    assert(  cut.consistent() ); 
    assert( !cut.consistent(*imP) );
    delete imP;
  }
  {
    //Test infeasible(im) method
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    
    OsiRowCut cut;   
    const int ne = 3;
    int inx[ne] = { 3, 5, 4 };
    double el[ne] = { 1., 1., 1. };
    cut.setRow(ne,inx,el);
    assert( !cut.infeasible(*imP) );
    
    OsiRowCut c1;
    assert( !c1.infeasible(*imP) );
    assert( c1.consistent() );
    assert( c1.consistent(*imP) );
    delete imP;
  }
  {
    //Test violation

    double solution[]={1.0,1.0,1.0,2.0,2.0,2.0};
    OsiRowCut cut;   
    const int ne = 3;
    int inx[ne] = { 3, 5, 4 };
    double el[ne] = { 1., 1., 1. };
    cut.setRow(ne,inx,el);
    cut.setUb(5.);
    assert( cut.violated(solution) );

    cut.setLb(5.);
    cut.setUb(10.);
    assert( !cut.violated(solution) );

    cut.setLb(6.1);
    assert( cut.violated(solution) );


  }
}

