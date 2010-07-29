// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifdef NDEBUG
#undef NDEBUG
#endif

//#include <cstdlib>
#include <cassert>
//#include <cmath>

#include "CoinPragma.hpp"

#include "OsiUnitTests.hpp"

#include "OsiColCut.hpp"

//--------------------------------------------------------------------------
void
OsiColCutUnitTest(const OsiSolverInterface * baseSiP, 
		  const std::string & mpsDir)
{

  // Test default constructor
  {
    OsiColCut r;
    assert( r.lbs_.getIndices()==NULL );
    assert( r.lbs_.getElements()==NULL );
    assert( r.ubs_.getIndices()==NULL );
    assert( r.ubs_.getElements()==NULL );
    assert( r.lbs_.getNumElements()==0);
    assert( r.ubs_.getNumElements()==0);
  }

  // Test set and get methods
  const int ne = 4;
  int inx[ne] = { 1, 3, 4, 7 };
  double el[ne] = { 1.2, 3.4, 5.6, 7.8 };
  const int ne3 = 0;
  int * inx3=NULL;
  double * el3= NULL;
  {
    OsiColCut r;    
    assert( r.lbs().getNumElements()==0 );
    assert( r.ubs().getNumElements()==0 );
        
    // Test setting/getting bounds
    r.setLbs( ne, inx, el );
    r.setEffectiveness(222.);
    assert( r.lbs().getNumElements()==ne );
    for ( int i=0; i<ne; i++ ) {
      assert( r.lbs().getIndices()[i] == inx[i] );
      assert( r.lbs().getElements()[i] == el[i] );
    }
    assert( r.effectiveness()==222. );
    
    r.setUbs( ne3, inx3, el3 );
    assert( r.ubs().getNumElements()==0 );
    assert( r.ubs().getIndices()==NULL );
    assert( r.ubs().getElements()==NULL );
    
  } 
  
  // Test copy constructor and assignment operator
  {
    OsiColCut rhs;
    {
      OsiColCut r;
      OsiColCut rC1(r);
      assert( rC1.lbs().getNumElements()==r.lbs().getNumElements() );
      assert( rC1.ubs().getNumElements()==r.ubs().getNumElements() );
  
      r.setLbs( ne, inx, el );
      r.setUbs( ne, inx, el );
      r.setEffectiveness(121.);

      assert( rC1.lbs().getNumElements()!=r.lbs().getNumElements() );
      assert( rC1.ubs().getNumElements()!=r.ubs().getNumElements() );
                
      OsiColCut rC2(r);
      assert( rC2.lbs().getNumElements()==r.lbs().getNumElements() );
      assert( rC2.ubs().getNumElements()==r.ubs().getNumElements() );
      assert( rC2.lbs().getNumElements()==ne );
      assert( rC2.ubs().getNumElements()==ne );
      for ( int i=0; i<ne; i++ ) {
        assert( rC2.lbs().getIndices()[i] == inx[i] );
        assert( rC2.lbs().getElements()[i] == el[i] );
        assert( rC2.ubs().getIndices()[i] == inx[i] );
        assert( rC2.ubs().getElements()[i] == el[i] );
      }
      assert( rC2.effectiveness()==121. );
      
      rhs=rC2;
    }
    // Test that rhs has correct values even though lhs has gone out of scope
    assert( rhs.lbs().getNumElements()==ne );
    assert( rhs.ubs().getNumElements()==ne );
    for ( int i=0; i<ne; i++ ) {
      assert( rhs.lbs().getIndices()[i] == inx[i] );
      assert( rhs.lbs().getElements()[i] == el[i] );
      assert( rhs.ubs().getIndices()[i] == inx[i] );
      assert( rhs.ubs().getElements()[i] == el[i] );
    }  
    assert( rhs.effectiveness()==121. );
  }

  // Test setting bounds with packed vector and operator==
  {    
    const int ne1 = 4;
    int inx1[ne] = { 1, 3, 4, 7 };
    double el1[ne] = { 1.2, 3.4, 5.6, 7.8 };
    const int ne2 = 2;
    int inx2[ne2]={ 1, 3 };
    double el2[ne2]= { 1.2, 3.4 };
    CoinPackedVector v1,v2;
    v1.setVector(ne1,inx1,el1);
    v2.setVector(ne2,inx2,el2);
    
    OsiColCut c1,c2;
    assert(   c1==c2  );
    assert( !(c1!=c2) );
    
    c1.setLbs(v1);
    assert(   c1 != c2  );
    assert( !(c1 == c2) );
    
    assert( c1.lbs()==v1 );
    c1.setUbs(v2);
    assert( c1.ubs()==v2 );
    c1.setEffectiveness(3.);
    assert( c1.effectiveness()==3. );

    {
      OsiColCut c3(c1);
      assert(   c3==c1 );
      assert( !(c3!=c1) );
    }
    {
      OsiColCut c3(c1);
      c3.setLbs(v2);
      assert(   c3!=c1 );
      assert( !(c3==c1) );
    }
    {
      OsiColCut c3(c1);
      c3.setUbs(v1);
      assert(   c3!=c1 );
      assert( !(c3==c1) );
    }
    {
      OsiColCut c3(c1);
      c3.setEffectiveness(5.);
      assert(   c3!=c1 );
      assert( !(c3==c1) );
    }
  }

  // internal consistency 
  {
    const int ne = 1;
    int inx[ne] = { -3 };
    double el[ne] = { 1.2 };
    OsiColCut r; 
    r.setLbs( ne, inx, el );
    assert( !r.consistent() );
  }
  {
    const int ne = 1;
    int inx[ne] = { -3 };
    double el[ne] = { 1.2 };
    OsiColCut r; 
    r.setUbs( ne, inx, el );
    assert( !r.consistent() );
  }
  {
    const int ne = 1;
    int inx[ne] = { 100 };
    double el[ne] = { 1.2 };
    const int ne1 = 2;
    int inx1[ne1] = { 50, 100 };
    double el1[ne1] = { 100., 100. };
    OsiColCut r; 
    r.setUbs( ne, inx, el );
    r.setLbs( ne1, inx1, el1 );
    assert( r.consistent() );

    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    assert( !r.consistent(*imP) );
    delete imP;
  }
  {
    const int ne = 1;
    int inx[ne] = { 100 };
    double el[ne] = { 1.2 };
    const int ne1 = 2;
    int inx1[ne1] = { 50, 100 };
    double el1[ne1] = { 100., 1. };
    OsiColCut r; 
    r.setUbs( ne, inx, el );
    r.setLbs( ne1, inx1, el1 );
    assert( r.consistent() );
  }
  {
    // Test consistent(IntegerModel) method.
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    
    OsiColCut cut;
    const int ne=1;
    int inx[ne]={20};
    double el[ne]={0.25};
    cut.setLbs(ne,inx,el);
    assert( !cut.consistent(*imP) );
    
    cut.setLbs(0,NULL,NULL);
    cut.setUbs(ne,inx,el); 
    assert( !cut.consistent(*imP) );
    
    inx[0]=4;
    cut.setLbs(ne,inx,el);
    cut.setUbs(0,NULL,NULL); 
    assert( cut.consistent(*imP) );
    
    el[0]=4.5;
    cut.setLbs(0,NULL,NULL);
    cut.setUbs(ne,inx,el); 
    assert( cut.consistent(*imP) );

    cut.setLbs(ne,inx,el);
    cut.setUbs(0,NULL,NULL); 
    assert( cut.consistent(*imP) );  // Vailid but infeasible
    assert( cut.infeasible(*imP) );
    
    el[0]=3.0;
    cut.setLbs(ne,inx,el);
    cut.setUbs(ne,inx,el); 
    assert( cut.consistent(*imP) ); 
    delete imP;
  }
  {
    //Test infeasible(im) method
    // Test consistent(IntegerModel) method.
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    
    OsiColCut cut;
    const int ne=1;
    int inx[ne]={4};
    double el[ne]={4.5};
    cut.setLbs(ne,inx,el);
    assert( cut.infeasible(*imP) );

    el[0]=0.25;
    cut.setLbs(0,NULL,NULL);
    cut.setUbs(ne,inx,el);
    assert( cut.infeasible(*imP) ); 
    
    el[0]=3.0;
    cut.setLbs(ne,inx,el);
    cut.setUbs(ne,inx,el); 
    assert( !cut.infeasible(*imP) );

    delete imP;
  }
  {
    //Test violation

    double solution[]={1.0};
    OsiColCut cut;
    const int ne=1;
    int inx[ne]={0};
    double el[ne]={4.5};
    cut.setLbs(ne,inx,el);
    assert( cut.violated(solution) );

    el[0]=0.25;
    cut.setLbs(0,NULL,NULL);
    cut.setUbs(ne,inx,el);
    assert( cut.violated(solution) );
    
    el[0]=1.0;
    cut.setLbs(ne,inx,el);
    cut.setUbs(ne,inx,el); 
    assert( !cut.violated(solution) );

  }
}
