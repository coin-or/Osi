// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiIndexedVector.hpp"
#include "OsiShallowPackedVector.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

//--------------------------------------------------------------------------
void
OsiIndexedVectorUnitTest()
{
  
  int i;
  // Test default constructor
  {
    OsiIndexedVector r;
    assert( r.indices_==NULL );
    assert( r.packedElements_==NULL );
    assert( r.elements_==NULL );
    assert( r.getNumElements()==0 );
    assert( r.capacity_==0);
  }
  
  // Test set and get methods
  const int ne = 4;
  int inx[ne] = { 1, 3, 4, 7 };
  double el[ne] = { 1.2, 3.4, 5.6, 7.8 };
  {
    OsiIndexedVector r;    
    assert( r.getNumElements()==0 );
    
    // Test setting/getting elements with int* & float* vectors
    r.setVector( ne, inx, el );
    assert( r.getNumElements()==ne );
    for ( i=0; i<ne; i++ ) {
      assert( r.getIndices()[i]  == inx[i] );
      assert( r[inx[i]]  == el[i] );
      assert( r.getElements()[i] == el[i]  );
    }
    assert ( r.getMaxIndex()==7 );
    assert ( r.getMinIndex()==1 );
    
    // Test setting/getting elements with indices out of order  
    const int ne2 = 5;
    int inx2[ne2] = { 2, 4, 8, 14, 3 };
    double el2[ne2] = { 2.2, 4.4, 6.6, 8.8, 3.3 };
    
    r.setVector(ne2,inx2,el2);
    
    assert( r.getNumElements()==ne2 );    
    
    assert( r.getIndices()[0]==inx2[0] );
    assert( r.getElements()[0]==el2[0] );
    
    assert( r.getIndices()[1]==inx2[1] );
    assert( r.getElements()[1]==el2[1] );
    
    assert( r.getIndices()[2]==inx2[2] );
    assert( r.getElements()[2]==el2[2] );
    
    assert( r.getIndices()[3]==inx2[3] );
    assert( r.getElements()[3]==el2[3] );
    
    assert( r.getIndices()[4]==inx2[4] );
    assert( r.getElements()[4]==el2[4] );
    
    assert ( r.getMaxIndex()==14 );
    assert ( r.getMinIndex()==2 );
    assert ( r.getMaxIndex()==14 );
    assert ( r.getMinIndex()==2 );
    {
      bool errorThrown = false;
      try {
        r.duplicateIndex();
      }
      catch (CoinError e) {
        errorThrown = true;
      }
      assert( !errorThrown );
    }

    OsiIndexedVector r1(ne2,inx2,el2);
    assert( r == r1 );   
  }    
  OsiIndexedVector r;
  
  // Test operator[] where index is duplicated 
  // Causes exception to be thrown
  {
    const int ne3 = 4;
    int inx3[ne3] = { 2, 4, 2, 3 };
    double el3[ne3] = { 2.2, 4.4, 8.8, 6.6 };
    bool errorThrown = false;
    try {
      r.setVector(ne3,inx3,el3);
    }
    catch (CoinError e) {
      errorThrown = true;
    }
    assert( errorThrown );
    
  } 
  
  
  
  {
    OsiIndexedVector r;
    const int ne = 3;
    int inx[ne] = { 1, 2, 3 };
    double el[ne] = { 2.2, 4.4, 8.8};
    r.setVector(ne,inx,el);
    int c = r.capacity();
    int max = r.getMaxIndex();
    int min = r.getMinIndex();
    // Test swap function
    r.swap(0,2);
    assert( r.getIndices()[0]==3 );
    assert( r.getIndices()[1]==2 );
    assert( r.getIndices()[2]==1 );
    assert( r.getElements()[0]==8.8 );
    assert( r.getElements()[1]==4.4 );
    assert( r.getElements()[2]==2.2 );
    assert( r.getMaxIndex() == max );
    assert( r.getMinIndex() == min );
    assert( r.capacity() == c );
    
    // Test the append function
    OsiIndexedVector s;
    const int nes = 4;
    int inxs[nes] = { 11, 12, 13, 14 };
    double els[nes] = { .122, 14.4, 18.8, 19.9};
    s.setVector(nes,inxs,els);
    r.append(s);
    assert( r.getNumElements()==7 );
    assert( r.getIndices()[0]==3 );
    assert( r.getIndices()[1]==2 );
    assert( r.getIndices()[2]==1 );
    assert( r.getIndices()[3]==11 );
    assert( r.getIndices()[4]==12 );
    assert( r.getIndices()[5]==13 );
    assert( r.getIndices()[6]==14 );
    assert( r.getElements()[0]==8.8 );
    assert( r.getElements()[1]==4.4 );
    assert( r.getElements()[2]==2.2 );
    assert( r.getElements()[3]==.122 );
    assert( r.getElements()[4]==14.4 );    
    assert( r.getElements()[5]==18.8 );
    assert( r.getElements()[6]==19.9 );
    assert( r.getMaxIndex() == 14 );
    assert( r.getMinIndex() == 1 );
    
    // Test the resize function
    c = r.capacity();
    r.truncate(4);
    // we will lose 11
    assert( r.getNumElements()==3 );
    assert( r.getIndices()[0]==3 );
    assert( r.getIndices()[1]==2 );
    assert( r.getIndices()[2]==1 );
    assert( r.getElements()[0]==8.8 );
    assert( r.getElements()[1]==4.4 );
    assert( r.getElements()[2]==2.2 );
    assert( r.getMaxIndex() == 3 );
    assert( r.getMinIndex() == 1 );
    assert( r.capacity() == c );
  }
  
  // Test borrow and return vector
  {
    OsiIndexedVector r,r2;
    const int ne = 3;
    int inx[ne] = { 1, 2, 3 };
    double el[ne] = { 2.2, 4.4, 8.8};
    double els[4] = { 0.0,2.2, 4.4, 8.8};
    r.setVector(ne,inx,el);
    r2.borrowVector(4,ne,inx,els);
    assert (r==r2);
    r2.returnVector();
    assert (!r2.capacity());
    assert (!r2.getNumElements());
    assert (!r2.denseVector());
    assert (!r2.getIndices());
  }
  
  // Test copy constructor and assignment operator
  {
    OsiIndexedVector rhs;
    {
      OsiIndexedVector r;
      {
	OsiIndexedVector rC1(r);      
	assert( 0==r.getNumElements() );
	assert( 0==rC1.getNumElements() );
	
	
	r.setVector( ne, inx, el ); 
	
	assert( ne==r.getNumElements() );
	assert( 0==rC1.getNumElements() ); 
      }
      
      OsiIndexedVector rC2(r);   
      
      assert( ne==r.getNumElements() );
      assert( ne==rC2.getNumElements() );
      
      for ( i=0; i<ne; i++ ) {
	assert( r.getIndices()[i] == rC2.getIndices()[i] );
	assert( r[r.getIndices()[i]] == rC2.getElements()[i]);
	assert( r.getElements()[i] == rC2.getElements()[i] );
      }
      
      rhs=rC2;
    }
    // Test that rhs has correct values even though lhs has gone out of scope
    assert( rhs.getNumElements()==ne );
    
    for ( i=0; i<ne; i++ ) {
      assert( inx[i] == rhs.getIndices()[i] );
      assert(  el[i] == rhs.getElements()[i] );
    } 
  }
  
  // Test operator==
  {
    OsiIndexedVector v1,v2;
    assert( v1==v2 );
    assert( v2==v1 );
    assert( v1==v1 );
    assert( !(v1!=v2) );
    assert( v1.isEquivalent(v2) );
    assert( v2.isEquivalent(v1) );
    assert( v1.isEquivalent(v1) );
    
    v1.setVector( ne, inx, el );
    assert ( !(v1==v2) );
    assert ( v1!=v2 );
    assert( !v1.isEquivalent(v2) );
    
    OsiIndexedVector v3(v1);
    assert( v3==v1 );
    assert( v3!=v2 );
    assert( v1.isEquivalent(v3) );
    assert( v3.isEquivalent(v1) );
    
    OsiIndexedVector v4(v2);
    assert( v4!=v1 );
    assert( v4==v2 );
  }
  
  {
    // Test sorting of indexed vectors    
    const int ne = 4;
    int inx[ne] = { 1, 4, 0, 2 };
    double el[ne] = { 10., 40., 1., 20. };
    OsiIndexedVector r;
    r.setVector(ne,inx,el);
    
    // Test that indices are in increasing order
    r.sort();
    for ( i=1; i<ne; i++ ) assert( r.getIndices()[i-1] < r.getIndices()[i] );

  }    
  {
    // Test operator[] and indexExists()
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiIndexedVector r;
    bool errorThrown = false;
    try {
      assert( r[1]==0. );
    }
    catch (CoinError e) {
      errorThrown = true;
    }
    assert( errorThrown );
    
    r.setVector(ne,inx,el);
    
    errorThrown = false;
    try {
      assert( r[-1]==0. );
    }
    catch (CoinError e) {
      errorThrown = true;
    }
    assert( errorThrown );
    
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);
    errorThrown = false;
    try {
      assert( r[5]==0. );
    }
    catch (CoinError e) {
      errorThrown = true;
    }
    assert( errorThrown );
    
    assert(  r.isExistingIndex(2) );
    assert( !r.isExistingIndex(3) );
    
    assert ( r.getMaxIndex()==4 );
    assert ( r.getMinIndex()==0 );
  }
  
  // Test that attemping to get min/max index of a 0,
  // length vector 
  {
    OsiIndexedVector nullVec;
    assert( nullVec.getMaxIndex() ==
	    /*std::numeric_limits<int>::max()*/INT_MIN/*0*/ );
    assert( nullVec.getMinIndex() ==
	    /*std::numeric_limits<int>::min()*/INT_MAX/*0*/ );
  } 
  
  // Test OsiFltEq with equivalent method
  {    
    const int ne = 4;
    int inx1[ne] = { 1, 3, 4, 7 };
    double el1[ne] = { 1.2, 3.4, 5.6, 7.8 };
    int inx2[ne] = { 7, 4, 3, 1 };
    double el2[ne] = { 7.8+.5, 5.6+.5, 3.4+.5, 1.2+.5 };
    OsiIndexedVector v1,v2;
    v1.setVector(ne,inx1,el1);
    v2.setVector(ne,inx2,el2);
    assert( !v1.isEquivalent(v2) );
    assert(  v1.isEquivalent(v2,OsiAbsFltEq(.6)) );
    assert(  v1.isEquivalent(v2,OsiRelFltEq(.6)) );
  }
  
  {
    // Test reserve
    OsiIndexedVector v1,v2;
    assert( v1.capacity()==0 );
    v1.reserve(6);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==0 );
    v2=v1;
    assert( v2.capacity() == 6 );
    assert( v2.getNumElements()==0 );
    assert( v2==v1 );
    v1.setVector(0,NULL,NULL);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==0 );
    assert( v2==v1 );
    v2=v1;
    assert( v2.capacity() == 6 );
    assert( v2.getNumElements()==0 );
    assert( v2==v1 );
    
    const int ne = 2;
    int inx[ne] = { 1, 3 };
    double el[ne] = { 1.2, 3.4 };
    v1.setVector(ne,inx,el);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==2 );
    v2=v1;
    assert( v2.capacity()==6 );
    assert( v2.getNumElements()==2 );
    assert( v2==v1 );
    
    const int ne1 = 5;
    int inx1[ne1] = { 1, 3, 4, 5, 6 };
    double el1[ne1] = { 1.2, 3.4, 5., 6., 7. };
    v1.setVector(ne1,inx1,el1);
    assert( v1.capacity()==7 );
    assert( v1.getNumElements()==5 );
    v2=v1;
    assert( v2.capacity()==7 );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );
    
    const int ne2 = 8;
    int inx2[ne2] = { 1, 3, 4, 5, 6, 7, 8, 9 };
    double el2[ne2] = { 1.2, 3.4, 5., 6., 7., 8., 9., 10. };
    v1.setVector(ne2,inx2,el2);
    assert( v1.capacity()==10 );
    assert( v1.getNumElements()==8 );
    v2=v1;
    assert( v2.getNumElements()==8 );    
    assert( v2==v1 );
    
    v1.setVector(ne1,inx1,el1);
    assert( v1.capacity()==10 );
    assert( v1.getNumElements()==5 );
    v2=v1;    
    assert( v2.capacity()==10 );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );
    
    v1.reserve(7);
    assert( v1.capacity()==10 );
    assert( v1.getNumElements()==5 );
    v2=v1;
    assert( v2.capacity()==10 );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );
    
  }
  
  // Test the insert method
  {
    OsiIndexedVector v1;
    assert( v1.getNumElements()==0 );
    assert( v1.capacity()==0 );
    
    assert( !v1.isExistingIndex(1) );
    v1.insert(1,1.);
    assert( v1.getNumElements()==1 );
    assert( v1.capacity()==2 );
    assert( v1.getIndices()[0] == 1 );
    assert( v1.getElements()[0] == 1. );
    assert( v1.isExistingIndex(1) );
    
    assert( !v1.isExistingIndex(10) );
    v1.insert(10,10.);
    assert( v1.getNumElements()==2 );
    assert( v1.capacity()==11 );
    assert( v1.getIndices()[1] == 10 );
    assert( v1.getElements()[1] == 10. );
    assert( v1.isExistingIndex(1) );
    assert( v1.isExistingIndex(10) );
    
    assert( !v1.isExistingIndex(20) );
    v1.insert(20,20.);
    assert( v1.getNumElements()==3 );
    assert( v1.capacity()==21 );
    assert( v1.getIndices()[2] == 20 );
    assert( v1.getElements()[2] == 20. );
    assert( v1.isExistingIndex(20) );
    
    assert( !v1.isExistingIndex(30) );
    v1.insert(30,30.);
    assert( v1.getNumElements()==4 );
    assert( v1.capacity()==31 );
    assert( v1.getIndices()[3] == 30 );
    assert( v1.getElements()[3] == 30. );
    assert( v1.isExistingIndex(30) );
    
    assert( !v1.isExistingIndex(40) );
    v1.insert(40,40.);
    assert( v1.getNumElements()==5 );
    assert( v1.capacity()==41 );
    assert( v1.getIndices()[4] == 40 );
    assert( v1.getElements()[4] == 40. );
    assert( v1.isExistingIndex(40) );
    
    assert( !v1.isExistingIndex(50) );
    v1.insert(50,50.);
    assert( v1.getNumElements()==6 );
    assert( v1.capacity()==51 );
    assert( v1.getIndices()[5] == 50 );
    assert( v1.getElements()[5] == 50. );
    assert( v1.isExistingIndex(50) );
    
    OsiIndexedVector v2;
    const int ne1 = 3;
    int inx1[ne1] = { 1, 3, 4 };
    double el1[ne1] = { 1.2, 3.4, 5. };
    v2.setVector(ne1,inx1,el1);    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==5 );

    // Test clean method - get rid of 1.2
    assert(v2.clean(3.0)==2);
    assert(v2.denseVector()[1]==0.0);

    // Below are purely for debug - so use assert
    // so we won't try with false
    // Test checkClean 
    v2.checkClean();
    assert( v2.getNumElements()==2 );

    // Get rid of all
    assert(v2.clean(10.0)==0);
    v2.checkClear();
    
  }
  
  {
    //Test setConstant and setElement     
    OsiIndexedVector v2;
    const int ne1 = 3;
    int inx1[ne1] = { 1, 3, 4 };
    v2.setConstant(ne1,inx1,3.14);    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==5 );
    assert( v2.getIndices()[0]==1 );
    assert( v2.getElements()[0]==3.14 );
    assert( v2.getIndices()[1]==3 );
    assert( v2.getElements()[1]==3.14 );
    assert( v2.getIndices()[2]==4 );
    assert( v2.getElements()[2]==3.14 );
    
    assert( v2[3] == 3.14 );
    
    OsiIndexedVector v2X(ne1,inx1,3.14);
    assert( v2 == v2X );
    
  }
  
  {
    //Test setFull 
    OsiIndexedVector v2;
    const int ne2 = 3;
    double el2[ne2] = { 1., 3., 4. };
    v2.setFull(ne2,el2);    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==3 );
    assert( v2.getIndices()[0]==0 );
    assert( v2.getElements()[0]==1. );
    assert( v2.getIndices()[1]==1 );
    assert( v2.getElements()[1]==3. );
    assert( v2.getIndices()[2]==2 );
    assert( v2.getElements()[2]==4. );
    
    assert( v2[1] == 3. );
    
    OsiIndexedVector v2X(ne2,el2);
    assert( v2 == v2X ); 
    
    v2.setFull(0,el2); 
    assert( v2[2] == 0. );  
  }
  
  
#if 0
  // what happens when someone sets 
  // the number of elements to be a negative number
  {    
    const int ne = 4;
    int inx1[ne] = { 1, 3, 4, 7 };
    double el1[ne] = { 1.2, 3.4, 5.6, 7.8 };
    OsiIndexedVector v1;
    v1.set(-ne,inx1,el1);
  }
#endif
  
  
  
  // Test sum
  { 
    OsiIndexedVector s;
    assert( s.sum() == 0 );
    
    s.insert(25,45.);
    assert(s.sum()==45.);
    
    const int ne1 = 5;
    int inx1[ne1]   = { 10,  3,  4,  7,  5  };
    double el1[ne1] = { 1., 5., 6., 2., 9. };
    s.setVector(ne1,inx1,el1);
    
    assert(s.sum()==1.+5.+6.+2.+9.);
  }
  
  // Just another interesting test
  {    
    // Create numerator vector
    const int ne1 = 2;
    int inx1[ne1]   = { 1,  4  };
    double el1[ne1] = { 1., 6. };
    OsiIndexedVector v1(ne1,inx1,el1);
    
    // create denominator vector
    const int ne2 = 3;
    int inx2[ne2] =   { 1,  2,  4 };
    double el2[ne2] = { 1., 7., 4.};
    OsiIndexedVector v2(ne2,inx2,el2);
    
    // Compute ratio
    OsiIndexedVector ratio = v1 / v2;
    
/*
  The original code here used sortIncrElement, ostensibly to test that the
  zero (nominally in ratio[2]) disappeared. In fact, it's never created, and
  sortIncr was a noop. Changed to sortDecr to see some action.  Then it
  turned out that the Sun CC compiler is math-challenged when it comes to
  optimization, declaring that 1.0/1.0 = 24! Working through a temp gets it
  over this block. To say I'm developing an active dislike for C++
  optimization is a considerable understatement. -- lh, 02.04.09 --
*/
    // Sort ratios
    ratio.sortDecrElement();
    
    // Test that the sort really worked
    assert( ratio.getNumElements() == 2);
    double temp = 1.0 ;
    temp /= 1.0 ;
    assert( ratio.getElements()[1] == temp );
    temp = 6.0 ;
    temp /= 4.0 ;
    assert( ratio.getElements()[0] == temp );
    
    // Get numerator of of sorted ratio vector
    assert( v1[ ratio.getIndices()[1] ] == 1.0 );
    assert( v1[ ratio.getIndices()[0] ] == 6.0 );
    
    // Get denominator of of sorted ratio vector
    assert( v2[ ratio.getIndices()[1] ] == 1.0 );
    assert( v2[ ratio.getIndices()[0] ] == 4.0 );
  }
  
  // Test copy constructor from ShallowPackedVector
  {
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiIndexedVector std(ne,inx,el);
    OsiShallowPackedVector * spvP = new OsiShallowPackedVector(ne,inx,el);
    OsiIndexedVector pv(*spvP);
    assert( pv == std );
    assert( pv.isEquivalent(std) );
    delete spvP;
    assert( pv == std );
    assert( pv.isEquivalent(std) );
    pv.sortIncrElement();
    assert( pv != std );
    assert( pv.isEquivalent(std) );
  }
  
  // Test assignment from ShallowPackedVector
  {
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiIndexedVector std(ne,inx,el);
    OsiShallowPackedVector * spvP = new OsiShallowPackedVector(ne,inx,el);
    OsiIndexedVector pv;
    pv = *spvP;
    assert( pv == std );
    assert( pv.isEquivalent(std) );
    delete spvP;
    assert( pv == std );
    assert( pv.isEquivalent(std) );
    pv.sortIncrElement();
    assert( pv != std );
    assert( pv.isEquivalent(std) );
  }
  
  {
    // Test that sample usage works
    
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiIndexedVector r(ne,inx,el);
    
    assert( r.getIndices()[0]== 1  );
    assert( r.getElements()[0]==10. );
    assert( r.getIndices()[1]== 4  );
    assert( r.getElements()[1]==40. );
    assert( r.getIndices()[2]== 0  );
    assert( r.getElements()[2]== 1. );
    assert( r.getIndices()[3]== 2  );
    assert( r.getElements()[3]==50. );
    
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);
    
    r.sortIncrElement();
    
    assert( r.getIndices()[0]== 0  );
    assert( r.getElements()[0]== 1. );
    assert( r.getIndices()[1]== 1  );
    assert( r.getElements()[1]==10. );
    assert( r.getIndices()[2]== 4  );
    assert( r.getElements()[2]==40. );
    assert( r.getIndices()[3]== 2  );
    assert( r.getElements()[3]==50. );    
    
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);
    
    OsiIndexedVector r1;
    r1=r;
    assert( r==r1 );
    assert( r.isEquivalent(r1) );
    r.sortDecrElement();
    assert( r!=r1 );
    assert( r.isEquivalent(r1) );
    
    OsiIndexedVector add = r + r1;
    assert( add[0] ==  1.+ 1. );
    assert( add[1] == 10.+10. );
    assert( add[2] == 50.+50. );
    assert( add[3] ==  0.+ 0. );
    assert( add[4] == 40.+40. );
    
    assert( r.sum() == 10.+40.+1.+50. );
  }
  
}
    


