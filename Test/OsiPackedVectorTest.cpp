// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiPackedVector.hpp"
#include "OsiShallowPackedVector.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

//--------------------------------------------------------------------------
void
OsiPackedVectorUnitTest()
{

  int i;
  // Test default constructor
  {
    OsiPackedVector r;
    assert( r.indices_==NULL );
    assert( r.origIndices_==NULL );
    assert( r.elements_==NULL );
    assert( r.getNumElements()==0 );
    assert( r.capacity_==0);
  }

  // Test set and get methods
  const int ne = 4;
  int inx[ne] = { 1, 3, 4, 7 };
  double el[ne] = { 1.2, 3.4, 5.6, 7.8 };
  {
    OsiPackedVector r;    
    assert( r.getNumElements()==0 );
    
    // Test setting/getting elements with int* & float* vectors
    r.setVector( ne, inx, el );
    assert( r.getNumElements()==ne );
    for ( i=0; i<ne; i++ ) {
      assert( r.getIndices()[i]  == inx[i] );
      assert( r.getOriginalPosition()[i]  == i );
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
    assert( r.getOriginalPosition()[0]==0 );
    assert( r.getElements()[0]==el2[0] );
    
    assert( r.getIndices()[1]==inx2[1] );
    assert( r.getOriginalPosition()[1]==1 );
    assert( r.getElements()[1]==el2[1] );
    
    assert( r.getIndices()[2]==inx2[2] );
    assert( r.getOriginalPosition()[2]==2 );
    assert( r.getElements()[2]==el2[2] );
    
    assert( r.getIndices()[3]==inx2[3] );
    assert( r.getOriginalPosition()[3]==3 );
    assert( r.getElements()[3]==el2[3] );
    
    assert( r.getIndices()[4]==inx2[4] );
    assert( r.getOriginalPosition()[4]==4 );
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

    OsiPackedVector r1(ne2,inx2,el2);
    assert( r == r1 );   
    
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

      errorThrown = false;
      try {
	 r.duplicateIndex();
      }
      catch (CoinError e) {
	 errorThrown = true;
      }
      assert( errorThrown );
    } 



    // Test sorting by increasing elements
    r.setVector(ne2,inx2,el2);
    bool incr=true;
    for ( i=1; i<ne2; i++ )
       if ( r.getElements()[i-1]>r.getElements()[i] )
	  incr=false;
    assert( !incr );
    r.sortIncrElement();
    incr = true;
    for ( i=1; i<ne2; i++ )
       if ( r.getElements()[i-1]>r.getElements()[i] )
	  incr=false;
    assert( incr );

  } 
  

  {
    OsiPackedVector r;
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
    OsiPackedVector s;
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
    assert( r.getNumElements()==4 );
    assert( r.getIndices()[0]==3 );
    assert( r.getIndices()[1]==2 );
    assert( r.getIndices()[2]==1 );
    assert( r.getIndices()[3]==11 );
    assert( r.getElements()[0]==8.8 );
    assert( r.getElements()[1]==4.4 );
    assert( r.getElements()[2]==2.2 );
    assert( r.getElements()[3]==.122 );
    assert( r.getMaxIndex() == 11 );
    assert( r.getMinIndex() == 1 );
    assert( r.capacity() == c );
  }


  // Test copy constructor and assignment operator
  {
    OsiPackedVector rhs;
    {
      OsiPackedVector r;
      {
        OsiPackedVector rC1(r);      
        assert( 0==r.getNumElements() );
        assert( 0==rC1.getNumElements() );
        
        
        r.setVector( ne, inx, el ); 
        
        assert( ne==r.getNumElements() );
        assert( 0==rC1.getNumElements() ); 
      }
      
      OsiPackedVector rC2(r);   
      
      assert( ne==r.getNumElements() );
      assert( ne==rC2.getNumElements() );
      
      for ( i=0; i<ne; i++ ) {
        assert( r.getIndices()[i] == rC2.getIndices()[i] );
        assert( r.getOriginalPosition()[i] == rC2.getOriginalPosition()[i] );
        assert( r.getElements()[i] == rC2.getElements()[i] );
      }
      
      rhs=rC2;
    }
    // Test that rhs has correct values even though lhs has gone out of scope
    assert( rhs.getNumElements()==ne );
    
    for ( i=0; i<ne; i++ ) {
      assert( inx[i] == rhs.getIndices()[i] );
      assert( i == rhs.getOriginalPosition()[i] );
      assert(  el[i] == rhs.getElements()[i] );
    } 
  }

  // Test operator==
  {
    OsiPackedVector v1,v2;
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

    OsiPackedVector v3(v1);
    assert( v3==v1 );
    assert( v3!=v2 );
    assert( v1.isEquivalent(v3) );
    assert( v3.isEquivalent(v1) );

    OsiPackedVector v4(v2);
    assert( v4!=v1 );
    assert( v4==v2 );
  }

  {
    // Test sorting of packed vectors    
    const int ne = 4;
    int inx[ne] = { 1, 4, 0, 2 };
    double el[ne] = { 10., 40., 1., 20. };
    double extSortKey[5] = { -20., 10., 5., 4., 20. };
    OsiPackedVector r;
    r.setVector(ne,inx,el);
    
    // Test that indices are in increasing order
    r.sortIncrIndex();
    for ( i=1; i<ne; i++ ) assert( r.getIndices()[i-1] < r.getIndices()[i] );
    
    // Sort by element values;
    r.sortIncrElement();  
    for ( i=1; i<ne; i++ ) assert( r.getElements()[i-1] < r.getElements()[i] );

    // Sort using indices into an external sort vector
    CoinIncrSolutionOrdered pcSo(extSortKey);
    r.sort(pcSo);   
    for ( i=1; i<ne; i++ ) 
      assert( extSortKey[r.getIndices()[i-1]] < extSortKey[r.getIndices()[i]] );


    // Now sort by indices.
    r.sortIncrIndex();
    for ( i=1; i<ne; i++ ) assert( r.getIndices()[i-1] < r.getIndices()[i] );
  }


  {
    // Test operator[] and indexExists()
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiPackedVector r;
    assert( r[1]==0. );

    r.setVector(ne,inx,el);

    assert( r[-1]==0. );
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);
    assert( r[ 5]==0. );
    assert(  r.isExistingIndex(2) );
    assert( !r.isExistingIndex(3) );

    r.sortIncrElement();
    
    assert( r[-1]==0. );
    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);
    assert( r[ 5]==0. );
    assert( !r.isExistingIndex(-1) );
    assert(  r.isExistingIndex(0) );
    assert( !r.isExistingIndex(3) );
    assert(  r.isExistingIndex(4) );
    assert( !r.isExistingIndex(5) );
    
    assert ( r.getMaxIndex()==4 );
    assert ( r.getMinIndex()==0 );
  }
  
  // Test that attemping to get min/max index of a 0,
  // length vector 
  {
    OsiPackedVector nullVec;
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
    OsiPackedVector v1,v2;
    v1.setVector(ne,inx1,el1);
    v2.setVector(ne,inx2,el2);
    assert( !v1.isEquivalent(v2) );
    assert(  v1.isEquivalent(v2,OsiAbsFltEq(.6)) );
    assert(  v1.isEquivalent(v2,OsiRelFltEq(.6)) );
  }

  {
    // Test reserve
    OsiPackedVector v1,v2;
    assert( v1.capacity()==0 );
    v1.reserve(6);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==0 );
    v2=v1;
    assert( v2.capacity() == v2.getNumElements() );
    assert( v2.getNumElements()==0 );
    assert( v2==v1 );
    v1.setVector(0,NULL,NULL);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==0 );
    assert( v2==v1 );
    v2=v1;
    assert( v2.capacity() == v2.getNumElements() );
    assert( v2.getNumElements()==0 );
    assert( v2==v1 );

    const int ne = 2;
    int inx[ne] = { 1, 3 };
    double el[ne] = { 1.2, 3.4 };
    v1.setVector(ne,inx,el);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==2 );
    v2=v1;
    assert( v2.capacity()==v2.getNumElements() );
    assert( v2.getNumElements()==2 );
    assert( v2==v1 );

    const int ne1 = 5;
    int inx1[ne1] = { 1, 3, 4, 5, 6 };
    double el1[ne1] = { 1.2, 3.4, 5., 6., 7. };
    v1.setVector(ne1,inx1,el1);
    assert( v1.capacity()==6 );
    assert( v1.getNumElements()==5 );
    v2=v1;
    assert( v2.capacity()==v2.getNumElements() );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );

    const int ne2 = 8;
    int inx2[ne2] = { 1, 3, 4, 5, 6, 7, 8, 9 };
    double el2[ne2] = { 1.2, 3.4, 5., 6., 7., 8., 9., 10. };
    v1.setVector(ne2,inx2,el2);
    assert( v1.capacity()==8 );
    assert( v1.getNumElements()==8 );
    v2=v1;
    assert( v2.capacity()==v2.getNumElements() );
    assert( v2.getNumElements()==8 );    
    assert( v2==v1 );
 
    v1.setVector(ne1,inx1,el1);
    assert( v1.capacity()==8 );
    assert( v1.getNumElements()==5 );
    v2=v1;    
    assert( v2.capacity()==8 );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );

    v1.reserve(7);
    assert( v1.capacity()==8 );
    assert( v1.getNumElements()==5 );
    v2=v1;
    assert( v2.capacity()==8 );
    assert( v2.getNumElements()==5 );
    assert( v2==v1 );

  }

  // Test the insert method
  {
    OsiPackedVector v1;
    assert( v1.getNumElements()==0 );
    assert( v1.capacity()==0 );

    assert( !v1.isExistingIndex(1) );
    v1.insert(1,1.);
    assert( v1.getNumElements()==1 );
    assert( v1.capacity()==5 );
    assert( v1.getIndices()[0] == 1 );
    assert( v1.getElements()[0] == 1. );
    assert( v1.isExistingIndex(1) );

    assert( !v1.isExistingIndex(10) );
    v1.insert(10,10.);
    assert( v1.getNumElements()==2 );
    assert( v1.capacity()==5 );
    assert( v1.getIndices()[1] == 10 );
    assert( v1.getElements()[1] == 10. );
    assert( v1.isExistingIndex(1) );
    assert( v1.isExistingIndex(10) );

    assert( !v1.isExistingIndex(20) );
    v1.insert(20,20.);
    assert( v1.getNumElements()==3 );
    assert( v1.capacity()==5 );
    assert( v1.getIndices()[2] == 20 );
    assert( v1.getElements()[2] == 20. );
    assert( v1.isExistingIndex(20) );

    assert( !v1.isExistingIndex(30) );
    v1.insert(30,30.);
    assert( v1.getNumElements()==4 );
    assert( v1.capacity()==5 );
    assert( v1.getIndices()[3] == 30 );
    assert( v1.getElements()[3] == 30. );
    assert( v1.isExistingIndex(30) );

    assert( !v1.isExistingIndex(40) );
    v1.insert(40,40.);
    assert( v1.getNumElements()==5 );
    assert( v1.capacity()==5 );
    assert( v1.getIndices()[4] == 40 );
    assert( v1.getElements()[4] == 40. );
    assert( v1.isExistingIndex(40) );
    
    assert( !v1.isExistingIndex(50) );
    v1.insert(50,50.);
    assert( v1.getNumElements()==6 );
    assert( v1.capacity()==10 );
    assert( v1.getIndices()[5] == 50 );
    assert( v1.getElements()[5] == 50. );
    assert( v1.isExistingIndex(50) );
    
    OsiPackedVector v2;
    const int ne1 = 3;
    int inx1[ne1] = { 1, 3, 4 };
    double el1[ne1] = { 1.2, 3.4, 5. };
    v2.setVector(ne1,inx1,el1);    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==3 );
    
  }

  {
    //Test setConstant and setElement     
    OsiPackedVector v2;
    const int ne1 = 3;
    int inx1[ne1] = { 1, 3, 4 };
    v2.setConstant(ne1,inx1,3.14);    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==3 );
    assert( v2.getIndices()[0]==1 );
    assert( v2.getElements()[0]==3.14 );
    assert( v2.getIndices()[1]==3 );
    assert( v2.getElements()[1]==3.14 );
    assert( v2.getIndices()[2]==4 );
    assert( v2.getElements()[2]==3.14 );

    assert( v2[3] == 3.14 );

    OsiPackedVector v2X(ne1,inx1,3.14);
    assert( v2 == v2X );

#if 0
    v2.setElement( 1 , 100. );    
    assert( v2.getNumElements()==3 );
    assert( v2.capacity()==3 );
    assert( v2.getIndices()[0]==1 );
    assert( v2.getElements()[0]==3.14 );
    assert( v2.getIndices()[1]==3 );
    assert( v2.getElements()[1]==100. );
    assert( v2.getIndices()[2]==4 );
    assert( v2.getElements()[2]==3.14 );
    
    assert( v2[3] == 100. );
    
    bool errorThrown = false;
    try {
      v2.setElement( 100, 100. );
    }
    catch (CoinError e) {
      errorThrown = true;
    }
    assert( errorThrown );
#endif
  }

  {
    //Test setFull 
    OsiPackedVector v2;
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

    OsiPackedVector v2X(ne2,el2);
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
    OsiPackedVector v1;
    v1.set(-ne,inx1,el1);
  }
#endif

  
  // Test adding vectors
  {    
    const int ne1 = 5;
    int inx1[ne1]   = { 1,  3,  4,  7,  5  };
    double el1[ne1] = { 1., 5., 6., 2., 9. };
    const int ne2 = 4;
    int inx2[ne2] =   { 7,  4,  2,  1  };
    double el2[ne2] = { 7., 4., 2., 1. };
    OsiPackedVector v1;
    v1.setVector(ne1,inx1,el1);
    OsiPackedVector v2;
    v2.setVector(ne2,inx2,el2);
    OsiPackedVector r = v1 + v2;

    const int ner = 6;
    int inxr[ner] =   {    1,     2,     3,     4,     5,     7  };
    double elr[ner] = { 1.+1., 0.+2., 5.+0., 6.+4., 9.+0., 2.+7. };
    OsiPackedVector rV;
    rV.setVector(ner,inxr,elr);
    assert( r.isEquivalent(rV) );
  } 
  
  // Test subtracting vectors
  {    
    const int ne1 = 5;
    int inx1[ne1]   = { 1,  3,  4,  7,  5  };
    double el1[ne1] = { 1., 5., 6., 2., 9. };
    const int ne2 = 4;
    int inx2[ne2] =   { 7,  4,  2,  1  };
    double el2[ne2] = { 7., 4., 2., 1. };
    OsiPackedVector v1;
    v1.setVector(ne1,inx1,el1);
    OsiPackedVector v2;
    v2.setVector(ne2,inx2,el2);
    OsiPackedVector r = v1 - v2;

    const int ner = 6;
    int inxr[ner] =   {    1,     2,     3,     4,     5,     7  };
    double elr[ner] = { 1.-1., 0.-2., 5.-0., 6.-4., 9.-0., 2.-7. };
    OsiPackedVector rV;
    rV.setVector(ner,inxr,elr);
    assert( r.isEquivalent(rV) );
  } 
  
  // Test multiplying vectors
  {    
    const int ne1 = 5;
    int inx1[ne1]   = { 1,  3,  4,  7,  5  };
    double el1[ne1] = { 1., 5., 6., 2., 9. };
    const int ne2 = 4;
    int inx2[ne2] =   { 7,  4,  2,  1  };
    double el2[ne2] = { 7., 4., 2., 1. };
    OsiPackedVector v1;
    v1.setVector(ne1,inx1,el1);
    OsiPackedVector v2;
    v2.setVector(ne2,inx2,el2);
    OsiPackedVector r = v1 * v2;

    const int ner = 6;
    int inxr[ner] =   {    1,     2,     3,     4,     5,     7  };
    double elr[ner] = { 1.*1., 0.*2., 5.*0., 6.*4., 9.*0., 2.*7. };
    OsiPackedVector rV;
    rV.setVector(ner,inxr,elr);
    assert( r.isEquivalent(rV) );
  } 
  
  // Test dividing vectors
  {    
    const int ne1 = 3;
    int inx1[ne1]   = { 1,  4,  7  };
    double el1[ne1] = { 1., 6., 2. };
    const int ne2 = 4;
    int inx2[ne2] =   { 7,  4,  2,  1  };
    double el2[ne2] = { 7., 4., 2., 1. };
    OsiPackedVector v1;
    v1.setVector(ne1,inx1,el1);
    OsiPackedVector v2;
    v2.setVector(ne2,inx2,el2);
    OsiPackedVector r = v1 / v2;

    const int ner = 4;
    int inxr[ner] =   {    1,     2,      4,     7  };
    double elr[ner] = { 1./1., 0./2.,  6./4., 2./7. };
    OsiPackedVector rV;
    rV.setVector(ner,inxr,elr);
    assert( r.isEquivalent(rV) );
  }
   
  // Test sum
  { 
    OsiPackedVector s;
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
    OsiPackedVector v1(ne1,inx1,el1);

    // create denominator vector
    const int ne2 = 3;
    int inx2[ne2] =   { 1,  2,  4 };
    double el2[ne2] = { 1., 7., 4.};
    OsiPackedVector v2(ne2,inx2,el2);

    // Compute ratio
    OsiPackedVector ratio = v1 / v2;

    // Sort ratios
    ratio.sortIncrElement();

    // Test that the sort really worked
    assert( ratio.getElements()[0] == 0.0/7.0 );
    assert( ratio.getElements()[1] == 1.0/1.0 );
    assert( ratio.getElements()[2] == 6.0/4.0 );

    // Get numerator of of sorted ratio vector
    assert( v1[ ratio.getIndices()[0] ] == 0.0 );
    assert( v1[ ratio.getIndices()[1] ] == 1.0 );
    assert( v1[ ratio.getIndices()[2] ] == 6.0 );

    // Get denominator of of sorted ratio vector
    assert( v2[ ratio.getIndices()[0] ] == 7.0 );
    assert( v2[ ratio.getIndices()[1] ] == 1.0 );
    assert( v2[ ratio.getIndices()[2] ] == 4.0 );
  }

  // Test copy constructor from ShallowPackedVector
  {
    const int ne = 4;
    int inx[ne] =   {  1,   4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiPackedVector std(ne,inx,el);
    OsiShallowPackedVector * spvP = new OsiShallowPackedVector(ne,inx,el);
    OsiPackedVector pv(*spvP);
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
    OsiPackedVector std(ne,inx,el);
    OsiShallowPackedVector * spvP = new OsiShallowPackedVector(ne,inx,el);
    OsiPackedVector pv;
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
    OsiPackedVector r(ne,inx,el);

    assert( r.getIndices()[0]== 1  );
    assert( r.getElements()[0]==10. );
    assert( r.getIndices()[1]== 4  );
    assert( r.getElements()[1]==40. );
    assert( r.getIndices()[2]== 0  );
    assert( r.getElements()[2]== 1. );
    assert( r.getIndices()[3]== 2  );
    assert( r.getElements()[3]==50. );

    assert( r.getOriginalPosition()[0]==0 );
    assert( r.getOriginalPosition()[1]==1 );
    assert( r.getOriginalPosition()[2]==2 );
    assert( r.getOriginalPosition()[3]==3 );

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
    
    assert( r.getOriginalPosition()[0]==2 );
    assert( r.getOriginalPosition()[1]==0 );
    assert( r.getOriginalPosition()[2]==1 );
    assert( r.getOriginalPosition()[3]==3 );

    assert( r[ 0]==1. );
    assert( r[ 1]==10.);
    assert( r[ 2]==50.);
    assert( r[ 3]==0. );
    assert( r[ 4]==40.);

    r.sortOriginalOrder();
    
    assert( r.getIndices()[0]== 1  );
    assert( r.getElements()[0]==10. );
    assert( r.getIndices()[1]== 4  );
    assert( r.getElements()[1]==40. );
    assert( r.getIndices()[2]== 0  );
    assert( r.getElements()[2]== 1. );
    assert( r.getIndices()[3]== 2  );
    assert( r.getElements()[3]==50. );

    OsiPackedVector r1;
    r1=r;
    assert( r==r1 );
    assert( r.isEquivalent(r1) );
    r.sortIncrElement();
    assert( r!=r1 );
    assert( r.isEquivalent(r1) );

    OsiPackedVector add = r + r1;
    assert( add[0] ==  1.+ 1. );
    assert( add[1] == 10.+10. );
    assert( add[2] == 50.+50. );
    assert( add[3] ==  0.+ 0. );
    assert( add[4] == 40.+40. );

    assert( r.sum() == 10.+40.+1.+50. );
  }
    
  {
    // Test findIndex
    const int ne = 4;
    int inx[ne] =   {  1,  -4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiPackedVector r(ne,inx,el);

    assert( r.findIndex(2)  == 3 );
    assert( r.findIndex(0)  == 2 );
    assert( r.findIndex(-4) == 1 );
    assert( r.findIndex(1)  == 0 );
    assert( r.findIndex(3)  == -1 );
  }  
#if 0
  {
    // Test construction with testing for duplicates as false
    const int ne = 4;
    int inx[ne] =   {  1,  -4,  0,   2 };
    double el[ne] = { 10., 40., 1., 50. };
    OsiPackedVector rT(ne,inx,el);
    OsiPackedVector r(false);

    assert( !r.isExistingIndex(1) );
    r.insert(1,10.);

    assert( !r.isExistingIndex(-4) );
    r.insert(-4,20.);

    assert( !r.isExistingIndex(0) );
    r.insert(0,1.);
    
    assert( r.isExistingIndex(-4) );  // This line is failing!!
                                      // If r is constructed with true,
                                      // then it does not fail
    int neg4Index = r.findIndex(-4);
    assert( neg4Index == 1 );
    r.setElement(neg4Index, r.getElements()[neg4Index] + 20);

    assert( !r.isExistingIndex(2) );
    r.insert(2,50.);

    assert( r == rT );
  }
#endif
}

