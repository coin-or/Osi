// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#include <cassert>
#include <utility>
#include "CoinSort.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif


void
tripleCompareUnitTest()
{
#if 0
  OsiTriple<int, int, double > p1_1(1,100,1.);
  OsiTriple<int, int, double > p1_10(1,100,10.);
  OsiTriple<int, int, double > p1_N1(1,100,-1.);
  //OsiTriple<int, int, double > p1_N10(1,100,-10.);
  OsiTriple<int, int, double > pN1_1(-1,200,1.);
  OsiTriple<int, int, double > pN1_10(-1,200,10.);
  OsiTriple<int, int, double > pN1_N10(-1,200,-10.);
  OsiTriple<int, int, double > p10_1(10,5,-1.);

  OsiIncrIndexOrdered pcIio;
  assert( !pcIio(pN1_1,pN1_10) );
  assert(  pcIio(pN1_1,p10_1) );
  assert(  pcIio(pN1_10,p10_1) );

  OsiDecrIndexOrdered pcDio;
  assert( !pcDio(pN1_1,pN1_10) );
  assert( !pcDio(pN1_1,p10_1) );
  assert( !pcDio(pN1_10,p10_1) );
  assert(  pcDio(p10_1,p1_10) );

  OsiIncrElementOrdered pcIeo;
  assert(  pcIeo(pN1_1,pN1_10) );
  assert( !pcIeo(pN1_1,p10_1) );
  assert( !pcIeo(pN1_10,p10_1) );

  OsiDecrElementOrdered pcDeo;
  assert( !pcDeo(pN1_1,pN1_10) );
  assert(  pcDeo(pN1_1,p10_1) );
  assert(  pcDeo(pN1_10,p10_1) );

  OsiIncrAbsElementOrdered pcIaeo;
  assert(  pcIaeo(pN1_1,pN1_10) );
  assert( !pcIaeo(pN1_1,p10_1) );
  assert( !pcIaeo(pN1_10,p10_1) );
  assert(  pcIaeo(pN1_1,pN1_N10) );
  assert( !pcIaeo(pN1_1,p1_N1) );
  assert( !pcIaeo(pN1_10,p1_N1) );

  OsiDecrAbsElementOrdered pcDaeo;
  assert( !pcDaeo(pN1_1,pN1_10) );
  assert( !pcDaeo(pN1_1,p10_1) );
  assert(  pcDaeo(pN1_10,p10_1) );
  assert( !pcDaeo(pN1_1,pN1_N10) );
  assert( !pcDaeo(pN1_1,p1_N1) );
  assert(  pcDaeo(pN1_10,p1_N1) );
  
  
  OsiIncreasingExternalVectorOrdered<int,int,double,double> pcEvA;

  const int numElements =11;
  double sortElements[numElements]={10., 9., 8., 7., 6., 5., 4., 3., 2., 1., 0. };
  OsiIncreasingExternalVectorOrdered<int,int,double,double> pcIevB(sortElements);
  assert( !pcIevB(p1_1,p1_10) );
  assert( !pcIevB(p1_1,p10_1) );
  assert( !pcIevB(p1_10,p10_1) );
  assert( !pcIevB(p1_10,p1_1) );
  assert(  pcIevB(p10_1,p1_1) );
  assert(  pcIevB(p10_1,p1_10) );

  OsiIncreasingExternalVectorOrdered<int,int,double,double> pcIevC;
  pcIevC=pcIevB;
  assert( !pcIevC(p1_1,p1_10) );
  assert( !pcIevC(p1_1,p10_1) );
  assert( !pcIevC(p1_10,p10_1) );
  assert( !pcIevC(p1_10,p1_1) );
  assert(  pcIevC(p10_1,p1_1) );
  assert(  pcIevC(p10_1,p1_10) );

  OsiIncrSolutionOrdered pciso(sortElements);
  assert( !pciso(p1_1,p1_10) );
  assert( !pciso(p1_1,p10_1) );
  assert( !pciso(p1_10,p10_1) );
  assert( !pciso(p1_10,p1_1) );
  assert(  pciso(p10_1,p1_1) );
  assert(  pciso(p10_1,p1_10) );

  OsiDecrSolutionOrdered pcdso(sortElements);
  assert( !pcdso(p1_1,p1_10) );
  assert(  pcdso(p1_1,p10_1) );
  assert(  pcdso(p1_10,p10_1) );
  assert( !pcdso(p1_10,p1_1) );
  assert( !pcdso(p10_1,p1_1) );
  assert( !pcdso(p10_1,p1_10) );
#endif
}
