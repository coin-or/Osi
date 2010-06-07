// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>

#include "CoinPragma.hpp"

#include "OsiUnitTests.hpp"

#include "OsiRowCutDebugger.hpp"

//--------------------------------------------------------------------------
// test cut debugger methods.
void
OsiRowCutDebuggerUnitTest(const OsiSolverInterface * baseSiP,
			  const std::string & mpsDir )
{
  
  CoinRelFltEq eq;
  
  // Test default constructor
  {
    OsiRowCutDebugger r;
    assert( r.integerVariable_==NULL );
    assert( r.optimalSolution_==NULL );
    assert( r.numberColumns_==0);
  }
  
  {
    // Get non trivial instance
    OsiSolverInterface * imP = baseSiP->clone();
    std::string fn = mpsDir+"exmip1";
    imP->readMps(fn.c_str(),"mps");
    //std::cerr <<imP->getNumRows() <<std::endl;
    assert( imP->getNumRows() == 5);
    
    // activate debugger
    imP->activateRowCutDebugger("ab cd /x/ /exmip1.asc");
    
    int i;
    
    // return debugger
    const OsiRowCutDebugger * debugger = imP->getRowCutDebugger();
    // check 
    assert (debugger!=NULL);
    
    assert (debugger->numberColumns_==8);
    
    const bool type[]={0,0,1,1,0,0,0,0};
    const double values[]= {2.5, 0, 1, 1, 0.5, 3, 0, 0.26315789473684253};
    CoinPackedVector objCoefs(8,imP->getObjCoefficients());
   
#if 0
    for (i=0;i<8;i++) {
      assert(type[i]==debugger->integerVariable_[i]);
      std::cerr <<i  <<" " <<values[i] <<" " <<debugger->optimalSolution_[i] <<std::endl;
    }
#endif
    
    double objValue = objCoefs.dotProduct(values);
    double debuggerObjValue = objCoefs.dotProduct(debugger->optimalSolution_);
    //std::cerr <<objValue <<" " <<debuggerObjValue <<std::endl;
    assert( eq(objValue,debuggerObjValue) );
    
    OsiRowCutDebugger rhs;
    {
      OsiRowCutDebugger rC1(*debugger);
      
      assert (rC1.numberColumns_==8);
      for (i=0;i<8;i++) {
        assert(type[i]==rC1.integerVariable_[i]);
      }
      assert( eq(objValue,objCoefs.dotProduct(rC1.optimalSolution_)) );
      
      rhs=rC1;
      assert (rhs.numberColumns_==8);
      for (i=0;i<8;i++) {
        assert(type[i]==rhs.integerVariable_[i]);
      }
      assert( eq(objValue,objCoefs.dotProduct(rhs.optimalSolution_)) );
    }
    // Test that rhs has correct values even though lhs has gone out of scope
    assert (rhs.numberColumns_==8);
    for (i=0;i<8;i++) {
      assert(type[i]==rhs.integerVariable_[i]);
    }
    assert( eq(objValue,objCoefs.dotProduct(rhs.optimalSolution_)) );
    OsiRowCut cut[2];
    
    const int ne = 3;
    int inx[ne] = { 0, 2, 3 };
    double el[ne] = { 1., 1., 1. };
    cut[0].setRow(ne,inx,el);
    cut[0].setUb(5.);
    
    el[1]=5;
    cut[1].setRow(ne,inx,el);
    cut[1].setUb(5);
    OsiCuts cs; 
    cs.insert(cut[0]);
    cs.insert(cut[1]);
    assert(!debugger->invalidCut(cut[0]));
    assert( debugger->invalidCut(cut[1]));
    assert(debugger->validateCuts(cs,0,2)==1);
    assert(debugger->validateCuts(cs,0,1)==0);
    delete imP;
  }
}
