// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// Test individual classes or groups of classes

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>

#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSimplexInterface.hpp"
#include "OsiRowCutDebugger.hpp"
#ifdef COIN_USE_OSL
#include "OsiOslSolverInterface.hpp"
#endif
#ifdef COIN_USE_XPR
#include "OsiXprSolverInterface.hpp"
#endif
#ifdef COIN_USE_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#ifdef COIN_USE_SPX
#include "OsiSpxSolverInterface.hpp"
#endif
#ifdef COIN_USE_VOL
#include "OsiVolSolverInterface.hpp"
#endif
#ifdef COIN_USE_DYLP
#include "OsiDylpSolverInterface.hpp"
#endif
#ifdef COIN_USE_GLPK
#include "OsiGlpkSolverInterface.hpp"
#endif
#ifdef COIN_USE_CLP
#include "OsiClpSolverInterface.hpp"
#endif
// Function Prototypes. Function definitions is in this file.
void testingMessage( const char * const msg );

//----------------------------------------------------------------
// unitTest [-mpsDir=V1] [-netlibDir=V2] [-skipOsiSolverInterface]
// 
// where:
//   -mpsDir: directory containing mps test files
//       Default value V1="../Mps/Sample"    
//   -netlibDir: directory containing netlib files
//       Default value V2="../Mps/Netlib"
//   -testOsiSolverInterface
//       If specified, then OsiSolveInterface::unitTest
//       is skipped over and not run.
//
// All parameters are optional.
//----------------------------------------------------------------

int main (int argc, const char *argv[])
{
  int i;

#ifdef COIN_USE_XPR
  OsiXprSolverInterface::setLogFileName("xprCallTrace.txt");
#endif

  // define valid parameter keywords
  std::set<std::string> definedKeyWords;
  definedKeyWords.insert("-mpsDir");
  definedKeyWords.insert("-netlibDir");
  definedKeyWords.insert("-testOsiSolverInterface");

  // Create a map of parmater keys and associated data
  std::map<std::string,std::string> parms;
  for ( i=1; i<argc; i++ ) {
    std::string parm(argv[i]);
    std::string key,value;
    unsigned int  eqPos = parm.find('=');

    // Does parm contain and '='
    if ( eqPos==std::string::npos ) {
      //Parm does not contain '='
      key = parm;
    }
    else {
      key=parm.substr(0,eqPos);
      value=parm.substr(eqPos+1);
    }

    // Is specifed key valid?
    if ( definedKeyWords.find(key) == definedKeyWords.end() ) {
      // invalid key word.
      // Write help text
      std::cerr <<"Undefined parameter \"" <<key <<"\".\n";
      std::cerr <<"Correct usage: \n";
      std::cerr <<"  unitTest [-mpsDir=V1] [-netlibDir=V2] [-testOsiSolverInterface]\n";
      std::cerr <<"  where:\n";
      std::cerr <<"    -mpsDir: directory containing mps test files\n";
      std::cerr <<"        Default value V1=\"../Mps/Sample\"\n";
      std::cerr <<"    -netlibDir: directory containing netlib files\n";
      std::cerr <<"        Default value V2=\"../Mps/Netlib\"\n";
      std::cerr <<"    -testOsiSolverInterface\n";
      std::cerr <<"        If specified, then OsiSolveInterface::unitTest\n";
      std::cerr <<"        is run.\n";
      return 1;
    }
    parms[key]=value;
  }
  
  const char dirsep =  CoinFindDirSeparator();
  // Set directory containing mps data files.
  std::string mpsDir;
  if (parms.find("-mpsDir") != parms.end())
    mpsDir=parms["-mpsDir"] + dirsep;
  else 
    mpsDir = dirsep == '/' ? "../Mps/Sample/" : "..\\Mps\\Sample\\";
 
  // Set directory containing netlib data files.
  std::string netlibDir;
  if (parms.find("-netlibDir") != parms.end())
    netlibDir=parms["-netlibDir"] + dirsep;
  else 
    netlibDir = dirsep == '/' ? "../Mps/Netlib/" : "..\\Mps\\Netlib\\";

#ifdef COIN_USE_OSL  
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing OsiRowCut with OsiOslSolverInterface\n" );
    OsiRowCutUnitTest(&oslSi,mpsDir);
  }
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing OsiColCut with OsiOslSolverInterface\n" );
    OsiColCutUnitTest(&oslSi,mpsDir);
  }
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiOslSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&oslSi,mpsDir);
  }
  {
    OsiOslSolverInterface oslSi;
    testingMessage( "Testing OsiSimplexInterface with OsiOslSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&oslSi,mpsDir);
  }

#endif

#ifdef COIN_USE_XPR  
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiRowCut with OsiXprSolverInterface\n" );
    OsiRowCutUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiColCut with OsiXprSolverInterface\n" );
    OsiColCutUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiXprSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&xprSi,mpsDir);
  }
  {
    OsiXprSolverInterface xprSi;
    testingMessage( "Testing OsiSimplexInterface with OsiXprSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&xprSi,mpsDir);
  }
#endif

#ifdef COIN_USE_CPX
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiRowCut with OsiCpxSolverInterface\n" );
    OsiRowCutUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiColCut with OsiCpxSolverInterface\n" );
    OsiColCutUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiCpxSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&cpxSi,mpsDir);
  }
  {
    OsiCpxSolverInterface cpxSi;
    testingMessage( "Testing OsiSimplexInterface with OsiCpxSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&cpxSi,mpsDir);
  }
#endif

#ifdef COIN_USE_SPX
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiRowCut with OsiSpxSolverInterface\n" );
    OsiRowCutUnitTest(&spxSi,mpsDir);
  }
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiColCut with OsiSpxSolverInterface\n" );
    OsiColCutUnitTest(&spxSi,mpsDir);
  }
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiSpxSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&spxSi,mpsDir);
  }
  {
    OsiSpxSolverInterface spxSi;
    testingMessage( "Testing OsiSimplexInterface with OsiSpxSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&spxSi,mpsDir);
  }
#endif

#ifdef COIN_USE_VOL
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing OsiRowCut with OsiVolSolverInterface\n" );
    OsiRowCutUnitTest(&volSi,mpsDir);
  }
  {
    OsiVolSolverInterface volSi;
    testingMessage( "Testing OsiColCut with OsiVolSolverInterface\n" );
    OsiColCutUnitTest(&volSi,mpsDir);
  }
#endif

#ifdef COIN_USE_DYLP
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing OsiRowCut with OsiDylpSolverInterface\n" );
    OsiRowCutUnitTest(&dylpSi,mpsDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing OsiColCut with OsiDylpSolverInterface\n" );
    OsiColCutUnitTest(&dylpSi,mpsDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiDylpSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&dylpSi,mpsDir);
  }
  {
    OsiDylpSolverInterface dylpSi;
    testingMessage( "Testing OsiSimplexInterface with OsiDylpSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&dylpSi,mpsDir);
  }
#endif

#ifdef COIN_USE_GLPK
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiRowCut with OsiGlpkSolverInterface\n" );
    OsiRowCutUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiColCut with OsiGlpkSolverInterface\n" );
    OsiColCutUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiGlpkSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&glpkSi,mpsDir);
  }
  {
    OsiGlpkSolverInterface glpkSi;
    testingMessage( "Testing OsiSimplexInterface with OsiGlpkSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&glpkSi,mpsDir);
  }
#endif

#ifdef COIN_USE_CLP  
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing OsiRowCut with OsiClpSolverInterface\n" );
    OsiRowCutUnitTest(&clpSi,mpsDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing OsiColCut with OsiClpSolverInterface\n" );
    OsiColCutUnitTest(&clpSi,mpsDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiClpSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&clpSi,mpsDir);
  }
  {
    OsiClpSolverInterface clpSi;
    testingMessage( "Testing OsiSimplexInterface with OsiClpSolverInterface\n" );
    OsiSimplexInterfaceCommonUnitTest(&clpSi,mpsDir);
  }
#endif

  testingMessage( "Testing OsiCuts\n" );
  OsiCutsUnitTest();

#ifdef COIN_USE_OSL
  testingMessage( "Testing OsiOslSolverInterface\n" );
  OsiOslSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_USE_XPR
  testingMessage( "Testing OsiXprSolverInterface\n" );
  OsiXprSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_USE_CPX
  testingMessage( "Testing OsiCpxSolverInterface\n" );
  OsiCpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_USE_SPX
  testingMessage( "Testing OsiSpxSolverInterface\n" );
  OsiSpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_USE_VOL
  testingMessage( "Testing OsiVolSolverInterface\n" );
  OsiVolSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_USE_DYLP
  testingMessage( "Testing OsiDylpSolverInterface\n" );
  OsiDylpSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_USE_GLPK
  testingMessage( "Testing OsiGlpkSolverInterface\n" );
  OsiGlpkSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_USE_CLP
  testingMessage( "Testing OsiClpSolverInterface\n" );
  OsiClpSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi;
#   if COIN_USE_OSL
    OsiSolverInterface * oslSi = new OsiOslSolverInterface;
    vecSi.push_back(oslSi);
#endif
#   if COIN_USE_XPR
    OsiSolverInterface * xprSi = new OsiXprSolverInterface;
    vecSi.push_back(xprSi);
#endif
#   if COIN_USE_CPX
    OsiSolverInterface * cpxSi = new OsiCpxSolverInterface;
    vecSi.push_back(cpxSi);
#endif
#   if COIN_USE_SPX
    OsiSolverInterface * spxSi = new OsiSpxSolverInterface;
    vecSi.push_back(spxSi);
#endif
#   if COIN_USE_CLP
    OsiSolverInterface * clpSi = new OsiClpSolverInterface;
    // Okay this is where John Forrest cheats by giving hints
    clpSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
    vecSi.push_back(clpSi);
#endif
#   if COIN_USE_DYLP
    OsiSolverInterface * dylpSi = new OsiDylpSolverInterface;
    vecSi.push_back(dylpSi);
#endif
#   if COIN_USE_GLPK
    OsiSolverInterface * glpkSi = new OsiGlpkSolverInterface;
    vecSi.push_back(glpkSi);
#endif
#   if COIN_USE_VOL
    OsiSolverInterface * volSi = new OsiVolSolverInterface;
    vecSi.push_back(volSi);
#endif

    testingMessage( "Testing OsiSolverInterface\n" );
    OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

    unsigned int i;
    for (i=0; i<vecSi.size(); i++)
      delete vecSi[i];
  }
  else {
    testingMessage( "***Skipped Testing of OsiSolverInterface    ***\n" );
    testingMessage( "***use -testOsiSolverInterface to test class***\n" );
  }

  testingMessage( "All tests completed successfully\n" );
  return 0;
}

 
// Display message on stdout and stderr
void testingMessage( const char * const msg )
{
  std::cerr <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}

