// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// Test individual classes or groups of classes

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>

#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "CoinSort.hpp"
#include "OsiShallowPackedVector.hpp"
#include "OsiPackedVector.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiRowCutDebugger.hpp"
#ifdef COIN_TEST_OSI_READER
#include "OsiMpsReader.hpp"
#endif
#ifdef COIN_USE_OSL
#include "OsiOslSolverInterface.hpp"
#endif
#ifdef COIN_USE_XPR
#include "OsiXprSolverInterface.hpp"
//#include <xpresso.h>
#endif
#ifdef COIN_USE_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#ifdef COIN_USE_VOL
#include "OsiVolSolverInterface.hpp"
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

  // *FIXME* : these tests should be written... 
  //  testingMessage( "Testing CoinHelperFunctions\n" );
  //  CoinHelperFunctionsUnitTest();
  //  testingMessage( "Testing CoinSort\n" );
  //  tripleCompareUnitTest();
  //  testingMessage( "Testing CoinError\n" );
  //  CoinErrorUnitTest();

  testingMessage( "Testing OsiShallowPackedVector\n" );
  OsiShallowPackedVectorUnitTest();

  testingMessage( "Testing OsiPackedVector\n" );
  OsiPackedVectorUnitTest();

  //testingMessage( "Testing OsiPackedMatrix\n" );
  //OsiPackedMatrixUnitTest();

#ifdef COIN_TEST_OSI_READER
  testingMessage( "Testing OsiMpsReader\n" );
  OsiMpsReaderUnitTest(mpsDir);
#endif

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

  testingMessage( "Testing OsiCuts\n" );
  OsiCutsUnitTest();

#ifdef COIN_USE_OSL
  testingMessage( "Testing OsiOslSolverInterface\n" );
  OsiOslSolverInterfaceUnitTest(mpsDir);
#endif

#ifdef COIN_USE_XPR
  testingMessage( "Testing OsiXprSolverInterface\n" );
  OsiXprSolverInterfaceUnitTest(mpsDir);
#endif

#ifdef COIN_USE_CPX
  testingMessage( "Testing OsiCpxSolverInterface\n" );
  OsiCpxSolverInterfaceUnitTest(mpsDir);
#endif

#ifdef COIN_USE_VOL
  testingMessage( "Testing OsiVolSolverInterface\n" );
  OsiVolSolverInterfaceUnitTest(mpsDir);
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

