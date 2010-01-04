// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// Test individual classes or groups of classes

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiConfig.h"

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <iostream>
#include <cstdio>

#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiRowCutDebugger.hpp"

/*
  Some convenient undef's, to make it easy to isolate a particular solver.
  Uncomment to disable a solver that's included in the build. Leave them
  commented if you're happy with running the unitTest for all solvers in the
  build.
*/

// #undef COIN_HAS_OSL
// #undef COIN_HAS_XPR
// #undef COIN_HAS_CPX
// #undef COIN_HAS_SPX
// #undef COIN_HAS_VOL
// #undef COIN_HAS_DYLP
// #undef COIN_HAS_GLPK
// #undef COIN_HAS_FMP
// #undef COIN_HAS_CLP
// #undef COIN_HAS_SYMPHONY
// #undef COIN_HAS_MSK
// #undef COIN_HAS_CBC
// #undef COIN_HAS_GRB

#ifdef COIN_HAS_OSL
#include "OsiOslSolverInterface.hpp"
#endif
#ifdef COIN_HAS_XPR
#include "OsiXprSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CPX
#include "OsiCpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SPX
#include "OsiSpxSolverInterface.hpp"
#endif
#ifdef COIN_HAS_VOL
#include "OsiVolSolverInterface.hpp"
#endif
#ifdef COIN_HAS_DYLP
#include "OsiDylpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_GLPK
#include "OsiGlpkSolverInterface.hpp"
#endif
#ifdef COIN_HAS_FMP
#include "OsiFmpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_SYMPHONY
#include "OsiSymSolverInterface.hpp"
#endif
#ifdef COIN_HAS_MSK
#include "OsiMskSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CBC
#include "OsiCbcSolverInterface.hpp"
#endif
#ifdef COIN_HAS_GRB
#include "OsiGrbSolverInterface.hpp"
#endif

namespace {

/*
  If anyone is feeling ambitious, it'd be a really good idea to handle i/o for
  the unittest by way of a standard CoinMessageHandler. Might require a bit of
  tweaking in CoinMessageHandler.
*/
 
// Display message on stdout and stderr. Flush cout buffer before printing the
// message, so that output comes out in order in spite of buffered cout.

void testingMessage( const char * const msg )
{
  std::cout.flush() ;
  std::cerr <<msg;
  //cout <<endl <<"*****************************************"
  //     <<endl <<msg <<endl;
}


/*
  Utility routine to process command line parameters. An unrecognised parameter
  will trigger the help message and a return value of false.
  
  This should be replaced with the one of the standard CoinUtils parameter
  mechanisms.
*/
bool processParameters (int argc, const char **argv,
			std::map<std::string,std::string> &parms)

{ 
/*
  Initialise the parameter keywords.
*/
  std::set<std::string> definedKeyWords;
  definedKeyWords.insert("-cerr2cout");
  definedKeyWords.insert("-mpsDir");
  definedKeyWords.insert("-netlibDir");
  definedKeyWords.insert("-testOsiSolverInterface");
  definedKeyWords.insert("-nobuf");
  definedKeyWords.insert("-cutsOnly");
/*
  Set default values for data directories.
*/
  const char dirsep =  CoinFindDirSeparator() ;
  std::string pathTmp ;

  pathTmp = ".." ;
  pathTmp += dirsep ;
  pathTmp += ".." ;
  pathTmp += dirsep ;
  pathTmp += "Data" ;
  pathTmp += dirsep ;

  parms["-mpsDir"] = pathTmp + "Sample"  ;
  parms["-netlibDir"] = pathTmp + "Netlib" ;

/*
  Read the command line parameters and fill a map of parameter keys and
  associated data. The parser allows for parameters which are only a keyword,
  or parameters of the form keyword=value (no spaces).
*/
  for (int i = 1 ; i < argc ; i++)
  { std::string parm(argv[i]) ;
    std::string key,value ;
    std::string::size_type eqPos = parm.find('=');

    if (eqPos == std::string::npos)
    { key = parm ; }
    else
    { key = parm.substr(0,eqPos) ;
      value = parm.substr(eqPos+1) ; }
/*
  Is the specifed key valid?
*/
    if (definedKeyWords.find(key) == definedKeyWords.end())
    { std::cerr << "Undefined parameter \"" << key << "\"." << std::endl ;
      std::cerr
	<< "Usage: "
	<< "unitTest [-nobuf] [-mpsDir=V1] [-netlibDir=V2] "
        << "[-testOsiSolverInterface] [-cutsOnly]" << std::endl ;
      std::cerr << "  where:" << std::endl ;
      std::cerr
	<< "    "
	<< "-cerr2cout: redirect cerr to cout; sometimes useful." << std::endl
	<< "\t" << "to synchronise cout & cerr." << std::endl ;
      std::cerr
	<< "    "
	<< "-mpsDir: directory containing mps test files." << std::endl
        << "\t" << "Default value V1=\"../../Data/Sample\"" << std::endl ;
      std::cerr
	<< "    "
	<< "-netlibDir: directory containing netlib files." << std::endl
        << "\t" << "Default value V2=\"../../Data/Netlib\"" << std::endl ;
      std::cerr
	<< "    "
	<< "-testOsiSolverInterface: "
        << "run each OSI on the netlib problem set." << std::endl
	<< "\t"
	<< "Default is to not run the netlib problem set." << std::endl ;
      std::cerr
	<< "    "
	<< "-cutsOnly: If specified, only OsiCut tests are run." << std::endl ;
      std::cerr
	<< "    "
        << "-nobuf: use unbuffered output." << std::endl
	<< "\t" << "Default is buffered output." << std::endl ;
      
      return (false) ; }
/*
  Valid keyword; stash the value for later reference.
*/
    parms[key]=value ; }
/*
  Tack the directory separator onto the data directories so we don't have to
  worry about it later.
*/
  parms["-mpsDir"] += dirsep ;
  parms["-netlibDir"] += dirsep ;
/*
  Did the user request unbuffered i/o? It seems we need to go after this
  through stdio --- using pubsetbuf(0,0) on the C++ streams has no
  discernible affect. Nor, for that matter, did setting the unitbuf flag on
  the streams. Why? At a guess, sync_with_stdio connects the streams to the
  stdio buffers, and the C++ side isn't programmed to change them?
*/
  if (parms.find("-nobuf") != parms.end())
  { // std::streambuf *coutBuf, *cerrBuf ;
    // coutBuf = std::cout.rdbuf() ;
    // coutBuf->pubsetbuf(0,0) ;
    // cerrBuf = std::cerr.rdbuf() ;
    // cerrBuf->pubsetbuf(0,0) ;
    setbuf(stderr,0) ;
    setbuf(stdout,0) ; }
/*
  Did the user request a redirect for cerr? This must occur before any i/o is
  performed.
*/
  if (parms.find("-cerr2cout") != parms.end())
  { std::cerr.rdbuf(std::cout.rdbuf()) ; }

  return (true) ; }


}	// end file-local namespace



//----------------------------------------------------------------
// unitTest [-nobuf] [-mpsDir=V1] [-netlibDir=V2] [-testOsiSolverInterface]
//	    [-cutsOnly]
// 
// where:
//   -nobuf: remove buffering on cout (stdout); useful to keep cout and cerr
//	 messages synchronised when redirecting output to a file or pipe.
//   -mpsDir: directory containing mps test files
//       Default value V1="../../Data/Sample"    
//   -netlibDir: directory containing netlib files
//       Default value V2="../../Data/Netlib"
//   -testOsiSolverInterface
//       If specified, then OsiSolveInterface::unitTest
//       is skipped over and not run.
//   -cutsOnly
//	 If specified, only OsiCut tests are run.
//
// All parameters are optional.
//----------------------------------------------------------------

int main (int argc, const char *argv[])

{ int totalErrCnt = 0;

/*
  Start off with various bits of initialisation that don't really belong
  anywhere else.

  First off, synchronise C++ stream i/o with C stdio. This makes debugging
  output a bit more comprehensible. It still suffers from interleave of cout
  (stdout) and cerr (stderr), but -nobuf deals with that.
*/
  std::ios::sync_with_stdio() ;
/*
  Suppress an popup window that Windows shows in response to a crash. See
  note at head of file.
*/
  WindowsErrorPopupBlocker();
/*
  Might as well make use of this convenient Xpress feature.
*/
#ifdef COIN_HAS_XPR
  OsiXprSolverInterface::setLogFileName("xprCallTrace.txt");
#endif

/*
  Process command line parameters.
*/
  std::map<std::string,std::string> parms ;

  if (processParameters(argc,argv,parms) == false)
  { return (1) ; }

  std::string mpsDir = parms["-mpsDir"] ;
  std::string netlibDir = parms["-netlibDir"] ;

try {
/*
  Test Osi{Row,Col}Cut routines.
*/
#ifdef COIN_HAS_OSL  
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

#ifdef COIN_HAS_XPR  
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

#ifdef COIN_HAS_CPX
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

#ifdef COIN_HAS_SPX
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
#endif

#ifdef COIN_HAS_VOL
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

#ifdef COIN_HAS_DYLP
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
#endif
  

#ifdef COIN_HAS_FMP
  {
    OsiFmpSolverInterface fmpSi;
    testingMessage( "Testing OsiRowCut with OsiFmpSolverInterface\n" );
    OsiRowCutUnitTest(&fmpSi,mpsDir);
  }
  {
    OsiFmpSolverInterface fmpSi;
    testingMessage( "Testing OsiColCut with OsiFmpSolverInterface\n" );
    OsiColCutUnitTest(&fmpSi,mpsDir);
  }
  // FortMP does not presently pass this test
  {
    OsiFmpSolverInterface fmpSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiFmpSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&fmpSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_GLPK
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
#endif

#ifdef COIN_HAS_CLP  
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
#endif

#ifdef COIN_HAS_SYMPHONY
  {
    OsiSymSolverInterface symSi;
    testingMessage( "Testing OsiRowCut with OsiSymSolverInterface\n" );
    OsiRowCutUnitTest(&symSi,mpsDir);
  }
  {
    OsiSymSolverInterface symSi;
    testingMessage( "Testing OsiColCut with OsiSymSolverInterface\n" );
    OsiColCutUnitTest(&symSi,mpsDir);
  }
  {
    OsiSymSolverInterface symSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiSymSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&symSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_MSK  
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiRowCut with OsiMskSolverInterface\n" );
    OsiRowCutUnitTest(&MskSi,mpsDir);
  }
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiColCut with OsiMskSolverInterface\n" );
    OsiColCutUnitTest(&MskSi,mpsDir);
  }
  {
    OsiMskSolverInterface MskSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiMskSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&MskSi,mpsDir);
  }
#endif
#ifdef COIN_HAS_CBC
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCut with OsiCbcSolverInterface\n" );
    OsiRowCutUnitTest(&cbcSi,mpsDir);
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiColCut with OsiCbcSolverInterface\n" );
    OsiColCutUnitTest(&cbcSi,mpsDir);
  }
  {
    OsiCbcSolverInterface cbcSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiCbcSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&cbcSi,mpsDir);
  }
#endif

#ifdef COIN_HAS_GRB
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiRowCut with OsiGrbSolverInterface\n" );
    OsiRowCutUnitTest(&grbSi,mpsDir);
  }
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiColCut with OsiGrbSolverInterface\n" );
    OsiColCutUnitTest(&grbSi,mpsDir);
  }
  {
    OsiGrbSolverInterface grbSi;
    testingMessage( "Testing OsiRowCutDebugger with OsiGrbSolverInterface\n" );
    OsiRowCutDebuggerUnitTest(&grbSi,mpsDir);
  }
#endif

  testingMessage( "Testing OsiCuts\n" );
  OsiCutsUnitTest();

/*
  Testing OsiCuts only? A useful option when doing memory access and leak
  checks. Keeps the run time to something reasonable.
*/
  if (parms.find("-cutsOnly") != parms.end()) {
    testingMessage( "Stopping after OsiCuts tests.\n" );
    return 0;
  }

/*
  Run the OsiXXX class test for each solver. It's up to the solver implementor
  to decide whether or not to run OsiSolverInterfaceCommonUnitTest. Arguably
  this should be required.
*/
#ifdef COIN_HAS_OSL
  testingMessage( "Testing OsiOslSolverInterface\n" );
  OsiOslSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_XPR
  testingMessage( "Testing OsiXprSolverInterface\n" );
  OsiXprSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_CPX
  testingMessage( "Testing OsiCpxSolverInterface\n" );
  OsiCpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_SPX
  testingMessage( "Testing OsiSpxSolverInterface\n" );
  OsiSpxSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_VOL
  testingMessage( "Testing OsiVolSolverInterface\n" );
  totalErrCnt += OsiVolSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_DYLP
  testingMessage( "Testing OsiDylpSolverInterface\n" );
  int errCnt = OsiDylpSolverInterfaceUnitTest(mpsDir,netlibDir);
  if (errCnt)
  { std::cerr
      << "OsiDylpSolverInterface testing issue: "
      << errCnt << " errors reported by OsiDylp unit test." << std::endl ;
    totalErrCnt += errCnt ; }
#endif
  
#ifdef COIN_HAS_GLPK
  testingMessage( "Testing OsiGlpkSolverInterface\n" );
  totalErrCnt += OsiGlpkSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_HAS_FMP
  testingMessage( "Testing OsiFmpSolverInterface\n" );
  OsiFmpSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif
  
#ifdef COIN_HAS_CLP
  testingMessage( "Testing OsiClpSolverInterface\n" );
  totalErrCnt += OsiClpSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_MSK
  testingMessage( "Testing OsiMskSolverInterface\n" );
  OsiMskSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_CBC
  testingMessage( "Testing OsiCbcSolverInterface\n" );
  totalErrCnt += OsiCbcSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_SYMPHONY
  testingMessage( "Testing OsiSymSolverInterface\n" );
  totalErrCnt += OsiSymSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

#ifdef COIN_HAS_GRB
  testingMessage( "Testing OsiGrbSolverInterface\n" );
  OsiGrbSolverInterfaceUnitTest(mpsDir,netlibDir);
#endif

/*
  Each solver has run its specialised unit test. Check now to see if we need to
  run through the Netlib problems.
*/
  if (parms.find("-testOsiSolverInterface") != parms.end())
  {
    // Create vector of solver interfaces
    std::vector<OsiSolverInterface*> vecSi;
#   if COIN_HAS_OSL
    OsiSolverInterface * oslSi = new OsiOslSolverInterface;
    vecSi.push_back(oslSi);
#endif
#   if COIN_HAS_XPR
    OsiSolverInterface * xprSi = new OsiXprSolverInterface;
    vecSi.push_back(xprSi);
#endif
#   if COIN_HAS_CPX
    OsiSolverInterface * cpxSi = new OsiCpxSolverInterface;
    vecSi.push_back(cpxSi);
#endif
#   if COIN_HAS_SPX
    OsiSolverInterface * spxSi = new OsiSpxSolverInterface;
    vecSi.push_back(spxSi);
#endif
#   if COIN_HAS_CLP
    OsiSolverInterface * clpSi = new OsiClpSolverInterface;
    // Okay this is where John Forrest cheats by giving hints
    clpSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
    clpSi->setHintParam(OsiDoReducePrint,true,OsiHintDo);
    vecSi.push_back(clpSi);
#endif
#   if COIN_HAS_SYMPHONY
    OsiSolverInterface * symSi = new OsiSymSolverInterface;
    { OsiSymSolverInterface *reallySymSi =
	  dynamic_cast<OsiSymSolverInterface *>(symSi) ;
      reallySymSi->setSymParam(OsiSymVerbosity, -2) ; }
    vecSi.push_back(symSi);
#endif
#   if COIN_HAS_DYLP
    OsiSolverInterface * dylpSi = new OsiDylpSolverInterface;
    // Heh, if it's good enough for John ...
    dylpSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry) ;
    dylpSi->setHintParam(OsiDoReducePrint,true,OsiHintDo) ;
    vecSi.push_back(dylpSi);
#endif
#   if COIN_HAS_GLPK
    OsiSolverInterface * glpkSi = new OsiGlpkSolverInterface;
    glpkSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry) ;
    glpkSi->setHintParam(OsiDoReducePrint,true,OsiHintDo) ;
    vecSi.push_back(glpkSi);
#endif
#   if COIN_HAS_FMP
    OsiSolverInterface * fmpSi = new OsiFmpSolverInterface;
    vecSi.push_back(fmpSi);
#endif
#   if COIN_HAS_MSK
    OsiSolverInterface * MskSi = new OsiMskSolverInterface;
    vecSi.push_back(MskSi);
#endif
#   if COIN_HAS_CBC
    OsiSolverInterface * cbcSi = new OsiCbcSolverInterface;
    // Okay this is where John Forrest cheats by giving hints
    cbcSi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
    cbcSi->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    vecSi.push_back(cbcSi);
#endif
#   if COIN_HAS_VOL
    OsiSolverInterface * volSi = new OsiVolSolverInterface;
    vecSi.push_back(volSi);
#endif
#   if COIN_HAS_GRB
    OsiSolverInterface * grbSi = new OsiGrbSolverInterface;
    vecSi.push_back(grbSi);
#endif

    testingMessage( "Testing OsiSolverInterface on Netlib problems.\n" );
    OsiSolverInterfaceMpsUnitTest(vecSi,netlibDir);

    unsigned int i;
    for (i=0; i<vecSi.size(); i++)
      delete vecSi[i];
  }
  else {
    testingMessage( "***Skipped Testing of OsiSolverInterface on Netlib problems***\n" );
    testingMessage( "***use -testOsiSolverInterface to run them.***\n" );
  }
} catch (CoinError& error) {
  std::cout.flush();
  std::cerr << "Caught CoinError exception: ";
  error.print(true);
  totalErrCnt += 1;
}
/*
  We're done. Report on the results.
*/
  if (totalErrCnt)
  { std::cout.flush() ;
    std::cerr
      << "Tests completed with " << totalErrCnt << " errors." << std::endl ; 
  } else
  { testingMessage("All tests completed successfully\n") ; }
  return totalErrCnt;
}

