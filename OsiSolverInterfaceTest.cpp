//  LAST EDIT: Tue Aug 28 17:40:30 2001 by Tobias Pfender (opt14!bzfpfend) 
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
 
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>
#include <cstdio>

#include "OsiSolverInterface.hpp"
#include "OsiVolSolverInterface.hpp"
#include "OsiFloatEqual.hpp"
#include "OsiPackedVector.hpp"
#include "OsiPackedMatrix.hpp"
#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

//--------------------------------------------------------------------------
// A helper function to compare the equivalence of two vectors 
static bool
equivalentVectors(const OsiSolverInterface * si1,
		  const OsiSolverInterface * si2,
		  double tol,
		  const double * v1,
		  const double * v2,
		  int size)
{
  bool retVal = true;
  OsiRelFltEq eq(tol);
  int i;
  for ( i=0; i<size; i++ ) {
    
    // If both are equal to infinity then iterate
    if ( fabs(v1[i])==si1->getInfinity() && fabs(v2[i])==si2->getInfinity() )
       continue;
    
    // Test to see if equal
    if ( !eq(v1[i],v2[i]) ) {
      std::cerr <<"eq " <<i <<" " <<v1[i] <<" " <<v2[i] <<std::endl;
      //const OsiOslSolverInterface * oslSi = dynamic_cast<const OsiOslSolverInterface*>(si1);
      //assert(oslSi != NULL);
      //EKKModel * m = ((OsiOslSolverInterface*)oslSi)->modelPtr();
      //ekk_printSolution(m);
      retVal = false;
      break;
    }

  }
  return retVal;
}

//--------------------------------------------------------------------------

void
OsiSolverInterfaceMpsUnitTest(
  const std::vector<OsiSolverInterface*> & vecEmptySiP,
  const std::string & mpsDir)
{
  int i;
  unsigned int m;

  // Establish that the solver interfaces obtain good solutions
  // for a set of problems (netlib)

  // Define test problems: 
  //   mps names, 
  //   maximization or minimization, 
  //   Number of rows and columns in problem, and
  //   objective function value
  std::vector<std::string> mpsName;
  std::vector<bool> min;
  std::vector<int> nRows;
  std::vector<int> nCols;
  std::vector<double> objValue;
  std::vector<double> objValueTol;
  
  mpsName.push_back("25fv47");
  min.push_back(true);
  nRows.push_back(822);
  nCols.push_back(1571);
  objValueTol.push_back(1.E-10);
  objValue.push_back(5.5018458883E+03);
  
  mpsName.push_back("80bau3b");min.push_back(true);nRows.push_back(2263);nCols.push_back(9799);objValueTol.push_back(1.e-10);objValue.push_back(9.8722419241E+05);
  mpsName.push_back("adlittle");min.push_back(true);nRows.push_back(57);nCols.push_back(97);objValueTol.push_back(1.e-10);objValue.push_back(2.2549496316e+05);
  mpsName.push_back("afiro");min.push_back(true);nRows.push_back(28);nCols.push_back(32);objValueTol.push_back(1.e-10);objValue.push_back(-4.6475314286e+02);
  mpsName.push_back("agg");min.push_back(true);nRows.push_back(489);nCols.push_back(163);objValueTol.push_back(1.e-10);objValue.push_back(-3.5991767287e+07);
  mpsName.push_back("agg2");min.push_back(true);nRows.push_back(517);nCols.push_back(302);objValueTol.push_back(1.e-10);objValue.push_back(-2.0239252356e+07);
  mpsName.push_back("agg3");min.push_back(true);nRows.push_back(517);nCols.push_back(302);objValueTol.push_back(1.e-10);objValue.push_back(1.0312115935e+07);
  mpsName.push_back("bandm");min.push_back(true);nRows.push_back(306);nCols.push_back(472);objValueTol.push_back(1.e-10);objValue.push_back(-1.5862801845e+02);
  mpsName.push_back("beaconfd");min.push_back(true);nRows.push_back(174);nCols.push_back(262);objValueTol.push_back(1.e-10);objValue.push_back(3.3592485807e+04);
  // CPLEX READ ERROR  mpsName.push_back("blend");min.push_back(true);nRows.push_back(75);nCols.push_back(83);objValueTol.push_back(1.e-10);objValue.push_back(-3.0812149846e+01);
  mpsName.push_back("bnl1");min.push_back(true);nRows.push_back(644);nCols.push_back(1175);objValueTol.push_back(1.e-10);objValue.push_back(1.9776295615E+03);
  mpsName.push_back("bnl2");min.push_back(true);nRows.push_back(2325);nCols.push_back(3489);objValueTol.push_back(1.e-10);objValue.push_back(1.8112365404e+03);
  mpsName.push_back("boeing1");min.push_back(true);nRows.push_back(/*351*/352);nCols.push_back(384);objValueTol.push_back(1.e-10);objValue.push_back(-3.3521356751e+02);
  mpsName.push_back("boeing2");min.push_back(true);nRows.push_back(167);nCols.push_back(143);objValueTol.push_back(1.e-10);objValue.push_back(-3.1501872802e+02);
  mpsName.push_back("bore3d");min.push_back(true);nRows.push_back(234);nCols.push_back(315);objValueTol.push_back(1.e-10);objValue.push_back(1.3730803942e+03);
  mpsName.push_back("brandy");min.push_back(true);nRows.push_back(221);nCols.push_back(249);objValueTol.push_back(1.e-10);objValue.push_back(1.5185098965e+03);
  mpsName.push_back("capri");min.push_back(true);nRows.push_back(272);nCols.push_back(353);objValueTol.push_back(1.e-10);objValue.push_back(2.6900129138e+03);
  mpsName.push_back("cycle");min.push_back(true);nRows.push_back(1904);nCols.push_back(2857);objValueTol.push_back(1.e-9);objValue.push_back(-5.2263930249e+00);
  mpsName.push_back("czprob");min.push_back(true);nRows.push_back(930);nCols.push_back(3523);objValueTol.push_back(1.e-10);objValue.push_back(2.1851966989e+06);
  mpsName.push_back("d2q06c");min.push_back(true);nRows.push_back(2172);nCols.push_back(5167);objValueTol.push_back(1.e-7);objValue.push_back(122784.21557456);
  mpsName.push_back("d6cube");min.push_back(true);nRows.push_back(416);nCols.push_back(6184);objValueTol.push_back(1.e-10);objValue.push_back(3.1549166667e+02);
  mpsName.push_back("degen2");min.push_back(true);nRows.push_back(445);nCols.push_back(534);objValueTol.push_back(1.e-10);objValue.push_back(-1.4351780000e+03);
  mpsName.push_back("degen3");min.push_back(true);nRows.push_back(1504);nCols.push_back(1818);objValueTol.push_back(1.e-10);objValue.push_back(-9.8729400000e+02);
  // CPLEX READ ERROR  mpsName.push_back("dfl001");min.push_back(true);nRows.push_back(6072);nCols.push_back(12230);objValueTol.push_back(1.e-5);objValue.push_back(1.1266396047E+07);
  mpsName.push_back("e226");min.push_back(true);nRows.push_back(224);nCols.push_back(282);objValueTol.push_back(1.e-10);objValue.push_back(/*-1.8751929066e+01*/-11.638929066371);
  mpsName.push_back("etamacro");min.push_back(true);nRows.push_back(401);nCols.push_back(688);objValueTol.push_back(1.e-6);objValue.push_back(-7.5571521774e+02 );
  mpsName.push_back("fffff800");min.push_back(true);nRows.push_back(525);nCols.push_back(854);objValueTol.push_back(1.e-6);objValue.push_back(5.5567961165e+05);
  mpsName.push_back("finnis");min.push_back(true);nRows.push_back(498);nCols.push_back(614);objValueTol.push_back(1.e-6);objValue.push_back(1.7279096547e+05);
  mpsName.push_back("fit1d");min.push_back(true);nRows.push_back(25);nCols.push_back(1026);objValueTol.push_back(1.e-10);objValue.push_back(-9.1463780924e+03);
  mpsName.push_back("fit1p");min.push_back(true);nRows.push_back(628);nCols.push_back(1677);objValueTol.push_back(1.e-10);objValue.push_back(9.1463780924e+03);
  mpsName.push_back("fit2d");min.push_back(true);nRows.push_back(26);nCols.push_back(10500);objValueTol.push_back(1.e-10);objValue.push_back(-6.8464293294e+04);
  mpsName.push_back("fit2p");min.push_back(true);nRows.push_back(3001);nCols.push_back(13525);objValueTol.push_back(1.e-9);objValue.push_back(6.8464293232e+04);
  // CPLEX READ ERROR  mpsName.push_back("forplan");min.push_back(true);nRows.push_back(162);nCols.push_back(421);objValueTol.push_back(1.e-6);objValue.push_back(-6.6421873953e+02);
  mpsName.push_back("ganges");min.push_back(true);nRows.push_back(1310);nCols.push_back(1681);objValueTol.push_back(1.e-5);objValue.push_back(-1.0958636356e+05);
  // CPLEX READ ERROR  mpsName.push_back("gfrd-pnc");min.push_back(true);nRows.push_back(617);nCols.push_back(1092);objValueTol.push_back(1.e-10);objValue.push_back(6.9022359995e+06);
  mpsName.push_back("greenbea");min.push_back(true);nRows.push_back(2393);nCols.push_back(5405);objValueTol.push_back(1.e-10);objValue.push_back(/*-7.2462405908e+07*/-72555248.129846);
  mpsName.push_back("greenbeb");min.push_back(true);nRows.push_back(2393);nCols.push_back(5405);objValueTol.push_back(1.e-10);objValue.push_back(/*-4.3021476065e+06*/-4302260.2612066);
  mpsName.push_back("grow15");min.push_back(true);nRows.push_back(301);nCols.push_back(645);objValueTol.push_back(1.e-10);objValue.push_back(-1.0687094129e+08);
  mpsName.push_back("grow22");min.push_back(true);nRows.push_back(441);nCols.push_back(946);objValueTol.push_back(1.e-10);objValue.push_back(-1.6083433648e+08);
  mpsName.push_back("grow7");min.push_back(true);nRows.push_back(141);nCols.push_back(301);objValueTol.push_back(1.e-10);objValue.push_back(-4.7787811815e+07);
  mpsName.push_back("israel");min.push_back(true);nRows.push_back(175);nCols.push_back(142);objValueTol.push_back(1.e-10);objValue.push_back(-8.9664482186e+05);
  mpsName.push_back("kb2");min.push_back(true);nRows.push_back(44);nCols.push_back(41);objValueTol.push_back(1.e-10);objValue.push_back(-1.7499001299e+03);
  mpsName.push_back("lotfi");min.push_back(true);nRows.push_back(154);nCols.push_back(308);objValueTol.push_back(1.e-10);objValue.push_back(-2.5264706062e+01);
  mpsName.push_back("maros");min.push_back(true);nRows.push_back(847);nCols.push_back(1443);objValueTol.push_back(1.e-10);objValue.push_back(-5.8063743701e+04);
  mpsName.push_back("maros-r7");min.push_back(true);nRows.push_back(3137);nCols.push_back(9408);objValueTol.push_back(1.e-10);objValue.push_back(1.4971851665e+06);
  mpsName.push_back("modszk1");min.push_back(true);nRows.push_back(688);nCols.push_back(1620);objValueTol.push_back(1.e-10);objValue.push_back(3.2061972906e+02);
  mpsName.push_back("nesm");min.push_back(true);nRows.push_back(663);nCols.push_back(2923);objValueTol.push_back(1.e-5);objValue.push_back(1.4076073035e+07);
  mpsName.push_back("perold");min.push_back(true);nRows.push_back(626);nCols.push_back(1376);objValueTol.push_back(1.e-6);objValue.push_back(-9.3807580773e+03);
  mpsName.push_back("pilot");min.push_back(true);nRows.push_back(1442);nCols.push_back(3652);objValueTol.push_back(1.e-5);objValue.push_back(/*-5.5740430007e+02*/-557.48972927292);
  mpsName.push_back("pilot4");min.push_back(true);nRows.push_back(411);nCols.push_back(1000);objValueTol.push_back(1.e-8);objValue.push_back(-2.5811392641e+03);
  mpsName.push_back("pilot87");min.push_back(true);nRows.push_back(2031);nCols.push_back(4883);objValueTol.push_back(1.e-5);objValue.push_back(3.0171072827e+02);
  mpsName.push_back("pilotnov");min.push_back(true);nRows.push_back(976);nCols.push_back(2172);objValueTol.push_back(1.e-10);objValue.push_back(-4.4972761882e+03);
  //mpsName.push_back("qap8");min.push_back(true);nRows.push_back(913);nCols.push_back(1632);objValueTol.push_back(1.e-10);objValue.push_back(2.0350000000e+02);
  //mpsName.push_back("qap12");min.push_back(true);nRows.push_back(3193);nCols.push_back(8856);objValueTol.push_back(1.e-10);objValue.push_back(5.2289435056e+02);
  //mpsName.push_back("qap15");min.push_back(true);nRows.push_back(6331);nCols.push_back(22275);objValueTol.push_back(1.e-10);objValue.push_back(1.0409940410e+03);
  mpsName.push_back("recipe");min.push_back(true);nRows.push_back(92);nCols.push_back(180);objValueTol.push_back(1.e-10);objValue.push_back(-2.6661600000e+02);
  mpsName.push_back("sc105");min.push_back(true);nRows.push_back(106);nCols.push_back(103);objValueTol.push_back(1.e-10);objValue.push_back(-5.2202061212e+01);
  mpsName.push_back("sc205");min.push_back(true);nRows.push_back(206);nCols.push_back(203);objValueTol.push_back(1.e-10);objValue.push_back(-5.2202061212e+01);
  mpsName.push_back("sc50a");min.push_back(true);nRows.push_back(51);nCols.push_back(48);objValueTol.push_back(1.e-10);objValue.push_back(-6.4575077059e+01);
  mpsName.push_back("sc50b");min.push_back(true);nRows.push_back(51);nCols.push_back(48);objValueTol.push_back(1.e-10);objValue.push_back(-7.0000000000e+01);
  mpsName.push_back("scagr25");min.push_back(true);nRows.push_back(472);nCols.push_back(500);objValueTol.push_back(1.e-10);objValue.push_back(-1.4753433061e+07);
  mpsName.push_back("scagr7");min.push_back(true);nRows.push_back(130);nCols.push_back(140);objValueTol.push_back(1.e-6);objValue.push_back(-2.3313892548e+06);
  mpsName.push_back("scfxm1");min.push_back(true);nRows.push_back(331);nCols.push_back(457);objValueTol.push_back(1.e-10);objValue.push_back(1.8416759028e+04);
  mpsName.push_back("scfxm2");min.push_back(true);nRows.push_back(661);nCols.push_back(914);objValueTol.push_back(1.e-10);objValue.push_back(3.6660261565e+04);
  mpsName.push_back("scfxm3");min.push_back(true);nRows.push_back(991);nCols.push_back(1371);objValueTol.push_back(1.e-10);objValue.push_back(5.4901254550e+04);
  mpsName.push_back("scorpion");min.push_back(true);nRows.push_back(389);nCols.push_back(358);objValueTol.push_back(1.e-10);objValue.push_back(1.8781248227e+03);
  mpsName.push_back("scrs8");min.push_back(true);nRows.push_back(491);nCols.push_back(1169);objValueTol.push_back(1.e-5);objValue.push_back(9.0429998619e+02);
  mpsName.push_back("scsd1");min.push_back(true);nRows.push_back(78);nCols.push_back(760);objValueTol.push_back(1.e-10);objValue.push_back(8.6666666743e+00);
  mpsName.push_back("scsd6");min.push_back(true);nRows.push_back(148);nCols.push_back(1350);objValueTol.push_back(1.e-10);objValue.push_back(5.0500000078e+01);
  mpsName.push_back("scsd8");min.push_back(true);nRows.push_back(398);nCols.push_back(2750);objValueTol.push_back(1.e-10);objValue.push_back(9.0499999993e+02);
  mpsName.push_back("sctap1");min.push_back(true);nRows.push_back(301);nCols.push_back(480);objValueTol.push_back(1.e-10);objValue.push_back(1.4122500000e+03);
  mpsName.push_back("sctap2");min.push_back(true);nRows.push_back(1091);nCols.push_back(1880);objValueTol.push_back(1.e-10);objValue.push_back(1.7248071429e+03);
  mpsName.push_back("sctap3");min.push_back(true);nRows.push_back(1481);nCols.push_back(2480);objValueTol.push_back(1.e-10);objValue.push_back(1.4240000000e+03);
  mpsName.push_back("seba");min.push_back(true);nRows.push_back(516);nCols.push_back(1028);objValueTol.push_back(1.e-10);objValue.push_back(1.5711600000e+04);
  mpsName.push_back("share1b");min.push_back(true);nRows.push_back(118);nCols.push_back(225);objValueTol.push_back(1.e-10);objValue.push_back(-7.6589318579e+04);
  mpsName.push_back("share2b");min.push_back(true);nRows.push_back(97);nCols.push_back(79);objValueTol.push_back(1.e-10);objValue.push_back(-4.1573224074e+02);
  mpsName.push_back("shell");min.push_back(true);nRows.push_back(537);nCols.push_back(1775);objValueTol.push_back(1.e-10);objValue.push_back(1.2088253460e+09);
  mpsName.push_back("ship04l");min.push_back(true);nRows.push_back(403);nCols.push_back(2118);objValueTol.push_back(1.e-10);objValue.push_back(1.7933245380e+06);
  mpsName.push_back("ship04s");min.push_back(true);nRows.push_back(403);nCols.push_back(1458);objValueTol.push_back(1.e-10);objValue.push_back(1.7987147004e+06);
  mpsName.push_back("ship08l");min.push_back(true);nRows.push_back(779);nCols.push_back(4283);objValueTol.push_back(1.e-10);objValue.push_back(1.9090552114e+06);
  mpsName.push_back("ship08s");min.push_back(true);nRows.push_back(779);nCols.push_back(2387);objValueTol.push_back(1.e-10);objValue.push_back(1.9200982105e+06);
  mpsName.push_back("ship12l");min.push_back(true);nRows.push_back(1152);nCols.push_back(5427);objValueTol.push_back(1.e-10);objValue.push_back(1.4701879193e+06);
  mpsName.push_back("ship12s");min.push_back(true);nRows.push_back(1152);nCols.push_back(2763);objValueTol.push_back(1.e-10);objValue.push_back(1.4892361344e+06);
  // CPLEX READ ERROR  mpsName.push_back("sierra");min.push_back(true);nRows.push_back(1228);nCols.push_back(2036);objValueTol.push_back(1.e-10);objValue.push_back(1.5394362184e+07);
  mpsName.push_back("stair");min.push_back(true);nRows.push_back(357);nCols.push_back(467);objValueTol.push_back(1.e-10);objValue.push_back(-2.5126695119e+02);
  mpsName.push_back("standata");min.push_back(true);nRows.push_back(360);nCols.push_back(1075);objValueTol.push_back(1.e-10);objValue.push_back(1.2576995000e+03);
  //mpsName.push_back("standgub");min.push_back(true);nRows.push_back(362);nCols.push_back(1184);objValueTol.push_back(1.e-10);objValue.push_back(1257.6995); 
  mpsName.push_back("standmps");min.push_back(true);nRows.push_back(468);nCols.push_back(1075);objValueTol.push_back(1.e-10);objValue.push_back(1.4060175000E+03); 
  mpsName.push_back("stocfor1");min.push_back(true);nRows.push_back(118);nCols.push_back(111);objValueTol.push_back(1.e-10);objValue.push_back(-4.1131976219E+04);
  mpsName.push_back("stocfor2");min.push_back(true);nRows.push_back(2158);nCols.push_back(2031);objValueTol.push_back(1.e-10);objValue.push_back(-3.9024408538e+04);
  //mpsName.push_back("stocfor3");min.push_back(true);nRows.push_back(16676);nCols.push_back(15695);objValueTol.push_back(1.e-10);objValue.push_back(-3.9976661576e+04);
  //mpsName.push_back("truss");min.push_back(true);nRows.push_back(1001);nCols.push_back(8806);objValueTol.push_back(1.e-10);objValue.push_back(4.5881584719e+05);
  mpsName.push_back("tuff");min.push_back(true);nRows.push_back(334);nCols.push_back(587);objValueTol.push_back(1.e-10);objValue.push_back(2.9214776509e-01);
  mpsName.push_back("vtpbase");min.push_back(true);nRows.push_back(199);nCols.push_back(203);objValueTol.push_back(1.e-10);objValue.push_back(1.2983146246e+05);
  mpsName.push_back("wood1p");min.push_back(true);nRows.push_back(245);nCols.push_back(2594);objValueTol.push_back(1.e-10);objValue.push_back(1.4429024116e+00);
  mpsName.push_back("woodw");min.push_back(true);nRows.push_back(1099);nCols.push_back(8405);objValueTol.push_back(1.e-10);objValue.push_back(1.3044763331E+00);
  
  // Create vector of solver interfaces
  std::vector<OsiSolverInterface*> vecSiP(vecEmptySiP.size());
  
  //assert(dynamic_cast<OsiOslSolverInterface*>(vecEmptySiP[0]) != NULL );
  
  // Loop once for each Mps File
  for (m=0; m<mpsName.size(); m++ ) {
    std::cerr <<"  processing mps file: " <<mpsName[m] 
      <<" (" <<m+1 <<" out of " <<mpsName.size() <<")" <<std::endl;
    
    // Read problem, looping once for each solver interface
    for( i=vecSiP.size()-1; i>=0; --i ) {
      
      // Populate Local SolverInterface vector
      vecSiP[i] = vecEmptySiP[i]->clone();   
      
      // Read data mps file,
      std::string fn = mpsDir+mpsName[m];
      vecSiP[i]->readMps(fn.c_str(),"mps");
      
      // set objective sense,
      if(min[m]) vecSiP[i]->setObjSense(1.0);
      else       vecSiP[i]->setObjSense(-1.0);
      
      // test size of matrix
      int nr=vecSiP[i]->getNumRows();
      int nc=vecSiP[i]->getNumCols();      
      assert( nr==nRows[m]-1 /*|| nr==nRows[m]*/);
      assert(nc==nCols[m]);      
    } 
    
    // Test that Solver Interfaces are returning the same values
    for( i=vecSiP.size()-1; i>0; --i ) {
      OsiPackedVector vim1,vi;
      
      // Test col lowerbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getColLower(),vecSiP[i  ]->getColLower(),
        vecSiP[i  ]->getNumCols() )
        );
      
      // Test col upperbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getColUpper(),vecSiP[i  ]->getColUpper(),
        vecSiP[i  ]->getNumCols() )
        );
      
      // Test obj coefficients
      assert (
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getObjCoefficients(),vecSiP[i  ]->getObjCoefficients(),
        vecSiP[i  ]->getNumCols() )
        );
      
      // Test row lowerbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowLower(),vecSiP[i  ]->getRowLower(),
        vecSiP[i  ]->getNumRows() )
        );
      
      // Test row upperbounds
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowUpper(),vecSiP[i  ]->getRowUpper(),
        vecSiP[i  ]->getNumRows() )
        );
            
      //Test rowsense
      {
        const char * rsm1 = vecSiP[i-1]->getRowSense();
        const char * rs   = vecSiP[i  ]->getRowSense();
        int nr = vecSiP[i]->getNumRows();
        int r;
        for ( r=0; r<nr; r++ )
          assert ( rsm1[r] == rs[r] );
      }
      
      // Test rhs
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRightHandSide(),vecSiP[i  ]->getRightHandSide(),
        vecSiP[i  ]->getNumRows() )
        );
      
      // Test range
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowRange(),vecSiP[i  ]->getRowRange(),
        vecSiP[i  ]->getNumRows() )
        );
      
           
      // Test getObjCoefficients
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getObjCoefficients(),vecSiP[i  ]->getObjCoefficients(),
        vecSiP[i  ]->getNumCols() )
        );

      // Test numels
      assert( vecSiP[i-1]->getNumElements() == vecSiP[i]->getNumElements() );
      
      // Objective Sense
      assert( vecSiP[i-1]->getObjSense() == vecSiP[i  ]->getObjSense() );

      //Test constraint matrix
      {
        const OsiPackedMatrix * rmm1=vecSiP[i-1]->getMatrixByRow();
        const OsiPackedMatrix * rm  =vecSiP[i  ]->getMatrixByRow();
        assert( rmm1->isEquivalent(*rm) );
        
        const OsiPackedMatrix * cmm1=vecSiP[i-1]->getMatrixByCol();
        const OsiPackedMatrix * cm  =vecSiP[i  ]->getMatrixByCol();
        assert( cmm1->isEquivalent(*cm) );
      }
      
    }

    //----------------
          
    
    // Before solving test:
    // isContinuous, isBinary, isInteger, isIngegerNonBinary, isFreeBinary
    // getFractionalIndices
    for( i=vecSiP.size()-1; i>0; --i ) {
      OsiPackedVector vim1,vi;
      
      int c;
      for ( c=0; c<vecSiP[i  ]->getNumCols(); c++ ) {
        assert( !(vecSiP[i-1]->isContinuous(c) ^ vecSiP[i]->isContinuous(c)) );
      }      
      
      //test getFractionalIndices
      {
        OsiVectorInt sm1 = vecSiP[i-1]->getFractionalIndices();
        OsiVectorInt s   = vecSiP[i  ]->getFractionalIndices();
        assert( sm1.size() == s.size() );
        int c;
        for ( c=s.size()-1; c>=0; --c )
          assert ( sm1[c] == s[c] );
      }

      // Test isContinuous, isBinary, isInteger, isIngegerNonBinary, isFreeBinary
      {
        int c;
        int nc = vecSiP[i]->getNumCols();
        for( c=0; c<nc; c++ ){
          assert( ( vecSiP[i-1]->isContinuous(c) &&  vecSiP[i]->isContinuous(c))
            ||    (!vecSiP[i-1]->isContinuous(c) && !vecSiP[i]->isContinuous(c)) );
          assert( ( vecSiP[i-1]->isBinary(c) &&  vecSiP[i]->isBinary(c))
            ||    (!vecSiP[i-1]->isBinary(c) && !vecSiP[i]->isBinary(c)) );
          assert( ( vecSiP[i-1]->isIntegerNonBinary(c) &&  vecSiP[i]->isIntegerNonBinary(c))
            ||    (!vecSiP[i-1]->isIntegerNonBinary(c) && !vecSiP[i]->isIntegerNonBinary(c)) );
          assert( ( vecSiP[i-1]->isFreeBinary(c) &&  vecSiP[i]->isFreeBinary(c))
            ||    (!vecSiP[i-1]->isFreeBinary(c) && !vecSiP[i]->isFreeBinary(c)) );
          assert( ( vecSiP[i-1]->isInteger(c) &&  vecSiP[i]->isInteger(c))
            ||    (!vecSiP[i-1]->isInteger(c) && !vecSiP[i]->isInteger(c)) );
        }
      }

    }

    
    
    // Solve problem, looping once for each solver interface
    for( i=vecSiP.size()-1; i>=0; --i ) { 

#ifdef COIN_USE_VOL
      // Volume Solver Interface can not solve netlib cases, so
      // if volume is being tested we are done.
      {
        OsiVolSolverInterface * vsi = dynamic_cast<OsiVolSolverInterface *>(vecSiP[i]);
        if ( vsi != NULL ) {
          // This one is VolSi, make sure it is the last one.
          // If it isn't then the order of solver interfaces in
          // vecSiP should be changed to ensure that the VolSi is last.
          assert( i == static_cast<int>(vecSiP.size()-1) );
          break;
        }
      }
#endif

      // Solve LP
      vecSiP[i]->initialSolve();
      
      // Test objective solution value
      {
        double soln = vecSiP[i]->getObjValue();       
        OsiRelFltEq eq(objValueTol[m]);
        //std::cerr <<soln <<std::endl;
        assert(eq(soln,objValue[m]));
      }
    }
    
    // Test that Solver Interfaces are returning the same solution values
    for( i=vecSiP.size()-1; i>0; --i ) {
      OsiPackedVector vim1,vi;

#ifdef COIN_USE_VOL
      // Volume Solver Interface can not solve netlib cases, so
      // if volume is being tested we are done.
      {
        OsiVolSolverInterface * vsi = dynamic_cast<OsiVolSolverInterface *>(vecSiP[i]);
        if ( vsi != NULL ) {
          // This one is VolSi, make sure it is the last one.
          // If it isn't then the order of solver interfaces in
          // vecSiP should be changed to ensure that the VolSi is last.
          assert( i == static_cast<int>(vecSiP.size()-1) );
          break;
        }
      }
#endif
 
#if 0
      // Solvers don't seem to give same col solutions.
      // Test col solution
      assert(equivalentVectors(vecSiP[i-1],vecSiP[i], .1 /*objValueTol[m]*/,
			       vecSiP[i-1]->colsol(),vecSiP[i  ]->colsol(),
			       vecSiP[i  ]->getNumCols() ) );
#endif

#if 0
      // Solvers don't seem to return the same row solution
      // Test row solution
      assert(equivalentVectors(vecSiP[i-1],vecSiP[i], .1 /*objValueTol[m]*/,
			       vecSiP[i-1]->rowprice(),vecSiP[i  ]->rowprice(),
			       vecSiP[i  ]->getNumRows() ) );
#endif
      
      int c;
      for ( c=0; c<vecSiP[i  ]->getNumCols(); c++ ) {
        //std::cerr <<c <<" " <<vecSiP[i-1]->isContinuous(c) <<" " <<vecSiP[i  ]->isContinuous(c) <<std::endl;
        assert( !(vecSiP[i-1]->isContinuous(c) ^ vecSiP[i]->isContinuous(c)) );
      }
      
      
      //test getFractionalIndices
      {
        OsiVectorInt sm1 = vecSiP[i-1]->getFractionalIndices();
        OsiVectorInt s   = vecSiP[i  ]->getFractionalIndices();
        assert( sm1.size() == s.size() );
        int c;
        for ( c=s.size()-1; c>=0; --c )
          assert ( sm1[c] == s[c] );
      }

      // Test isContinuous, isBinary, isInteger, isIngegerNonBinary, isFreeBinary
      {
        int c;
        int nc = vecSiP[i]->getNumCols();
        for( c=0; c<nc; c++ ){
          assert( ( vecSiP[i-1]->isContinuous(c) &&  vecSiP[i]->isContinuous(c))
            ||    (!vecSiP[i-1]->isContinuous(c) && !vecSiP[i]->isContinuous(c)) );
          assert( ( vecSiP[i-1]->isBinary(c) &&  vecSiP[i]->isBinary(c))
            ||    (!vecSiP[i-1]->isBinary(c) && !vecSiP[i]->isBinary(c)) );
          assert( ( vecSiP[i-1]->isIntegerNonBinary(c) &&  vecSiP[i]->isIntegerNonBinary(c))
            ||    (!vecSiP[i-1]->isIntegerNonBinary(c) && !vecSiP[i]->isIntegerNonBinary(c)) );
          assert( ( vecSiP[i-1]->isFreeBinary(c) &&  vecSiP[i]->isFreeBinary(c))
            ||    (!vecSiP[i-1]->isFreeBinary(c) && !vecSiP[i]->isFreeBinary(c)) );
          assert( ( vecSiP[i-1]->isInteger(c) &&  vecSiP[i]->isInteger(c))
            ||    (!vecSiP[i-1]->isInteger(c) && !vecSiP[i]->isInteger(c)) );
        }
      }

    }
    
    
    // Delete locally created solver interfaces
    for( i=vecSiP.size()-1; i>=0; --i)
      delete vecSiP[i];
  }
  
}

//#############################################################################
//#############################################################################

static bool testIntParam(OsiSolverInterface * si, int k, int val)
{
  int i = 123456789, orig = 123456789;
  bool ret;
  OsiIntParam key = static_cast<OsiIntParam>(k);
  si->getIntParam(key, orig);
  if (si->setIntParam(key, val)) {
    ret = (si->getIntParam(key, i) == true) && (i == val);
  } else {
    ret = (si->getIntParam(key, i) == true) && (i == orig);
  }
  return ret;
}
static bool testDblParam(OsiSolverInterface * si, int k, double val)
{
  double d = 123456789.0, orig = 123456789.0;
  bool ret;
  OsiDblParam key = static_cast<OsiDblParam>(k);
  si->getDblParam(key, orig);
  if (si->setDblParam(key, val)) {
    ret = (si->getDblParam(key, d) == true) && (d == val);
  } else {
    ret = (si->getDblParam(key, d) == true) && (d == orig);
  }
  return ret;
}

//#############################################################################
//#############################################################################

void
OsiSolverInterfaceCommonUnitTest(const OsiSolverInterface* emptySi,
				 const std::string & mpsDir)
{
  
  int i;
  OsiRelFltEq eq;
  std::string fn = mpsDir+"exmip1";
  OsiSolverInterface * exmip1Si = emptySi->clone(); 
  exmip1Si->readMps(fn.c_str(),"mps");

#if 0
  // Solver Interfaces do not presently pass this set of tests.
  // Something should be done to make them consistent.
  // When there is an empty solver interface should
  // getColLower() return NULL or a zero length vector??
  // similar question for other methods.

  // Test that values returned from an empty solverInterface
  {
    OsiSolverInterface * si = emptySi->clone();
    assert( si->getNumRows()==0 );
    assert( si->getNumCols()==0 );
    assert( si->getNumElements()==0 );
    assert( si->getColLower()==NULL );
    assert( si->getColUpper()==NULL );
    assert( si->colsol()==NULL );
    assert( si->getObjCoefficients()==NULL );
    assert( si->getRowRange()==NULL );
    assert( si->getRightHandSide()!=NULL );
    assert( si->getRowSense()==NULL );
    assert( si->getRowLower()==NULL );
    assert( si->getRowUpper()==NULL );
    delete si;
  }
#endif
  
  // Test that problem was loaded correctly
  {
    const char   * exmip1Sirs  = exmip1Si->getRowSense();
    assert( exmip1Sirs[0]=='G' );
    assert( exmip1Sirs[1]=='L' );
    assert( exmip1Sirs[2]=='E' );
    assert( exmip1Sirs[3]=='R' );
    assert( exmip1Sirs[4]=='R' );
    
    const double * exmip1Sirhs = exmip1Si->getRightHandSide();
    assert( eq(exmip1Sirhs[0],2.5) );
    assert( eq(exmip1Sirhs[1],2.1) );
    assert( eq(exmip1Sirhs[2],4.0) );
    assert( eq(exmip1Sirhs[3],5.0) );
    assert( eq(exmip1Sirhs[4],15.) ); 
    
    const double * exmip1Sirr  = exmip1Si->getRowRange();
    assert( eq(exmip1Sirr[0],0.0) );
    assert( eq(exmip1Sirr[1],0.0) );
    assert( eq(exmip1Sirr[2],0.0) );
    assert( eq(exmip1Sirr[3],5.0-1.8) );
    assert( eq(exmip1Sirr[4],15.0-3.0) );
    
    OsiPackedMatrix pm;
    pm.setExtraGap(0.0);
    pm.setExtraMajor(0.0);
    pm = *exmip1Si->getMatrixByRow();
    
    const double * ev = pm.getElements();
    assert( eq(ev[0],   3.0) );
    assert( eq(ev[1],   1.0) );
    assert( eq(ev[2],  -2.0) );
    assert( eq(ev[3],  -1.0) );
    assert( eq(ev[4],  -1.0) );
    assert( eq(ev[5],   2.0) );
    assert( eq(ev[6],   1.1) );
    assert( eq(ev[7],   1.0) );
    assert( eq(ev[8],   1.0) );
    assert( eq(ev[9],   2.8) );
    assert( eq(ev[10], -1.2) );
    assert( eq(ev[11],  5.6) );
    assert( eq(ev[12],  1.0) );
    assert( eq(ev[13],  1.9) );
    
    const int * mi = pm.getVectorStarts();
    assert( mi[0]==0 );
    assert( mi[1]==5 );
    assert( mi[2]==7 );
    assert( mi[3]==9 );
    assert( mi[4]==11 );
    assert( mi[5]==14 );
    
    const int * ei = pm.getIndices();
    assert( ei[0]  ==  0 );
    assert( ei[1]  ==  1 );
    assert( ei[2]  ==  3 );
    assert( ei[3]  ==  4 );
    assert( ei[4]  ==  7 );
    assert( ei[5]  ==  1 );
    assert( ei[6]  ==  2 );
    assert( ei[7]  ==  2 );
    assert( ei[8]  ==  5 );
    assert( ei[9]  ==  3 );
    assert( ei[10] ==  6 );
    assert( ei[11] ==  0 );
    assert( ei[12] ==  4 );
    assert( ei[13] ==  7 );    
    
    assert( pm.getMajorDim() == 5 ); 
    assert( pm.getNumElements() == 14 );
    
    int nc = exmip1Si->getNumCols();
    int nr = exmip1Si->getNumRows();
    const double * cl = exmip1Si->getColLower();
    const double * cu = exmip1Si->getColUpper();
    const double * rl = exmip1Si->getRowLower();
    const double * ru = exmip1Si->getRowUpper();
    assert( nc == 8 );
    assert( nr == 5 );
    assert( eq(cl[0],2.5) );
    assert( eq(cl[1],0.0) );
    assert( eq(cl[2],0.0) );
    assert( eq(cl[3],0.0) );
    assert( eq(cl[4],0.5) );
    assert( eq(cl[5],0.0) );
    assert( eq(cl[6],0.0) );
    assert( eq(cl[7],0.0) );
    assert( eq(cu[0],exmip1Si->getInfinity()) );
    assert( eq(cu[1],4.1) );
    assert( eq(cu[2],1.0) );
    assert( eq(cu[3],1.0) );
    assert( eq(cu[4],4.0) );
    assert( eq(cu[5],exmip1Si->getInfinity()) );
    assert( eq(cu[6],exmip1Si->getInfinity()) );
    assert( eq(cu[7],4.3) );

    assert( eq(rl[0],2.5) );
    assert( eq(rl[1],-exmip1Si->getInfinity()) );
    assert( eq(rl[2],4.0) );
    assert( eq(rl[3],1.8) );
    assert( eq(rl[4],3.0) );
    assert( eq(ru[0],exmip1Si->getInfinity()) );
    assert( eq(ru[1],2.1) );
    assert( eq(ru[2],4.0) );
    assert( eq(ru[3],5.0) );
    assert( eq(ru[4],15.0) );
    
    //const double * cs = exmip1Si->colsol();
    //assert( eq(cs[0],0.0) );
    //assert( eq(cs[7],0.0) );
    //assert( eq(exmip1Si->getObjValue(),0.0) );
    
    assert( eq( exmip1Si->getObjCoefficients()[0],  1.0) );
    assert( eq( exmip1Si->getObjCoefficients()[1],  0.0) );
    assert( eq( exmip1Si->getObjCoefficients()[2],  0.0) );
    assert( eq( exmip1Si->getObjCoefficients()[3],  0.0) );
    assert( eq( exmip1Si->getObjCoefficients()[4],  2.0) );
    assert( eq( exmip1Si->getObjCoefficients()[5],  0.0) );
    assert( eq( exmip1Si->getObjCoefficients()[6],  0.0) );
    assert( eq( exmip1Si->getObjCoefficients()[7], -1.0) );
  }
    
  // Test load and assign problem
  {   
    {
      OsiSolverInterface * base = exmip1Si->clone();
      OsiSolverInterface *  si1 = emptySi->clone(); 
      OsiSolverInterface *  si2 = emptySi->clone(); 
      OsiSolverInterface *  si3 = emptySi->clone(); 
      OsiSolverInterface *  si4 = emptySi->clone();  
      OsiSolverInterface *  si5 = emptySi->clone();  
      OsiSolverInterface *  si6 = emptySi->clone(); 
      OsiSolverInterface *  si7 = emptySi->clone();  
      OsiSolverInterface *  si8 = emptySi->clone(); 
        
      si1->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowSense(),base->getRightHandSide(),
		       base->getRowRange());
      si2->loadProblem(*base->getMatrixByRow(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowSense(),base->getRightHandSide(),
		       base->getRowRange());
      si3->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowLower(),base->getRowUpper() );
      si4->loadProblem(*base->getMatrixByCol(),
		       base->getColLower(),base->getColUpper(),
		       base->getObjCoefficients(),
		       base->getRowLower(),base->getRowUpper() );
        
      OsiPackedMatrix * pm = new OsiPackedMatrix(*base->getMatrixByCol());
      double * clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      double * cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      double * objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      double * rlb = new double[base->getNumRows()];
      std::copy(base->getRowLower(),
		base->getRowLower()+base->getNumRows(),rlb);
      double * rub = new double[base->getNumRows()];
      std::copy(base->getRowUpper(),
		base->getRowUpper()+base->getNumRows(),rub);
      si5->assignProblem(pm,clb,cub,objc,rlb,rub);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rlb==NULL);
      assert(rub==NULL);
        
      pm = new OsiPackedMatrix(*base->getMatrixByRow());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      rlb = new double[base->getNumRows()];
      std::copy(base->getRowLower(),
		base->getRowLower()+base->getNumRows(),rlb);
      rub = new double[base->getNumRows()];
      std::copy(base->getRowUpper(),
		base->getRowUpper()+base->getNumRows(),rub);
      si6->assignProblem(pm,clb,cub,objc,rlb,rub);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rlb==NULL);
      assert(rub==NULL);      
        
      pm = new OsiPackedMatrix(*base->getMatrixByCol());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      char * rsen = new char[base->getNumRows()];
      std::copy(base->getRowSense(),
		base->getRowSense()+base->getNumRows(),rsen);
      double * rhs = new double[base->getNumRows()];
      std::copy(base->getRightHandSide(),
		base->getRightHandSide()+base->getNumRows(),rhs);
      double * rng = new double[base->getNumRows()];
      std::copy(base->getRowRange(),
		base->getRowRange()+base->getNumRows(),rng);
      si7->assignProblem(pm,clb,cub,objc,rsen,rhs,rng);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rsen==NULL);
      assert(rhs==NULL);
      assert(rng==NULL);
        
      pm = new OsiPackedMatrix(*base->getMatrixByCol());
      clb = new double[base->getNumCols()];
      std::copy(base->getColLower(),
		base->getColLower()+base->getNumCols(),clb);
      cub = new double[base->getNumCols()];
      std::copy(base->getColUpper(),
		base->getColUpper()+base->getNumCols(),cub);
      objc = new double[base->getNumCols()];
      std::copy(base->getObjCoefficients(),
		base->getObjCoefficients()+base->getNumCols(),objc);
      rsen = new char[base->getNumRows()];
      std::copy(base->getRowSense(),
		base->getRowSense()+base->getNumRows(),rsen);
      rhs = new double[base->getNumRows()];
      std::copy(base->getRightHandSide(),
		base->getRightHandSide()+base->getNumRows(),rhs);
      rng = new double[base->getNumRows()];
      std::copy(base->getRowRange(),
		base->getRowRange()+base->getNumRows(),rng);
      si8->assignProblem(pm,clb,cub,objc,rsen,rhs,rng);
      assert(pm==NULL);
      assert(clb==NULL);
      assert(cub==NULL);
      assert(objc==NULL);
      assert(rsen==NULL);
      assert(rhs==NULL);
      assert(rng==NULL);
     
        
      // Create an indices vector        
      OsiPackedVector basePv,pv;
      assert(base->getNumCols()<10);
      assert(base->getNumRows()<10);
      int indices[10];
      int i;
      for (i=0; i<10; i++) indices[i]=i;
        
      // Test solve methods.
      try {
	base->initialSolve();
	si1->initialSolve();
	si2->initialSolve();
	si3->initialSolve();
	si4->initialSolve();
	si5->initialSolve();
	si6->initialSolve();
	si7->initialSolve();
	si8->initialSolve();
      }        
      catch (CoinError e) {
#ifdef COIN_USE_VOL
	// Vol solver interface is expected to throw
	// an error if the data has a ranged row.
          
	// Check that using Vol SI
	OsiVolSolverInterface * vsi =
	  dynamic_cast<OsiVolSolverInterface *>(base);
	assert( vsi != NULL );        
          
	// Test that there is non-zero range
	basePv.setFull(base->getNumRows(),base->getRowRange());
	pv.setConstant( base->getNumRows(), indices, 0.0 );
	assert(!basePv.isEquivalent(pv)); 
#else
	assert(0==1);
#endif
      }  
        
      // Test collower
      basePv.setVector(base->getNumCols(),indices,base->getColLower());
      pv.setVector( si1->getNumCols(),indices, si1->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getColLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getColLower());
      assert(basePv.isEquivalent(pv));
        
      // Test colupper
      basePv.setVector(base->getNumCols(),indices,base->getColUpper());
      pv.setVector( si1->getNumCols(),indices, si1->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getColUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getColUpper());
      assert(basePv.isEquivalent(pv));
        
      // Test getObjCoefficients
      basePv.setVector(base->getNumCols(),indices,base->getObjCoefficients());
      pv.setVector( si1->getNumCols(),indices, si1->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumCols(),indices, si2->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumCols(),indices, si3->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumCols(),indices, si4->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumCols(),indices, si5->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumCols(),indices, si6->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumCols(),indices, si7->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumCols(),indices, si8->getObjCoefficients());
      assert(basePv.isEquivalent(pv));
                
      // Test rowrhs
      basePv.setFull(base->getNumRows(),base->getRightHandSide());
      pv.setFull( si1->getNumRows(), si1->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si2->getNumRows(), si2->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si3->getNumRows(), si3->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si4->getNumRows(), si4->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si5->getNumRows(), si5->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si6->getNumRows(), si6->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si7->getNumRows(), si7->getRightHandSide());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si8->getNumRows(), si8->getRightHandSide());
      assert(basePv.isEquivalent(pv));
        
      // Test rowrange
      basePv.setFull(base->getNumRows(),base->getRowRange());
      pv.setFull( si1->getNumRows(), si1->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si2->getNumRows(), si2->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si3->getNumRows(), si3->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si4->getNumRows(), si4->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si5->getNumRows(), si5->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si6->getNumRows(), si6->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si7->getNumRows(), si7->getRowRange());
      assert(basePv.isEquivalent(pv));
      pv.setFull( si8->getNumRows(), si8->getRowRange());
      assert(basePv.isEquivalent(pv));
        
      // Test row sense
      {
	const char * cb = base->getRowSense();
	const char * c1 = si1->getRowSense();
	const char * c2 = si2->getRowSense();
	const char * c3 = si3->getRowSense();
	const char * c4 = si4->getRowSense();
	const char * c5 = si5->getRowSense();
	const char * c6 = si6->getRowSense();
	const char * c7 = si7->getRowSense();
	const char * c8 = si8->getRowSense();
	int nr = base->getNumRows();
	for ( i=0; i<nr; i++ ) {
	  assert( cb[i]==c1[i] );
	  assert( cb[i]==c2[i] );
	  assert( cb[i]==c3[i] );
	  assert( cb[i]==c4[i] );
	  assert( cb[i]==c5[i] );
	  assert( cb[i]==c6[i] );
	  assert( cb[i]==c7[i] );
	  assert( cb[i]==c8[i] );
	}
      }
        
      // Test rowlower
      basePv.setVector(base->getNumRows(),indices,base->getRowLower());
      pv.setVector( si1->getNumRows(),indices, si1->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumRows(),indices, si2->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumRows(),indices, si3->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumRows(),indices, si4->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumRows(),indices, si5->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumRows(),indices, si6->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumRows(),indices, si7->getRowLower());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumRows(),indices, si8->getRowLower());
      assert(basePv.isEquivalent(pv));
        
      // Test rowupper
      basePv.setVector(base->getNumRows(),indices,base->getRowUpper());
      pv.setVector( si1->getNumRows(),indices, si1->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si2->getNumRows(),indices, si2->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si3->getNumRows(),indices, si3->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si4->getNumRows(),indices, si4->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si5->getNumRows(),indices, si5->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si6->getNumRows(),indices, si6->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si7->getNumRows(),indices, si7->getRowUpper());
      assert(basePv.isEquivalent(pv));
      pv.setVector( si8->getNumRows(),indices, si8->getRowUpper());
      assert(basePv.isEquivalent(pv)); 
        
      // Test Constraint Matrix
      assert( base->getMatrixByCol()->isEquivalent(*si1->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si1->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si2->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si2->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si3->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si3->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si4->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si4->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si5->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si5->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si6->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si6->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si7->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si7->getMatrixByRow()) );
      assert( base->getMatrixByCol()->isEquivalent(*si8->getMatrixByCol()) );
      assert( base->getMatrixByRow()->isEquivalent(*si8->getMatrixByRow()) );
        
      // Test Objective Value
      assert( eq(base->getObjValue(),si1->getObjValue()) );
      assert( eq(base->getObjValue(),si2->getObjValue()) );
      assert( eq(base->getObjValue(),si3->getObjValue()) );
      assert( eq(base->getObjValue(),si4->getObjValue()) );
      assert( eq(base->getObjValue(),si5->getObjValue()) );
      assert( eq(base->getObjValue(),si6->getObjValue()) );
      assert( eq(base->getObjValue(),si7->getObjValue()) );
      assert( eq(base->getObjValue(),si8->getObjValue()) );
        
      // Clean-up
      delete si8;
      delete si7;
      delete si6;
      delete si5;
      delete si4;
      delete si3;
      delete si2;
      delete si1;
      delete base;
    }     
    // Test load/assign with null parms
    {
      //Load problem with row bounds and all rims at defaults
      {
	OsiSolverInterface *  si = emptySi->clone(); 
          
	si->loadProblem(*exmip1Si->getMatrixByCol(),NULL,NULL,NULL,NULL,NULL);
          
	// Test column settings
	assert(si->getNumCols()==exmip1Si->getNumCols() );
	for ( i=0; i<si->getNumCols(); i++ ) {
	  assert( eq(si->getColLower()[i],0.0) );
	  assert( eq(si->getColUpper()[i],si->getInfinity()) );
	  assert( eq(si->getObjCoefficients()[i],0.0) );
	}
	// Test row settings
	assert(si->getNumRows()==exmip1Si->getNumRows() );
	const double * rh = si->getRightHandSide();
	const double * rr = si->getRowRange();
	const char * rs = si->getRowSense();
	const double * rl = si->getRowLower();
	const double * ru = si->getRowUpper();
	for ( i=0; i<si->getNumRows(); i++ ) {
	  assert( eq(rh[i],0.0) );
	  assert( eq(rr[i],0.0) );
	  assert( 'N'==rs[i] );
	  assert( eq(rl[i],-si->getInfinity()) );
	  assert( eq(ru[i], si->getInfinity()) );
	}
          
	delete si;
      }
      //Load problem with row rhs and all rims at defaults
      {
	OsiSolverInterface *  si = emptySi->clone(); 
          
	si->loadProblem(*exmip1Si->getMatrixByRow(),
			NULL,NULL,NULL,
			exmip1Si->getRowSense(),
			exmip1Si->getRightHandSide(),
			exmip1Si->getRowRange());
	// Test column settings
	assert(si->getNumCols()==exmip1Si->getNumCols() );
	for ( i=0; i<si->getNumCols(); i++ ) {
	  assert( eq(si->getColLower()[i],0.0) );
	  assert( eq(si->getColUpper()[i],si->getInfinity()) );
	  assert( eq(si->getObjCoefficients()[i],0.0) );
	}
	// Test row settings
	assert(si->getNumRows()==exmip1Si->getNumRows() );
	for ( i=0; i<si->getNumRows(); i++ ) {
	  char s = si->getRowSense()[i];
	  assert( eq(si->getRightHandSide()[i],
		     exmip1Si->getRightHandSide()[i]) );
	  assert( eq(si->getRowRange()[i],
		     exmip1Si->getRowRange()[i]) );
	  assert( s==exmip1Si->getRowSense()[i] );
            
	  if ( s=='G' ) {
	    assert( eq(si->getRowLower()[i],
		       exmip1Si->getRightHandSide()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       si->getInfinity()) );
	  }
	  else if ( s=='L' ) {
	    assert( eq(si->getRowLower()[i],
		       -si->getInfinity()) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	  else if ( s=='E' ) {
	    assert( eq(si->getRowLower()[i],
		       si->getRowUpper()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	  else if ( s=='N' ) {
	    assert( eq(si->getRowLower()[i], -si->getInfinity()) );
	    assert( eq(si->getRowUpper()[i],  si->getInfinity()) );
	  }
	  else if ( s=='R' ) {
	    assert( eq(si->getRowLower()[i],
		       exmip1Si->getRightHandSide()[i] -
		       exmip1Si->getRowRange()[i]) );
	    assert( eq(si->getRowUpper()[i],
		       exmip1Si->getRightHandSide()[i]) );
	  }
	}
          
	delete si;
      }
    }
  }

  // Add a Laci suggested test case
  // Load in a problem as column ordered matrix, 
  // extract the row ordered copy, 
  // add a row, 
  // extract the row ordered copy again and test whether it's ok. 
  // (the same can be done with reversing the role
  //  of row and column ordered.)
  {
    OsiSolverInterface *  si = emptySi->clone(); 
      
    si->loadProblem(
		    *(exmip1Si->getMatrixByCol()),
		    exmip1Si->getColLower(),
		    exmip1Si->getColUpper(),
		    exmip1Si->getObjCoefficients(),
		    exmip1Si->getRowSense(),
		    exmip1Si->getRightHandSide(),
		    exmip1Si->getRowRange() );

    OsiPackedMatrix pm1 = *(si->getMatrixByRow());

    // Get a row of the matrix to make a cut
    OsiPackedVector pv =exmip1Si->getMatrixByRow()->getVector(1);
    pv.setElement(0,3.14*pv.getElements()[0]);

    OsiRowCut rc;
    rc.setRow( pv );
    rc.setLb( exmip1Si->getRowLower()[1]-0.5 );
    rc.setUb( exmip1Si->getRowUpper()[1]-0.5 );

    OsiCuts cuts;
    cuts.insert(rc);

    si->applyCuts(cuts);
      
    OsiPackedMatrix pm2 = *(si->getMatrixByRow());

    assert(pm1.getNumRows()==pm2.getNumRows()-1);
    int i;
    for( i=0; i<pm1.getNumRows(); ++i ) {
      assert( pm1.getVector(i) == pm2.getVector(i) );
    }
    // Test that last row of pm2 is same as added cut
    assert( pm2.getVector(pm2.getNumRows()-1).isEquivalent(pv) );

    delete si;
  }
  {
    OsiSolverInterface *  si = emptySi->clone(); 
      
    si->loadProblem(
		    *(exmip1Si->getMatrixByRow()),
		    exmip1Si->getColLower(),
		    exmip1Si->getColUpper(),
		    exmip1Si->getObjCoefficients(),
		    exmip1Si->getRowLower(),
		    exmip1Si->getRowUpper() );

    OsiPackedMatrix pm1 = *(si->getMatrixByCol());

    // Get a row of the matrix to make a cut
    OsiPackedVector pv =exmip1Si->getMatrixByRow()->getVector(1);
    pv.setElement(0,3.14*pv.getElements()[0]);

    OsiRowCut rc;
    rc.setRow( pv );
    rc.setLb( exmip1Si->getRowLower()[1]-0.5 );
    rc.setUb( exmip1Si->getRowUpper()[1]-0.5 );

    OsiCuts cuts;
    cuts.insert(rc);

    si->applyCuts(cuts);
      
    OsiPackedMatrix pm2 = *(si->getMatrixByCol());

    assert( pm1.isColOrdered() );
    assert( pm2.isColOrdered() );
    assert( pm1.getNumRows()==pm2.getNumRows()-1 );

    OsiPackedMatrix pm1ByRow;
    pm1ByRow.reverseOrderedCopyOf(pm1);
    OsiPackedMatrix pm2ByRow;
    pm2ByRow.reverseOrderedCopyOf(pm2);

    assert( !pm1ByRow.isColOrdered() );
    assert( !pm2ByRow.isColOrdered() );      
    assert( pm1ByRow.getNumRows()==pm2ByRow.getNumRows()-1 );
    assert( pm1.getNumRows() == pm1ByRow.getNumRows() );
    assert( pm2.getNumRows() == pm2ByRow.getNumRows() );

    int i;
    for( i=0; i<pm1ByRow.getNumRows(); ++i ) {
      assert( pm1ByRow.getVector(i) == pm2ByRow.getVector(i) );
    }
    // Test that last row of pm2 is same as added cut
    assert( pm2ByRow.getVector(pm2ByRow.getNumRows()-1).isEquivalent(pv) );

    delete si;
  }
    
  delete exmip1Si;

  {
    // Testing parameter settings
    OsiSolverInterface *  si = emptySi->clone();
    int i;
    int ival;
    double dval;
    assert(si->getIntParam(OsiLastIntParam, ival) == false);
    assert(si->getDblParam(OsiLastDblParam, dval) == false);
    assert(si->setIntParam(OsiLastIntParam, 0) == false);
    assert(si->setDblParam(OsiLastDblParam, 0) == false);

    for (i = 0; i < OsiLastIntParam; ++i) {
      const bool exists = si->getIntParam(static_cast<OsiIntParam>(i), ival);
      // existence and test should result in the same
      assert(!exists ^ testIntParam(si, i, -1));
      assert(!exists ^ testIntParam(si, i, 0));
      assert(!exists ^ testIntParam(si, i, 1));
      assert(!exists ^ testIntParam(si, i, 9999999));
      assert(!exists ^ testIntParam(si, i, INT_MAX));
      if (exists)
	assert(si->getIntParam(static_cast<OsiIntParam>(i), ival));
    }
    for (i = 0; i < OsiLastDblParam; ++i) {
      const bool exists = si->getDblParam(static_cast<OsiDblParam>(i), dval);
      // existence and test should result in the same
      assert(!exists ^ testDblParam(si, i, -1e50));
      assert(!exists ^ testDblParam(si, i, -1e10));
      assert(!exists ^ testDblParam(si, i, -1));
      assert(!exists ^ testDblParam(si, i, -1e-4));
      assert(!exists ^ testDblParam(si, i, -1e-15));
      assert(!exists ^ testDblParam(si, i, 1e50));
      assert(!exists ^ testDblParam(si, i, 1e10));
      assert(!exists ^ testDblParam(si, i, 1));
      assert(!exists ^ testDblParam(si, i, 1e-4));
      assert(!exists ^ testDblParam(si, i, 1e-15));
      if (exists)
	assert(si->setDblParam(static_cast<OsiDblParam>(i), dval));
	
    }
    delete si;
  }

  /*
    With this matrix we have a primal/dual infeas problem. Leaving the first
    row makes it primal feas, leaving the first col makes it dual feas.
    All vars are >= 0

    obj: -1  2 -3  4 -5 (min)

          0 -1  0  0 -2  >=  1
          1  0 -3  0  4  >= -2
          0  3  0 -5  0  >=  3
          0  0  5  0 -6  >= -4
          2 -4  0  6  0  >=  5
  */
}

