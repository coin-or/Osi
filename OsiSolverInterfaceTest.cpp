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
#ifdef COIN_USE_DYLP
#  include "OsiDyLpSolverInterface.hpp"
#endif
#include "OsiXprSolverInterface.hpp"
#include "OsiCpxSolverInterface.hpp"
#include "OsiOslSolverInterface.hpp"
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


/*! \brief Run solvers on NetLib problems.

  The routine creates a vector of NetLib problems (problem name, objective,
  various other characteristics), and a vector of solvers to be tested.
  
  Each solver is run on each problem. The run is deemed successful if the
  solver reports the correct problem size after loading and returns the
  correct objective value after optimization.

  If multiple solvers are available, the results are compared pairwise against
  the results reported by adjacent solvers in the solver vector. Due to
  limitations of the volume solver, it must be the last solver in vecEmptySiP.
*/

void OsiSolverInterfaceMpsUnitTest
  (const std::vector<OsiSolverInterface*> & vecEmptySiP,
   const std::string & mpsDir)

{ int i ;
  unsigned int m ;

/*
  Vectors to hold test problem names and characteristics. The objective value
  after optimization (objValue) must agree to the specified tolerance
  (objValueTol).
*/
  std::vector<std::string> mpsName ;
  std::vector<bool> min ;
  std::vector<int> nRows ;
  std::vector<int> nCols ;
  std::vector<double> objValue ;
  std::vector<double> objValueTol ;
/*
  And a macro to make the vector creation marginally readable.
*/
#define PUSH_MPS(zz_mpsName_zz,zz_min_zz,\
		 zz_nRows_zz,zz_nCols_zz,zz_objValue_zz,zz_objValueTol_zz) \
  mpsName.push_back(zz_mpsName_zz) ; \
  min.push_back(zz_min_zz) ; \
  nRows.push_back(zz_nRows_zz) ; \
  nCols.push_back(zz_nCols_zz) ; \
  objValueTol.push_back(zz_objValueTol_zz) ; \
  objValue.push_back(zz_objValue_zz) ;

/*
  Load up the problem vector. Note that the row counts here include the
  objective function.
*/
  PUSH_MPS("25fv47",true,822,1571,5.5018458883E+03,1.0e-10)
  PUSH_MPS("80bau3b",true,2263,9799,9.8722419241E+05,1.e-10)
  PUSH_MPS("adlittle",true,57,97,2.2549496316e+05,1.e-10)
  PUSH_MPS("afiro",true,28,32,-4.6475314286e+02,1.e-10)
  PUSH_MPS("agg",true,489,163,-3.5991767287e+07,1.e-10)
  PUSH_MPS("agg2",true,517,302,-2.0239252356e+07,1.e-10)
  PUSH_MPS("agg3",true,517,302,1.0312115935e+07,1.e-10)
  PUSH_MPS("bandm",true,306,472,-1.5862801845e+02,1.e-10)
  PUSH_MPS("beaconfd",true,174,262,3.3592485807e+04,1.e-10)
  PUSH_MPS("blend",true,75,83,-3.0812149846e+01,1.e-10)
  PUSH_MPS("bnl1",true,644,1175,1.9776295615E+03,1.e-10)
  PUSH_MPS("bnl2",true,2325,3489,1.8112365404e+03,1.e-10)
  PUSH_MPS("boeing1",true,/*351*/352,384,-3.3521356751e+02,1.e-10)
  PUSH_MPS("boeing2",true,167,143,-3.1501872802e+02,1.e-10)
  PUSH_MPS("bore3d",true,234,315,1.3730803942e+03,1.e-10)
  PUSH_MPS("brandy",true,221,249,1.5185098965e+03,1.e-10)
  PUSH_MPS("capri",true,272,353,2.6900129138e+03,1.e-10)
  PUSH_MPS("cycle",true,1904,2857,-5.2263930249e+00,1.e-9)
  PUSH_MPS("czprob",true,930,3523,2.1851966989e+06,1.e-10)
  PUSH_MPS("d2q06c",true,2172,5167,122784.21557456,1.e-7)
  PUSH_MPS("d6cube",true,416,6184,3.1549166667e+02,1.e-10)
  PUSH_MPS("degen2",true,445,534,-1.4351780000e+03,1.e-10)
  PUSH_MPS("degen3",true,1504,1818,-9.8729400000e+02,1.e-10)
  PUSH_MPS("dfl001",true,6072,12230,1.1266396047E+07,1.e-5)
  PUSH_MPS("e226",true,224,282,(-18.751929066+7.113),1.e-10) // NOTE: Objective function has constant of 7.113
  PUSH_MPS("etamacro",true,401,688,-7.5571521774e+02 ,1.e-6)
  PUSH_MPS("fffff800",true,525,854,5.5567961165e+05,1.e-6)
  PUSH_MPS("finnis",true,498,614,1.7279096547e+05,1.e-6)
  PUSH_MPS("fit1d",true,25,1026,-9.1463780924e+03,1.e-10)
  PUSH_MPS("fit1p",true,628,1677,9.1463780924e+03,1.e-10)
  PUSH_MPS("fit2d",true,26,10500,-6.8464293294e+04,1.e-10)
  PUSH_MPS("fit2p",true,3001,13525,6.8464293232e+04,1.e-9)
  PUSH_MPS("forplan",true,162,421,-6.6421873953e+02,1.e-6)
  PUSH_MPS("ganges",true,1310,1681,-1.0958636356e+05,1.e-5)
  PUSH_MPS("gfrd-pnc",true,617,1092,6.9022359995e+06,1.e-10)
  PUSH_MPS("greenbea",true,2393,5405,/*-7.2462405908e+07*/-72555248.129846,1.e-10)
  PUSH_MPS("greenbeb",true,2393,5405,/*-4.3021476065e+06*/-4302260.2612066,1.e-10)
  PUSH_MPS("grow15",true,301,645,-1.0687094129e+08,1.e-10)
  PUSH_MPS("grow22",true,441,946,-1.6083433648e+08,1.e-10)
  PUSH_MPS("grow7",true,141,301,-4.7787811815e+07,1.e-10)
  PUSH_MPS("israel",true,175,142,-8.9664482186e+05,1.e-10)
  PUSH_MPS("kb2",true,44,41,-1.7499001299e+03,1.e-10)
  PUSH_MPS("lotfi",true,154,308,-2.5264706062e+01,1.e-10)
  PUSH_MPS("maros",true,847,1443,-5.8063743701e+04,1.e-10)
  PUSH_MPS("maros-r7",true,3137,9408,1.4971851665e+06,1.e-10)
  PUSH_MPS("modszk1",true,688,1620,3.2061972906e+02,1.e-10)
  PUSH_MPS("nesm",true,663,2923,1.4076073035e+07,1.e-5)
  PUSH_MPS("perold",true,626,1376,-9.3807580773e+03,1.e-6)
  PUSH_MPS("pilot",true,1442,3652,/*-5.5740430007e+02*/-557.48972927292,1.e-5)
  PUSH_MPS("pilot4",true,411,1000,-2.5811392641e+03,1.e-8)
  PUSH_MPS("pilot87",true,2031,4883,3.0171072827e+02,1.e-4)
  PUSH_MPS("pilotnov",true,976,2172,-4.4972761882e+03,1.e-10)
  // ?? PUSH_MPS("qap8",true,913,1632,2.0350000000e+02,1.e-10)
  // ?? PUSH_MPS("qap12",true,3193,8856,5.2289435056e+02,1.e-10)
  // ?? PUSH_MPS("qap15",true,6331,22275,1.0409940410e+03,1.e-10)
  PUSH_MPS("recipe",true,92,180,-2.6661600000e+02,1.e-10)
  PUSH_MPS("sc105",true,106,103,-5.2202061212e+01,1.e-10)
  PUSH_MPS("sc205",true,206,203,-5.2202061212e+01,1.e-10)
  PUSH_MPS("sc50a",true,51,48,-6.4575077059e+01,1.e-10)
  PUSH_MPS("sc50b",true,51,48,-7.0000000000e+01,1.e-10)
  PUSH_MPS("scagr25",true,472,500,-1.4753433061e+07,1.e-10)
  PUSH_MPS("scagr7",true,130,140,-2.3313892548e+06,1.e-6)
  PUSH_MPS("scfxm1",true,331,457,1.8416759028e+04,1.e-10)
  PUSH_MPS("scfxm2",true,661,914,3.6660261565e+04,1.e-10)
  PUSH_MPS("scfxm3",true,991,1371,5.4901254550e+04,1.e-10)
  PUSH_MPS("scorpion",true,389,358,1.8781248227e+03,1.e-10)
  PUSH_MPS("scrs8",true,491,1169,9.0429998619e+02,1.e-5)
  PUSH_MPS("scsd1",true,78,760,8.6666666743e+00,1.e-10)
  PUSH_MPS("scsd6",true,148,1350,5.0500000078e+01,1.e-10)
  PUSH_MPS("scsd8",true,398,2750,9.0499999993e+02,1.e-10)
  PUSH_MPS("sctap1",true,301,480,1.4122500000e+03,1.e-10)
  PUSH_MPS("sctap2",true,1091,1880,1.7248071429e+03,1.e-10)
  PUSH_MPS("sctap3",true,1481,2480,1.4240000000e+03,1.e-10)
  PUSH_MPS("seba",true,516,1028,1.5711600000e+04,1.e-10)
  PUSH_MPS("share1b",true,118,225,-7.6589318579e+04,1.e-10)
  PUSH_MPS("share2b",true,97,79,-4.1573224074e+02,1.e-10)
  PUSH_MPS("shell",true,537,1775,1.2088253460e+09,1.e-10)
  PUSH_MPS("ship04l",true,403,2118,1.7933245380e+06,1.e-10)
  PUSH_MPS("ship04s",true,403,1458,1.7987147004e+06,1.e-10)
  PUSH_MPS("ship08l",true,779,4283,1.9090552114e+06,1.e-10)
  PUSH_MPS("ship08s",true,779,2387,1.9200982105e+06,1.e-10)
  PUSH_MPS("ship12l",true,1152,5427,1.4701879193e+06,1.e-10)
  PUSH_MPS("ship12s",true,1152,2763,1.4892361344e+06,1.e-10)
  PUSH_MPS("sierra",true,1228,2036,1.5394362184e+07,1.e-10)
  PUSH_MPS("stair",true,357,467,-2.5126695119e+02,1.e-10)
  PUSH_MPS("standata",true,360,1075,1.2576995000e+03,1.e-10)
  // GUB PUSH_MPS("standgub",true,362,1184,1257.6995,1.e-10) 
  PUSH_MPS("standmps",true,468,1075,1.4060175000E+03,1.e-10) 
  PUSH_MPS("stocfor1",true,118,111,-4.1131976219E+04,1.e-10)
  PUSH_MPS("stocfor2",true,2158,2031,-3.9024408538e+04,1.e-10)
  // ?? PUSH_MPS("stocfor3",true,16676,15695,-3.9976661576e+04,1.e-10)
  // ?? PUSH_MPS("truss",true,1001,8806,4.5881584719e+05,1.e-10)
  PUSH_MPS("tuff",true,334,587,2.9214776509e-01,1.e-10)
  PUSH_MPS("vtpbase",true,199,203,1.2983146246e+05,1.e-10)
  PUSH_MPS("wood1p",true,245,2594,1.4429024116e+00,1.e-10)
  PUSH_MPS("woodw",true,1099,8405,1.3044763331E+00,1.e-10)

#undef PUSH_MPS

/*
  Create a vector of solver interfaces that we can use to run the test
  problems. The strategy is to create a fresh clone of the `empty' solvers
  from vecEmptySiP for each problem, then proceed in stages: read the MPS
  file, solve the problem, check the solution. If there are multiple
  solvers in vecSiP, the results of each solver are compared with its
  neighbors in the vector.
*/
  std::vector<OsiSolverInterface*> vecSiP(vecEmptySiP.size()) ;

  // Create vector to store a name for each solver interface
  // and a count on the number of problems the solver intface solved.
  std::vector<std::string> siName;
  std::vector<int> numProbSolved;
  for ( i=0; i<vecSiP.size(); i++ ) {
    siName.push_back("");
    numProbSolved.push_back(0);
  }

/*
  Open the main loop to step through the MPS problems.
*/
  for (m = 0 ; m < mpsName.size() ; m++)
  { std::cerr << "  processing mps file: " << mpsName[m] 
      << " (" << m+1 << " out of " << mpsName.size() << ")" << std::endl ;
/*
  Stage 1: Read the MPS file into each solver interface.

  Fill vecSiP with fresh clones of the solvers and read in the MPS file. As
  a basic check, make sure the size of the constraint matrix is correct.
*/
    for (i = vecSiP.size()-1 ; i >= 0 ; --i)
    { vecSiP[i] = vecEmptySiP[i]->clone() ;
      
      std::string fn = mpsDir+mpsName[m] ;
      vecSiP[i]->readMps(fn.c_str(),"mps") ;
      
      if (min[m])
	vecSiP[i]->setObjSense(1.0) ;
      else
	vecSiP[i]->setObjSense(-1.0) ;
      
      int nr = vecSiP[i]->getNumRows() ;
      int nc = vecSiP[i]->getNumCols() ;
      assert(nr == nRows[m]-1) ;
      assert(nc == nCols[m]) ; } 
/*
  If we have multiple solvers, compare the representations.
*/
    for (i = vecSiP.size()-1 ; i > 0 ; --i)
    { OsiPackedVector vim1,vi ;
      
      // Compare col lowerbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getColLower(),vecSiP[i  ]->getColLower(),
        vecSiP[i  ]->getNumCols() )
        ) ;
      
      // Compare col upperbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getColUpper(),vecSiP[i  ]->getColUpper(),
        vecSiP[i  ]->getNumCols() )
        ) ;
      
      // Compare row lowerbounds
      assert(
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowLower(),vecSiP[i  ]->getRowLower(),
        vecSiP[i  ]->getNumRows() )
        ) ;
      
      // Compare row upperbounds
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowUpper(),vecSiP[i  ]->getRowUpper(),
        vecSiP[i  ]->getNumRows() )
        ) ;
            
      // Compare row sense
      {
        const char * rsm1 = vecSiP[i-1]->getRowSense() ;
        const char * rs   = vecSiP[i  ]->getRowSense() ;
        int nr = vecSiP[i]->getNumRows() ;
        int r ;
        for (r = 0 ; r < nr ; r++) assert (rsm1[r] == rs[r]) ;
      }
      
      // Compare rhs
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRightHandSide(),vecSiP[i  ]->getRightHandSide(),
        vecSiP[i  ]->getNumRows() )
        ) ;
      
      // Compare range
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getRowRange(),vecSiP[i  ]->getRowRange(),
        vecSiP[i  ]->getNumRows() )
        ) ;
      
      // Compare objective sense
      assert( vecSiP[i-1]->getObjSense() == vecSiP[i  ]->getObjSense() ) ;
           
      // Compare objective coefficients
      assert( 
        equivalentVectors(vecSiP[i-1],vecSiP[i], 1.e-10,
        vecSiP[i-1]->getObjCoefficients(),vecSiP[i  ]->getObjCoefficients(),
        vecSiP[i  ]->getNumCols() )
        ) ;

      // Compare number of elements
      assert( vecSiP[i-1]->getNumElements() == vecSiP[i]->getNumElements() ) ;
      
      // Compare constraint matrix
      { const OsiPackedMatrix * rmm1=vecSiP[i-1]->getMatrixByRow() ;
        const OsiPackedMatrix * rm  =vecSiP[i  ]->getMatrixByRow() ;
        assert( rmm1->isEquivalent(*rm) ) ;
        
        const OsiPackedMatrix * cmm1=vecSiP[i-1]->getMatrixByCol() ;
        const OsiPackedMatrix * cm  =vecSiP[i  ]->getMatrixByCol() ;
        assert( cmm1->isEquivalent(*cm) ) ; } }
/*
  If we have multiple solvers, compare the variable type information
*/
    for (i = vecSiP.size()-1 ; i > 0 ; --i)
    { OsiPackedVector vim1,vi ;
      int c ;
      
      { OsiVectorInt sm1 = vecSiP[i-1]->getFractionalIndices() ;
        OsiVectorInt s   = vecSiP[i  ]->getFractionalIndices() ;
        assert( sm1.size() == s.size() ) ;
        for (c = s.size()-1 ; c >= 0 ; --c) assert( sm1[c] == s[c] ) ; }

      { int nc = vecSiP[i]->getNumCols() ;
        for (c = 0 ; c < nc ; c++)
	{ assert(
	    vecSiP[i-1]->isContinuous(c) == vecSiP[i]->isContinuous(c)
	    ) ;
          assert(
	    vecSiP[i-1]->isBinary(c) == vecSiP[i]->isBinary(c)
	    ) ;
          assert(
	    vecSiP[i-1]->isIntegerNonBinary(c) ==
	    vecSiP[i  ]->isIntegerNonBinary(c)
	    ) ;
          assert(
	    vecSiP[i-1]->isFreeBinary(c) == vecSiP[i]->isFreeBinary(c)
	    ) ;
          assert(
	    vecSiP[i-1]->isInteger(c) == vecSiP[i]->isInteger(c)
	    ) ; } } }
/*
  Stage 2: Call each solver to solve the problem.

  We call each solver, then check the return code and objective.

  Note that the volume solver can't handle the Netlib cases. The strategy is
  to require that it be the last solver in vecSiP and then break out of the
  loop. This ensures that all previous solvers are run and compared to one
  another.
*/

    for (i = 0 ; i < static_cast<int>(vecSiP.size()) ; ++i)
    {
#     ifdef COIN_USE_VOL
      { 
        OsiVolSolverInterface * si =
          dynamic_cast<OsiVolSolverInterface *>(vecSiP[i]) ;
        if (si != NULL )  { 
          siName[i]="OsiVolSolverInterface";
          // VOL does not solve netlib cases so don't bother trying to solve
          break ; 
        }
      }
#     endif
#     ifdef COIN_USE_DYLP
      { 
        OsiDylpSolverInterface * si =
          dynamic_cast<OsiDylpSolverInterface *>(vecSiP[i]) ;
        if (si != NULL )  {    
          siName[i]="OsiDylpSolverInterface";
          // Is this an MPS file that OslDylpSolverInterface handles
          if ( mpsName[m]=="cycle" ||
               mpsName[m]=="d6cube" ||
               mpsName[m]=="fit1d" || 
               mpsName[m]=="grow15" || 
               mpsName[m]=="grow22" || 
               mpsName[m]=="maros" || 
               mpsName[m]=="pilot" || 
               mpsName[m]=="pilot4" || 
               mpsName[m]=="pilotnov" || 
               mpsName[m]=="wood1p" )break ; 
        }
      }
#     endif
#     ifdef COIN_USE_XPR
      { 
        OsiXprSolverInterface * si =
          dynamic_cast<OsiXprSolverInterface *>(vecSiP[i]) ;
        if (si != NULL )  {    
          siName[i]="OsiXprSolverInterface";
        }
      }
#     endif
#     ifdef COIN_USE_CPX
      { 
        OsiCpxSolverInterface * si =
          dynamic_cast<OsiCpxSolverInterface *>(vecSiP[i]) ;
        if (si != NULL )  {    
          siName[i]="OsiCpxSolverInterface";
        }
      }
#     endif
#     ifdef COIN_USE_OSL
      { 
        OsiOslSolverInterface * si =
          dynamic_cast<OsiOslSolverInterface *>(vecSiP[i]) ;
        if (si != NULL )  {    
          siName[i]="OsiOslSolverInterface";
        }
      }
#     endif
      
      vecSiP[i]->initialSolve() ;
      
      if (vecSiP[i]->isProvenOptimal()) { 
        double soln = vecSiP[i]->getObjValue();       
        OsiRelFltEq eq(objValueTol[m]) ;
        if (eq(soln,objValue[m])) { 
          //std::cerr << soln << " = " << objValue[m] << " ; ok." <<std::endl; 
          numProbSolved[i]++; 
        }
        else  { 
          std::cerr <<siName[i] <<" " <<soln << " != " <<objValue[m] << "; error=" ;
          std::cerr <<fabs(objValue[m] - soln) <<std::endl; 
        }
      }
      else
        if (vecSiP[i]->isProvenPrimalInfeasible()) { 
          std::cerr << "error; primal infeasible" << std::endl ; 
        }
        else
          if (vecSiP[i]->isProvenDualInfeasible()) { 
            std::cerr << "error; dual infeasible" << std::endl ; 
          }
          else
            if (vecSiP[i]->isIterationLimitReached()) { 
              std::cerr << "error; iteration limit" << std::endl ; 
            }
            else
              if (vecSiP[i]->isAbandoned())  { 
                std::cerr << "error; abandoned" << std::endl ; 
              }
              else  { 
                std::cerr << "error; unknown" << std::endl ; 
              } 
    }
    /*
    Delete the used solver interfaces so we can reload fresh clones for the
              next problem.
    */
    for (i = vecSiP.size()-1 ; i >= 0 ; --i) delete vecSiP[i] ;
  }

  for ( i=0; i<siName.size(); i++ ) {
    std::cerr 
      <<siName[i] 
      <<" solved " 
      <<numProbSolved[i]
      <<" out of "
      <<objValue.size()
      <<std::endl;
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
    assert( si->getColSolution()==NULL );
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

    // Test getting and setting of objective offset
    double objOffset;
    assert( exmip1Si->getDblParam(OsiObjOffset,objOffset) );
    assert( eq( objOffset, 0.0 ) );
    assert( exmip1Si->setDblParam(OsiObjOffset, 3.21) );
    assert( exmip1Si->getDblParam(OsiObjOffset,objOffset) );
    assert( eq( objOffset, 3.21 ) );
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
      {
        double objOffset;
        base->getDblParam(OsiObjOffset,objOffset);
        si1->setDblParam(OsiObjOffset,objOffset);
        si2->setDblParam(OsiObjOffset,objOffset);
        si3->setDblParam(OsiObjOffset,objOffset);
        si4->setDblParam(OsiObjOffset,objOffset);
        si5->setDblParam(OsiObjOffset,objOffset);
        si6->setDblParam(OsiObjOffset,objOffset);
        si7->setDblParam(OsiObjOffset,objOffset);
        si8->setDblParam(OsiObjOffset,objOffset);
      }
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

