// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <string>
#include <iostream>

#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"

#include "OsiRowCutDebugger.hpp"

/* If we are on the path to the optimal integer solution then
   check if any generated cuts cut off the optimal solution!
   
   If so then print offending cuts and return non-zero code
*/


int OsiRowCutDebugger::validateCuts(const OsiCuts & cs, 
				  int first, int last) const
{
  int nbad=0; 
  int i;
  const double epsilon=1.0e-8;
  const int nRowCuts = CoinMin(cs.sizeRowCuts(),last);
  
  for (i=first; i<nRowCuts; i++){
    
    OsiRowCut rcut = cs.rowCut(i);
    CoinPackedVector rpv = rcut.row();
    const int n = rpv.getNumElements();
    const int * indices = rpv.getIndices();
    const double * elements = rpv.getElements();
    int k;
    double lb=rcut.lb();
    double ub=rcut.ub();
    
    double sum=0.0;
    
    for (k=0; k<n; k++){
      int column=indices[k];
      sum += optimalSolution_[column]*elements[k];
    }
    // is it violated
    if (sum >ub + epsilon ||sum < lb - epsilon) {
      double violation=CoinMax(sum-ub,lb-sum);
      std::cout<<"Cut "<<i<<" with "<<n
	  <<" coefficients, cuts off optimal solutions by "<<violation
          <<", lo="<<lb<<", ub="<<ub<<std::endl;
      for (k=0; k<n; k++){
	int column=indices[k];
        std::cout<<"( "<<column<<" , "<<elements[k]<<" ) ";
      }
      std::cout <<std::endl;
      nbad++;
    }
  }
  return nbad;
}


/* If we are on the path to the optimal integer solution then
   check out if generated cut cuts off the optimal solution!
   
   If so then print offending cut and return non-zero code
*/

bool OsiRowCutDebugger::invalidCut(const OsiRowCut & rcut) const 
{
  bool bad=false; 
  const double epsilon=1.0e-8;
  
  CoinPackedVector rpv = rcut.row();
  const int n = rpv.getNumElements();
  const int * indices = rpv.getIndices();
  const double * elements = rpv.getElements();
  int k;
  
  double lb=rcut.lb();
  double ub=rcut.ub();
  double sum=0.0;
  
  for (k=0; k<n; k++){
    int column=indices[k];
    sum += optimalSolution_[column]*elements[k];
  }
  // is it violated
  if (sum >ub + epsilon ||sum < lb - epsilon) {
    double violation=CoinMax(sum-ub,lb-sum);
    std::cout<<"Cut with "<<n
	<<" coefficients, cuts off optimal solutions by "<<violation
        <<", lo="<<lb<<", ub="<<ub<<std::endl;
    for (k=0; k<n; k++){
      int column=indices[k];
      std::cout<<"( "<<column<<" , "<<elements[k]<<" ) ";
    }
    std::cout<<std::endl;
    bad=true;
  }
  return bad;
}

// Returns true if still on optimal path, false otherwise
bool OsiRowCutDebugger::onOptimalPath(const OsiSolverInterface & si) const
{
  if (integerVariable_) {
    int nCols=si.getNumCols(); 
    if (nCols!=numberColumns_)
      return false; // check user has not modified problem
    int i;
    const double * collower = si.getColLower();
    const double * colupper = si.getColUpper();
    bool onOptimalPath=true;
    for (i=0;i<numberColumns_;i++) {
      if (integerVariable_[i]) {
	// value of integer variable in solution
	double value=optimalSolution_[i]; 
	if (value>colupper[i]+1.0e-3 || value<collower[i]-1.0e-3) {
	  onOptimalPath=false;
	  break;
	}
      }
    }
    return onOptimalPath;
  } else {
    // no information
    return false;
  }
}
// Returns true if debugger is active 
bool OsiRowCutDebugger::active() const
{
  return (integerVariable_!=NULL);
}

// Activate using name of model
// returns whether debug activated
bool OsiRowCutDebugger::activate( const OsiSolverInterface & si, 
				   const char * model)
{
  int i;
  //get rid of any arrays
  delete [] integerVariable_;
  delete [] optimalSolution_;
  numberColumns_ = 0;
  int expectedNumberColumns = 0;


  enum {undefined, pure0_1, continuousWith0_1, generalMip } probType;


  // Convert input parameter model to be lowercase and 
  // only consider charcters between '/' and '.'
  std::string modelL; //name in lowercase 
  int iput=0;
  for (i=0;i<(int) strlen(model);i++) {
    char value=tolower(model[i]);
    if (value=='/') {
      iput=0;
      modelL.erase();
    } else if (value=='.') {
      break;
    } else {
      modelL.append(1,value);
    }
  }


  CoinPackedVector intSoln;
  probType = undefined;

  //--------------------------------------------------------
  //
  // Define additional problems by adding it as an additional
  // "else if ( modelL == '???' ) { ... }" 
  // stanza below.
  //
  // Assign values to probType and intSoln.
  //
  // probType - pure0_1, continuousWith0_1, or generalMip
  // 
  // intSoln -
  //    when probType is pure0_1
  //       intSoln contains the indices of the variables
  //       at 1 in the optimal solution
  //    when probType is continuousWith0_1
  //       intSoln contains the indices of integer
  //       variables at one in the optimal solution
  //    when probType is generalMip
  //       intSoln contains the the indices of the integer
  //       variables and their value in the optimal solution
  //--------------------------------------------------------

  // exmip1
  if ( modelL == "exmip1" ) {
    probType=continuousWith0_1;
    intSoln.insert(2,1.);
    intSoln.insert(3,1.);
    expectedNumberColumns=8;
  }

  // p0033
  else if ( modelL == "p0033" ) {
    probType=pure0_1;
    int intIndicesAt1[]={ 0,6,7,9,13,17,18,22,24,25,26,27,28,29 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=33;
  }

  // flugpl
  else if ( modelL == "flugpl" ) {
    probType=generalMip;
    int intIndicesV[] = { 1 , 3 , 4 , 6 , 7 , 9 ,10 ,12 ,13 ,15 };
    double intSolnV[] = { 6.,60., 6.,60.,16.,70., 7.,70.,12.,75.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=18;
  }

  // enigma
  else if ( modelL == "enigma" ) {
    probType=pure0_1;
    int intIndicesAt1[]={ 0,18,25,36,44,59,61,77,82,93 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=100;
  }

  // mod011
  else if ( modelL == "mod011" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={ 10,29,32,40,58,77,80,88 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=10958;
  }

  // mas76
  else if ( modelL == "mas76" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={ 4,11,13,18,42,46,48,52,85,93,114,119,123,128,147}; 
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=151;
  }

  // ltw3
  else if ( modelL == "ltw3" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={ 20,23,24,26,32,33,40,47 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=48;
  }

  // mod008
  else if ( modelL == "mod008" ) {
    probType=pure0_1;
    int intIndicesAt1[]={1,59,83,116,123};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=319;
  }

  // mod010
  else if ( modelL == "mod010" ) {
    probType=pure0_1;
    int intIndicesAt1[]={2,9,16,22,26,50,65,68,82,86,102,145,
	149,158,181,191,266,296,376,479,555,625,725,851,981,
	1030,1095,1260,1321,1339,1443,1459,1568,1602,1780,1856,
	1951,2332,2352,2380,2471,2555,2577,2610,2646,2647};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=2655;
  }

  // modglob
  else if ( modelL == "modglob" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={204,206,208,212,216,218,220,222,230,232,
	234,236,244,248,250,254,256,258,260,262,264,266,268,274,
	278,282,284,286,288};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=422;
  }

  // p0201
  else if ( modelL == "p0201" ) {
    probType=pure0_1;
    int intIndicesAt1[]={8,10,21,38,39,56,60,74,79,92,94,110,111,
	128,132,146,151,164,166,182,183,200};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=201;
  }

  // p0282
  else if ( modelL == "p0282" ) {
    probType=pure0_1;
    int intIndicesAt1[]={3,11,91,101,103,117,155,169,191,199,215,
	223,225,237,240,242,243,244,246,248,251,254,256,257,260,
	262,263,273,275,276,277,280,281};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=282;
  }

  // p0548
  else if ( modelL == "p0548" ) {
    probType=pure0_1;
    int intIndicesAt1[]={2,3,13,14,17,23,24,43,44,47,61,62,74,75,
	81,82,92,93,96,98,105,120,126,129,140,141,153,154,161,162,
	165,177,182,184,189,192,193,194,199,200,209,214,215,218,222,
	226,234,239,247,256,257,260,274,286,301,305,306,314,317,318,
	327,330,332,334,336,340,347,349,354,358,368,369,379,380,385,
	388,389,390,393,394,397,401,402,406,407,417,419,420,423,427,
	428,430,437,439,444,446,447,450,451,452,472,476,477,480,488,
	491,494,500,503,508,509,510,511,512,515,517,518,519,521,522,
	523,525,526,527,528,529,530,531,532,533,536,537,538,539,541,
	542,545,547};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=548;
  }

  // p2756
  else if ( modelL == "p2756" ) {
    probType=pure0_1;
    int intIndicesAt1[]={7,25,50,63,69,71,81,124,164,208,210,212,214,
	220,266,268,285,299,301,322,362,399,455,464,468,475,518,574,
	588,590,612,632,652,679,751,767,794,819,838,844,892,894,913,
	919,954,966,996,998,1021,1027,1044,1188,1230,1248,1315,1348,
	1366,1367,1420,1436,1473,1507,1509,1521,1555,1558,1607,1659,
	1715,1746,1761,1789,1800,1844,1885,1913,1916,1931,1992,2002,
	2050,2091,2155,2158,2159,2197,2198,2238,2264,2292,2318,2481,
	2496,2497,2522,2531,2573,2583,2587,2588,2596,2635,2637,2639,
	2643,2645,2651,2653,2672,2675,2680,2683,2708,2727,2730,2751};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=2756;
  }

  // nw04
  else if ( modelL == "nw04" ) {
    probType=pure0_1;
    int intIndicesAt1[]={
      231 ,1792 ,1980 ,7548 ,21051 ,28514 ,53087 ,53382 ,76917 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=87482;
  }

  // bell3a
  else if ( modelL == "bell3a" ) {
    probType=generalMip;
    int intIndicesV[]={61,62,65,66,67,68,69,70};
    double intSolnV[] = {4.,21.,4.,4.,6.,1.,25.,8.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=133;
  }

  // 10teams
  else if ( modelL == "10teams" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={236,298,339,379,443,462,520,576,616,646,690,
	749,778,850,878,918,986,996,1065,1102,1164,1177,1232,1281,1338,
	1358,1421,1474,1522,1533,1607,1621,1708,1714,1775,1835,1887,
	1892,1945,1989};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=2025;
  }

  // rentacar
  else if ( modelL == "rentacar" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      9502 ,9505 ,9507 ,9511 ,9512 ,9513 ,9514 ,9515 ,9516 ,9521 ,
      9522 ,9526 ,9534 ,9535 ,9536 ,9537 ,9542 ,9543 ,9544 ,9548 ,
      9550 ,9554 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=9557;
  }

  // qiu
  else if ( modelL == "qiu" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      0 ,5 ,8 ,9 ,11 ,13 ,16 ,17 ,19 ,20 ,
      24 ,28 ,32 ,33 ,35 ,37 ,40 ,47 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=840;
  }

  // pk1
  else if ( modelL == "pk1" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      1 ,4 ,5 ,6 ,7 ,11 ,13 ,16 ,17 ,23 ,
      24 ,27 ,28 ,34 ,35 ,37 ,43 ,44 ,45 ,46 ,
      47 ,51 ,52 ,54 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=86;
  }

  // pp08a
  else if ( modelL == "pp08a" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      177 ,179 ,181 ,183 ,185 ,190 ,193 ,195 ,197 ,199 ,
      202 ,204 ,206 ,208 ,216 ,220 ,222 ,229 ,235 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=240;
  }

  // pp08aCUTS
  else if ( modelL == "pp08acuts" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      177 ,179 ,181 ,183 ,185 ,190 ,193 ,195 ,197 ,199 ,
      202 ,204 ,206 ,208 ,216 ,220 ,222 ,229 ,235 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=240;
  }

  // danoint
  else if ( modelL == "danoint" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={3,5,8,11,15,21,24,25,31,34,37,42,46,48,51,56};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=521;
  }

  // dcmulti
  else if ( modelL == "dcmulti" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={2,3,11,14,15,16,21,24,28,34,35,36,39,40,41,42,
	45,52,53,60,61,64,65,66,67};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=548;
  }

  // egout
  else if ( modelL == "egout" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={0,3,5,7,8,9,11,12,13,15,16,17,18,20,21,22,
	23,24,25,26,27,28,29,32,34,36,37,38,39,40,42,43,44,45,46,47,
	48,49,52,53,54};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=141;
  }

  // fixnet6
  else if ( modelL == "fixnet6" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={1,16,23,31,37,51,64,179,200,220,243,287,
	375,413,423,533,537,574,688,690,693,712,753,773,778,783,847};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=878;
  }

  // khb05250
  else if ( modelL == "khb05250" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={1,3,8,11,12,15,16,17,18,21,22,23};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=1350;
  }

  // lseu
  else if ( modelL == "lseu" ) {
    probType=pure0_1;
    int intIndicesAt1[]={0,1,6,13,26,33,38,43,50,52,63,65,85};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=89;
  }

  // air03
  else if ( modelL == "air03" ) {
    probType=pure0_1;
    int intIndicesAt1[]={
      1, 3, 5, 13, 14, 28, 38, 49, 75, 76, 
      151, 185, 186, 271, 370, 466, 570, 614, 732, 819, 
      1151, 1257, 1490, 2303, 2524, 3301, 3616, 4129, 4390, 4712, 
      5013, 5457, 5673, 6436, 7623, 8122, 8929, 10689, 10694, 10741, 
      10751
    };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=10757;
  }

  // air04
  else if ( modelL == "air04" ) {
    probType=pure0_1;
    int intIndicesAt1[]={
      0, 1, 3, 4, 5, 6, 7, 9, 11, 12, 
      13, 17, 19, 20, 21, 25, 26, 27, 28, 29, 
      32, 35, 36, 39, 40, 42, 44, 45, 47, 48, 
      49, 50, 51, 52, 53, 56, 57, 58, 60, 63, 
      64, 66, 67, 68, 73, 74, 80, 81, 83, 85, 
      87, 92, 93, 94, 95, 99, 101, 102, 105, 472, 
      616, 680, 902, 1432, 1466, 1827, 2389, 2535, 2551, 2883, 
      3202, 3215, 3432, 3438, 3505, 3517, 3586, 3811, 3904, 4092, 
      4685, 4700, 4834, 4847, 4892, 5189, 5211, 5394, 5878, 6045, 
      6143, 6493, 6988, 7511, 7664, 7730, 7910, 8041, 8350, 8615, 
      8635, 8670
    };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=8904;
  }

  // air05
  else if ( modelL == "air05" ) {
    probType=pure0_1;
    int intIndicesAt1[]={
      2, 4, 5, 6, 7, 8, 9, 10, 14, 15, 
      19, 20, 25, 34, 35, 37, 39, 40, 41, 42, 
      43, 44, 45, 47, 48, 50, 52, 55, 57, 58, 
      66, 72, 105, 218, 254, 293, 381, 695, 1091, 1209, 
      1294, 1323, 1348, 1580, 1769, 2067, 2156, 2162, 2714, 2732, 
      3113, 3131, 3145, 3323, 3398, 3520, 3579, 4295, 5025, 5175, 
      5317, 5340, 6324, 6504, 6645, 6809
    };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=7195;
  }

  // seymour
  else if ( modelL == "seymour" ) {
    probType=pure0_1;
    int intIndicesAt1[]=
      {
	1, 2, 3, 5, 6, 7, 9, 11, 12, 16, 
	18, 22, 23, 25, 27, 31, 32, 34, 35, 36, 
	38, 39, 40, 42, 44, 45, 46, 49, 50, 51, 
	52, 54, 55, 56, 58, 61, 63, 65, 67, 68, 
	69, 70, 71, 75, 79, 81, 82, 84, 85, 86, 
	87, 88, 89, 91, 93, 95, 97, 98, 99, 100, 
	101, 102, 103, 106, 108, 112, 116, 118, 119, 120, 
	122, 123, 124, 125, 126, 129, 130, 132, 135, 137, 
	140, 141, 142, 143, 144, 148, 150, 151, 154, 156, 
	159, 160, 162, 163, 164, 165, 167, 169, 170, 174, 
	177, 178, 180, 181, 182, 183, 188, 189, 192, 194, 
	200, 201, 202, 203, 204, 211, 214, 218, 226, 227, 
	228, 231, 232, 237, 240, 242, 244, 247, 248, 249, 
	251, 253, 256, 257, 259, 261, 264, 265, 266, 268, 
	270, 272, 278, 280, 284, 286, 288, 289, 291, 292, 
	296, 299, 302, 305, 307, 308, 311, 312, 313, 314, 
	315, 316, 317, 319, 321, 325, 328, 332, 334, 335, 
	337, 338, 339, 340, 343, 346, 355, 357, 358, 365, 
	369, 372, 373, 374, 375, 376, 378, 381, 383, 386, 
	392, 396, 399, 402, 403, 412, 416, 419, 424, 425, 
	426, 427, 430, 431, 432, 436, 437, 438, 440, 441, 
	443, 450, 451, 452, 453, 456, 460, 461, 462, 467, 
	469, 475, 476, 477, 478, 479, 485, 486, 489, 491, 
	493, 498, 500, 501, 508, 513, 515, 516, 518, 519, 
	520, 524, 527, 541, 545, 547, 548, 559, 562, 563, 
	564, 566, 567, 570, 572, 575, 576, 582, 583, 587, 
	589, 595, 599, 602, 610, 611, 615, 622, 631, 646, 
	647, 649, 652, 658, 662, 665, 667, 671, 676, 679, 
	683, 685, 686, 688, 689, 691, 699, 705, 709, 711, 
	712, 716, 721, 722, 724, 726, 729, 732, 738, 739, 
	741, 745, 746, 747, 749, 752, 757, 765, 767, 768, 
	775, 779, 780, 791, 796, 798, 808, 809, 812, 813, 
	817, 819, 824, 825, 837, 839, 849, 851, 852, 857, 
	865, 874, 883, 885, 890, 897, 902, 907, 913, 915, 
	923, 924, 927, 931, 933, 936, 938, 941, 945, 949, 
	961, 970, 971, 978, 984, 985, 995, 997, 999, 1001, 
	1010, 1011, 1012, 1025, 1027, 1035, 1043, 1055, 1056, 1065, 
	1077, 1089, 1091, 1096, 1100, 1104, 1112, 1126, 1130, 1131, 
	1132, 1134, 1136, 1143, 1149, 1162, 1163, 1164, 1183, 1184, 
	1191, 1200, 1201, 1209, 1215, 1220, 1226, 1228, 1229, 1233, 
	1241, 1243, 1244, 1258, 1277, 1279, 1285, 1291, 1300, 1303, 
	1306, 1311, 1320, 1323, 1333, 1344, 1348, 1349, 1351, 1356, 
	1363, 1364, 1365, 1366};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=1372;
  }

  // stein27
  else if ( modelL == "stein27" ) {
    probType=pure0_1;
    int intIndicesAt1[]={0,1,3,4,5,6,7,8,9,11,13,16,17,19,21,22,25,26};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=27;
  }

  // stein45
  else if ( modelL == "stein45" ) {
    probType=pure0_1;
    int intIndicesAt1[]={0,1,4,5,6,7,8,9,10,11,14,17,18,19,21,23,24,25,26,28,
    31,32,33,36,37,39,40,42,43,44};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=45;
  }

  // misc03
  else if ( modelL == "misc03" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={4,40,62,75,99,114,127,134,147,148,150,
	152,154,155,157};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=160;
  }

  // misc06
  else if ( modelL == "misc06" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={
      1557 ,1560 ,1561 ,1580 ,1585 ,1588 ,1589 ,1614 ,1615 ,1616 ,
      1617 ,1626 ,1630 ,1631 ,1642 ,1643 ,1644 ,1645 ,1650 ,1654 ,
      1658 ,1659 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=1808;
  }

  // misc07
  else if ( modelL == "misc07" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={21,27,57,103,118,148,185,195,205,209,243,
	245,247,249,251,253,255,257};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=260;
  }

  // rgn
  else if ( modelL == "rgn" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={16 ,49 ,72 ,92 };
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=180;
  }

  // mitre
  else if ( modelL == "mitre" ) {
    probType=pure0_1;
    int intIndicesAt1[]={
      4,37,67,93,124,154,177,209,240,255,287,319,340,
      372,403,425,455,486,516,547,579,596,628,661,676,
      713,744,758,795,825,851,881,910,933,963,993,1021,
      1052,1082,1111,1141,1172,1182,1212,1242,1272,1303,
      1332,1351,1382,1414,1445,1478,1508,1516,1546,1576,
      1601,1632,1662,1693,1716,1749,1781,1795,1828,1860,
      1876,1909,1940,1962,1994,2027,2058,2091,2122,2128,
      2161,2192,2226,2261,2290,2304,2339,2369,2393,2426,
      2457,2465,2500,2529,2555,2590,2619,2633,2665,2696,
      2728,2760,2792,2808,2838,2871,2896,2928,2960,2981,
      3014,3045,3065,3098,3127,3139,3170,3200,3227,3260,
      3292,3310,3345,3375,3404,3437,3467,3482,3513,3543,
      3558,3593,3623,3653,3686,3717,3730,3762,3794,3814,
      3845,3877,3901,3936,3966,3988,4019,4049,4063,4096,
      4126,4153,4186,4216,4245,4276,4306,4318,4350,4383,
      4402,4435,4464,4486,4519,4550,4578,4611,4641,4663,
      4695,4726,4738,4768,4799,4830,4863,4892,4919,4950,
      4979,4991,5024,5054,5074,5107,5137,5165,5198,5228,
      5244,5275,5307,5325,5355,5384,5406,5436,5469,5508,
      5538,5568,5585,5615,5646,5675,5705,5734,5745,5774,
      5804,5836,5865,5895,5924,5954,5987,6001,6033,6064,
      6096,6126,6155,6172,6202,6232,6250,6280,6309,6328,
      6361,6392,6420,6450,6482,6500,6531,6561,6598,6629,
      6639,6669,6699,6731,6762,6784,6814,6844,6861,6894,
      6924,6955,6988,7018,7042,7075,7105,7116,7149,7179,
      7196,7229,7258,7282,7312,7345,7376,7409,7438,7457,
      7487,7520,7534,7563,7593,7624,7662,7692,7701,7738,
      7769,7794,7827,7857,7872,7904,7935,7960,7990,8022,
      8038,8071,8101,8137,8167,8199,8207,8240,8269,8301,
      8334,8363,8387,8420,8450,8470,8502,8534,8550,8580,
      8610,8639,8669,8699,8709,8741,8772,8803,8834,8867,
      8883,8912,8942,8973,9002,9032,9061,9094,9124,9128,
      9159,9201,9232,9251,9280,9310,9333,9338,9405,9419,
      9423,9428,9465,9472,9482,9526,9639,9644,9666,9673,
      9729,9746,9751,9819,9832,9833,9894,9911,9934,9990,
      10007,10012,10083,10090,10095,10137,10176,10177,
      10271,10279,10280,10288,10292,10298,10299,10319,
      10351,10490,10505,10553,10571,10579,10600,10612,
      10683,10688};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=10724;
  }

  // gen
  else if ( modelL == "gen" ) {
    probType=generalMip;
    int intIndicesV[]={15,34,35,36,37,38,39,40,41,42,43,44,45,57,58,
	59,60,61,62,63,64,65,66,67,68,69,84,85,86,87,88,89,90,91,92,
	93,107,108,109,110,111,112,113,114,120,121,122,123,124,125,
	126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,
		       141,142,143,432,433,434,435,436};
    double intSolnV[] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
	1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
	1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
	1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,23.,12.,11.,14.,16.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=870;
  }

  // noswot
  else if ( modelL == "noswot" ) {
    probType=generalMip;
    int intIndicesV[]={};
    double intSolnV[]={};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=128;
  }


  // qnet1
  else if ( modelL == "qnet1" ) {
    probType=generalMip;
    int intIndicesV[]={
      61 ,69 ,79 ,81 ,101 ,104 ,111 ,115 ,232 ,265 ,
      267 ,268 ,269 ,285 ,398 ,412 ,546 ,547 ,556 ,642 ,
      709 ,718 ,721 ,741 ,760 ,1073 ,1077 ,1084 ,1097 ,1100 ,
      1104 ,1106 ,1109 ,1248 ,1259 ,1260 ,1263 ,1265 ,1273 ,1274 ,
      1276 ,1286 ,1291 ,1302 ,1306 ,1307 ,1316 ,1350 ,1351 ,1363 ,
      1366 ,1368 ,1371 ,1372 ,1380 ,1381 ,1385 ,1409 };
    double intSolnV[]={
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,2. ,1. ,1. ,1. ,1. ,3. ,
      1. ,2. ,1. ,1. ,1. ,1. ,6. ,3. ,1. ,1. ,
      5. ,2. ,1. ,2. ,2. ,1. ,2. ,1. };
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=1541;
  }

  // qnet1_o (? marginally different from qnet1)
  else if ( modelL == "qnet1_o" ) {
    probType=generalMip;
    int intIndicesV[]={
      61 ,69 ,79 ,81 ,101 ,106 ,111 ,114 ,115 ,232 ,
      266 ,267 ,268 ,269 ,277 ,285 ,398 ,412 ,546 ,547 ,
      556 ,642 ,709 ,718 ,721 ,741 ,760 ,1073 ,1077 ,1084 ,
      1097 ,1100 ,1104 ,1106 ,1109 ,1248 ,1259 ,1260 ,1263 ,1265 ,
      1273 ,1274 ,1276 ,1286 ,1291 ,1302 ,1306 ,1307 ,1316 ,1350 ,
      1351 ,1363 ,1366 ,1368 ,1371 ,1372 ,1380 ,1381 ,1385 ,1409 };
    double intSolnV[]={
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,1. ,
      1. ,1. ,1. ,1. ,1. ,1. ,2. ,1. ,1. ,1. ,
      1. ,3. ,1. ,2. ,1. ,1. ,1. ,1. ,6. ,3. ,
      1. ,1. ,5. ,2. ,1. ,2. ,2. ,1. ,2. ,1. };
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=1541;
  }

  // gt2
  else if ( modelL == "gt2" ) {
    probType=generalMip;
    int intIndicesV[]={82,85,88,92,94,95,102,103,117,121,122,128,
		       141,146,151,152,165,166,176,179};
    double intSolnV[] = {1.,3.,1.,5.,2.,1.,1.,2.,2.,2.,1.,2.,1.,1.,
	2.,1.,1.,6.,1.,1.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=188;
  }

  // fiber
  else if ( modelL == "fiber" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]={36,111,190,214,235,270,338,346,372,386,
	421,424,441,470,473,483,484,498,580,594,597,660,689,735,
	742,761,762,776,779,817,860,1044,1067,1122,1238};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=1298;
  }

  // vpm1
  else if ( modelL == "vpm1" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]=
      { 180,181,182,185,195,211,214,226,231,232,244,251,263,269,
	285,294,306,307,314,  319};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=378;
  }

  // vpm2
  else if ( modelL == "vpm2" ) {
    probType=continuousWith0_1;
    int intIndicesAt1[]=
      {170,173,180,181,182,185,193,194,196,213,219,220,226,
       245,251,262,263,267,269,273,288,289,294,319,320};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=378;
  }

  // l152lav
  else if ( modelL == "l152lav" ) {
    probType=pure0_1;
    int intIndicesAt1[]={1,16,30,33,67,111,165,192,198,321,411,449,
	906,961,981,1052,1075,1107,1176,1231,1309,1415,1727,1847,
	1902,1917,1948,1950};
    int numIndices = sizeof(intIndicesAt1)/sizeof(int);
    intSoln.setConstant(numIndices,intIndicesAt1,1.0);
    expectedNumberColumns=1989;
  }

  // bell5
  else if ( modelL == "bell5" ) {
    probType=generalMip;
    int intIndicesV[]={0,1,2,3,4,6,12,13,15,33,34,36,47,48,49,50,
		       51,52,53,54,56};
    double intSolnV[] = {1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,11.,
	2.,38.,2.,498.,125.,10.,17.,41.,19.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=104;
  }

  // blend2
  else if ( modelL == "blend2" ) {
    probType=generalMip;
    int intIndicesV[]={24,35,44,45,46,52,63,64,70,71,76,84,85,
	132,134,151,152,159,164,172,173,289,300,309,310,311,
		       317,328,329,335,336,341,349,350};
    double intSolnV[] = {2.,1.,1.,1.,1.,1.,1.,1.,2.,1.,1.,1.,
	2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
	1.,1.,1.,1.,1.,1.,1.,1.,1.};
    int vecLen = sizeof(intIndicesV)/sizeof(int);
    intSoln.setVector(vecLen,intIndicesV,intSolnV);
    expectedNumberColumns=353;
  }


  // check to see if the model parameter is 
  // a known problem.
  if ( probType != undefined && si.getNumCols() == expectedNumberColumns) {

    // Specified model is a known problem
    
    numberColumns_ = si.getNumCols();

    integerVariable_= new bool[numberColumns_];
    optimalSolution_=new double[numberColumns_];
    //CoinFillN(integerVariable_, numberColumns_,0);
    //CoinFillN(optimalSolution_,numberColumns_,0.0);

    if ( probType == pure0_1 ) {
      
      // mark all variables as integer
      CoinFillN(integerVariable_,numberColumns_,true);
      // set solution to 0.0 for all not mentioned
      CoinFillN(optimalSolution_,numberColumns_,0.0);

      // mark column solution that have value 1
      for ( i=0; i<intSoln.getNumElements(); i++ ) {
        int col = intSoln.getIndices()[i];
        optimalSolution_[col] = intSoln.getElements()[i];
        assert( optimalSolution_[col]==1. );
      }

    }
    else {
      // Probtype is continuousWith0_1 or generalMip 
      assert( probType==continuousWith0_1 || probType==generalMip );

      OsiSolverInterface * siCopy = si.clone();
      assert(siCopy->getMatrixByCol()->isEquivalent(*si.getMatrixByCol()));

      // Loop once for each column looking for integer variables
      for (i=0;i<numberColumns_;i++) {

        // Is the this an integer variable?
        if(siCopy->isInteger(i)) {

          // integer variable found
          integerVariable_[i]=true;

    
          // Determine optimal solution value for integer i
          // from values saved in intSoln and probType.
          double soln;
          if ( probType==continuousWith0_1 ) {

            // Since 0_1 problem, had better be binary
            assert( siCopy->isBinary(i) );

            // intSoln only contains integers with optimal value of 1
            if ( intSoln.isExistingIndex(i) ) {
              soln = 1.0;
              assert( intSoln[i]==1. );
            }
            else {
              soln = 0.0;
            }

          } else {
            // intSoln only contains integers with nonzero optimal values
            if ( intSoln.isExistingIndex(i) ) {
              soln = intSoln[i];
              assert( intSoln[i]>=1. );
            }
            else {
              soln = 0.0;
            }
          }
          
          // Set bounds in copyied problem to fix variable to its solution     
          siCopy->setColUpper(i,soln);
          siCopy->setColLower(i,soln);
          
        }
        else {
          // this is not an integer variable
          integerVariable_[i]=false;
        }
      }
 
      // All integers have been fixed at optimal value.
      // Now solve to get continuous values
#if 0
      assert( siCopy->getNumRows()==5);        
      assert( siCopy->getNumCols()==8); 
      int r,c;
      for ( r=0; r<siCopy->getNumRows(); r++ ) {
        std::cerr <<"rhs[" <<r <<"]=" <<(si.rhs())[r] <<" " <<(siCopy->rhs())[r] <<std::endl;
      }
      for ( c=0; c<siCopy->getNumCols(); c++ ) {
        std::cerr <<"collower[" <<c <<"]=" <<(si.collower())[c] <<" " <<(siCopy->collower())[c] <<std::endl;
        std::cerr <<"colupper[" <<c <<"]=" <<(si.colupper())[c] <<" " <<(siCopy->colupper())[c] <<std::endl;
      }
#endif
      siCopy->initialSolve();
#if 0
      for ( c=0; c<siCopy->getNumCols(); c++ ) {
        std::cerr <<"colsol[" <<c <<"]=" <<optimalSolution_[c] <<" " <<(siCopy->colsol())[c] <<std::endl;
      }
      OsiRelFltEq eq;
      assert( eq(siCopy->getObjValue(),3.2368421052632));
#endif
      
      // Save column solution
      CoinCopyN(siCopy->getColSolution(),numberColumns_,optimalSolution_);

      delete siCopy;
    }
  }
 
  //if (integerVariable_!=NULL) si.rowCutDebugger_=this;

  return (integerVariable_!=NULL);
}


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiRowCutDebugger::OsiRowCutDebugger ()
: numberColumns_(0),
  integerVariable_(NULL),
  optimalSolution_(NULL)
{
  // nothing to do here
}

//-------------------------------------------------------------------
// Alternate Constructor with model name
//-------------------------------------------------------------------
// Constructor with name of model
OsiRowCutDebugger::OsiRowCutDebugger ( 
        const OsiSolverInterface & si, 
        const char * model)
: numberColumns_(0),
  integerVariable_(NULL),
  optimalSolution_(NULL)
{
  activate(si,model);
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiRowCutDebugger::OsiRowCutDebugger (const OsiRowCutDebugger & source)
{  
  // copy 
  numberColumns_=source.numberColumns_;
  integerVariable_=new bool[numberColumns_];
  optimalSolution_=new double[numberColumns_];
  CoinCopyN(source.integerVariable_,  numberColumns_, integerVariable_ );
  CoinCopyN(source.optimalSolution_, numberColumns_, optimalSolution_);
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiRowCutDebugger::~OsiRowCutDebugger ()
{
  // free memory
  delete [] integerVariable_;
  delete [] optimalSolution_;
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiRowCutDebugger &
OsiRowCutDebugger::operator=(const OsiRowCutDebugger& rhs)
{
  if (this != &rhs) {
    delete [] integerVariable_;
    delete [] optimalSolution_;
    // copy 
    numberColumns_=rhs.numberColumns_;
    integerVariable_=new bool[numberColumns_];
    optimalSolution_=new double[numberColumns_];
    CoinCopyN(rhs.integerVariable_,  numberColumns_, integerVariable_ );
    CoinCopyN(rhs.optimalSolution_, numberColumns_, optimalSolution_);
  }
  return *this;
}
