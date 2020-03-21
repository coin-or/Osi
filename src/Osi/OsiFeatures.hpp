// Copyright (C) 2020, COIN-OR Foundation
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file OsiFeatures.hpp
    \brief Defines problem features which can be used to perform algorithm
    parameter recommendation.

    This source contains a list of problem features that can be extracted from a Mixed-Integer Linear 
    Program (MIP) from an OsiSolverInterface object. All features are numeric and stored as double.
    All features are extracted in at most O(nz) time, where nz is the number of non-zeros in the 
    constraint matrix.

    Some details on algorithm recommendation for CLP can be found in:

    Vilas Boas, Matheus G.; Santos, Haroldo G.; Merschmann, Luiz H.C. and Vanden
      Berghe, Greet. Optimal Decision Trees for the Algorithm Selection Problem:
      Integer Programming Based Approaches. International Transactions in Operational
      Research, DOI 10.1111/itor.12724. 2019.

*/

class OsiSolverInterface;

/*! List of all features that are extracted */
enum OsiFeature
{   
  OFcols = 0,  //< number of columns (variables)
  OFrows = 1,  //< number of rows (constraints)
  OFcolsPerRow = 2, //< cols/rows
  OFequalities = 3, //< number of constraints =
  OFpercEqualities = 4, //< percentage of inequalities
  OFinequalities = 5, //< number of constraints <= or >=
  OFnz = 6, //< number of non-zero elements in the constraint matrix
  OFdensity = 7, //< density of the constraint matrix, i.e. (nz / (cols*rows))*100

  OFbin = 8, //< number of binary variables
  OFgenInt = 9, //< total number of integer (excluding binaries)
  OFinteger = 10, //< total number of integer (includes binaries) variables
  OFcontinuous = 11, //< number of continuous variables
  OFpercInteger = 12, //< percentage of integer variables
  OFpercBin = 13, //< percentage of binary variables

  /* row types */
  OFrowsPartitioning = 14, //< number of partitioning constraints, e.g. x1 + x2 ... = 1 (binary vars)
  OFpercRowsPartitioning = 15, //< percentage of partitioning constraints, e.g. x1 + x2 ... = 1 (binary vars)
  OFrowsPacking = 16, //< number of packing constraints, e.g. x1 + x2 ... <= 1 (binary vars)
  OFpercRowsPacking = 17, //< percentage of packing constraints, e.g. x1 + x2 ... <= 1 (binary vars)
  OFrowsCovering = 18, //< number of covering constraints, e.g. x1 + x2 ... >= 1 (binary vars)
  OFpercRowsCovering = 19, //< percentage of covering constraints, e.g. x1 + x2 ... >= 1 (binary vars)
  OFrowsCardinality = 20, //< number of cardinality constraints, e.g. x1 + x2 ... = k, k >= 2 (binary vars)
  OFpercRowsCardinality = 21, //< percentage of cardinality constraints, e.g. x1 + x2 ... = k, k >= 2 (binary vars)
  OFrowsKnapsack = 22, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, b integer >= 2 (binary vars)
  OFpercRowsKnapsack = 23, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, b integer >= 2 (binary vars)
  OFrowsIntegerKnapsack = 24, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2 (binary vars)
  OFpercRowsIntegerKnapsack = 25, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2 (binary vars)
  OFrowsInvKnapsack = 26, //< number of invariant knapsack constraints, e.g. x1 + x2 ... <= b, b >= 2 (binary vars)
  OFpercRowsInvKnapsack = 27, //< percentage of invariant knapsack constraints, e.g. x1 + x2 ... <= b, b >= 2 (binary vars)
  OFrowsSingleton = 28, //< number of constraints with only one variable
  OFpercRowsSingleton = 29, //< percentage of constraints with only one variable
  OFrowsAggr = 30, //< number of constraints in the form ax + by = c
  OFpercRowsAggr = 31, //< percentage of constraints in the form ax + by = c
  OFrowsPrec = 32, //< number of constraints in the form ax - ay <= b
  OFpercRowsPrec = 33, //< percentage of constraints in the form ax - ay <= b
  OFrowsVarBnd = 34, //< number of constraints with only one variable
  OFpercRowsVarBnd = 35, //< percentage of constraints with only one variable
  OFrowsBinPacking = 36, //< number of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2, at least one cj >= 2 (binary vars)
  OFpercRowsBinPacking = 37, //< percentage of knapsack constraints, e.g. c1*x1 + c2*x2 ... <= b, cj and b integer, b >= 2, at least one cj >= 2 (binary vars)
  OFrowsMixedBin = 38, //< constraint that involves binary and continuous variables
  OFpercRowsMixedBin = 39, //< percentage constraint that involves binary and continuous variables
  OFrowsGenInt = 40,  //< constraints with some general integer (not binary) variables
  OFpercRowsGenInt = 41,  //< percentage constraints with some general integer (not binary) variables
  OFrowsFlowBin = 42, //< equality, at least 2 positive and 2 negative coefficients, only bin vars
  OFpercRowsFlowBin = 43, //< percentage equality, at least 2 positive and 2 negative coefficients, only bin vars
  OFrowsFlowMx = 44, //< equality, at least 2 positive and 2 negative coefficients
  OFpercRowsFlowMx = 45, //< equality, at least 2 positive and 2 negative coefficients

  /* statistics constraint matrix */
  OFaMin = 46,
  OFaMax = 47,
  OFaAvg = 48,
  OFaStdDev = 49,
  OFaRatioLSA = 50,
  OFaAllInt = 51,
  OFaPercInt = 52,

  /* statistics objective function */
  OFobjMin = 53,
  OFobjMax = 54,
  OFobjAvg = 55,
  OFobjStdDev = 56,
  OFobjRatioLSA = 57,
  OFobjAllInt = 58,
  OFobjPercInt = 59,

  /* statistics right hand side */
  OFrhsMin = 60,
  OFrhsMax = 61,
  OFrhsAvg = 62,
  OFrhsStdDev = 63,
  OFrhsRatioLSA = 64,
  OFrhsAllInt = 65,
  OFrhsPercInt = 66,

  /* statistics non-zero distribution rows */
  OFrowNzMin = 67,
  OFrowNzMax = 68,
  OFrowNzAvg = 69,
  OFrowNzStdDev = 70,

  /* statistics non-zero distribution cols */
  OFcolNzMin = 71,
  OFcolNzMax = 72,
  OFcolNzAvg = 73,
  OFcolNzStdDev = 74,

  // constraints with nz less or equal
  OFrowsLess4Nz = 75,
  OFrowsLess8Nz = 76,
  OFrowsLess16Nz = 77,
  OFrowsLess32Nz = 78,
  OFrowsLess64Nz = 79,
  OFrowsLess128Nz = 80,
  OFrowsLess256Nz = 81,
  OFrowsLess512Nz = 82,
  OFrowsLess1024Nz = 83,
  OFpercRowsLess4Nz = 84,
  OFpercRowsLess8Nz = 85,
  OFpercRowsLess16Nz = 86,
  OFpercRowsLess32Nz = 87,
  OFpercRowsLess64Nz = 88,
  OFpercRowsLess128Nz = 89,
  OFpercRowsLess256Nz = 90,
  OFpercRowsLess512Nz = 91,
  OFpercRowsLess1024Nz = 92,

  // constraints nz at least
  OFrowsLeast4Nz = 93,
  OFrowsLeast8Nz = 94,
  OFrowsLeast16Nz = 95,
  OFrowsLeast32Nz = 96,
  OFrowsLeast64Nz = 97,
  OFrowsLeast128Nz = 98,
  OFrowsLeast256Nz = 99,
  OFrowsLeast512Nz = 100,
  OFrowsLeast1024Nz = 101,
  OFrowsLeast2048Nz = 102,
  OFrowsLeast4096Nz = 103,
  OFpercRowsLeast4Nz = 104,
  OFpercRowsLeast8Nz = 105,
  OFpercRowsLeast16Nz = 106,
  OFpercRowsLeast32Nz = 107,
  OFpercRowsLeast64Nz = 108,
  OFpercRowsLeast128Nz = 109,
  OFpercRowsLeast256Nz = 110,
  OFpercRowsLeast512Nz = 111,
  OFpercRowsLeast1024Nz = 112,
  OFpercRowsLeast2048Nz = 113,
  OFpercRowsLeast4096Nz = 114,

  // constraints with nz less or equal
  OFcolsLess4Nz = 115,
  OFcolsLess8Nz = 116,
  OFcolsLess16Nz = 117,
  OFcolsLess32Nz = 118,
  OFcolsLess64Nz = 119,
  OFcolsLess128Nz = 120,
  OFcolsLess256Nz = 121,
  OFcolsLess512Nz = 122,
  OFcolsLess1024Nz = 123,
  OFpercColsLess4Nz = 124,
  OFpercColsLess8Nz = 125,
  OFpercColsLess16Nz = 126,
  OFpercColsLess32Nz = 127,
  OFpercColsLess64Nz = 128,
  OFpercColsLess128Nz = 129,
  OFpercColsLess256Nz = 130,
  OFpercColsLess512Nz = 131,
  OFpercColsLess1024Nz = 132,

  // constraints nz at least
  OFcolsLeast4Nz = 133,
  OFcolsLeast8Nz = 134,
  OFcolsLeast16Nz = 135,
  OFcolsLeast32Nz = 136,
  OFcolsLeast64Nz = 137,
  OFcolsLeast128Nz = 138,
  OFcolsLeast256Nz = 139,
  OFcolsLeast512Nz = 140,
  OFcolsLeast1024Nz = 141,
  OFcolsLeast2048Nz = 142,
  OFcolsLeast4096Nz = 143,

  OFpercColsLeast4Nz = 144,
  OFpercColsLeast8Nz = 145,
  OFpercColsLeast16Nz = 146,
  OFpercColsLeast32Nz = 147,
  OFpercColsLeast64Nz = 148,
  OFpercColsLeast128Nz = 149,
  OFpercColsLeast256Nz = 150,
  OFpercColsLeast512Nz = 151,
  OFpercColsLeast1024Nz = 152,
  OFpercColsLeast2048Nz = 153,
  OFpercColsLeast4096Nz = 154,

  OFCount = 155 //< Number of features
};

class OsiFeatures {
public:
  /** @brief number of features */
  static int n; 

  /** @brief name of the i-th feature */
  static const char *name(int i);

  /** @brief name of an specific feature feature */
  static const char *name( const OsiFeature of );

  /** @brief computes all feature values, the size of this vector should be at least OFCount */
  static void compute(double *features, OsiSolverInterface *solver);
};


/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
