// Copyright (C) 2026, COIN-OR Foundation
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Tests for OsiFeatures constraint-type identification.
//
// Constraint types tested (all require binary-only variables):
//
//   Unit-coefficient (all ±1 coefs):
//     - Set packing      : +1 coefs, <= 1
//     - Set partitioning : +1 coefs, =  1
//     - Set covering     : +1 coefs, >= 1
//     - Cardinality      : +1 coefs, =  k  (k >= 2)
//     - Invariant knapsack: +1 coefs, <= k  (k >= 2)
//
//   Mixed-coefficient:
//     - Knapsack         : different positive weights, <= b (b >= ~1.1)
//     - Integer knapsack : same + all coefs integer
//     - Bin packing      : exactly 1 negative coef, nz >= 2, rhs >= ~1.1
//
// Each group is tested in canonical form, sign-flipped form where the
// code supports it, and negative cases.

#include "OsiClpSolverInterface.hpp"
#include "OsiFeatures.hpp"
#include "CoinPackedMatrix.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

static bool eq(double a, double b)
{
  return fabs(a - b) < 1e-8;
}

static int nErrors = 0;
static int nTests = 0;

static void check(const char *test, const char *feat, double got, double expected)
{
  ++nTests;
  if (!eq(got, expected)) {
    fprintf(stderr, "FAIL %s: %s = %g, expected %g\n", test, feat, got, expected);
    ++nErrors;
  }
}

// Helper: load a column-major model with all-binary variables.
static void loadBinary(OsiClpSolverInterface &si,
  int ncols, int nrows,
  const CoinBigIndex *start, const int *index, const double *value,
  const char *sense, const double *rhs)
{
  std::vector<double> collb(ncols, 0.0), colub(ncols, 1.0), obj(ncols, 0.0);
  std::vector<double> range(nrows, 0.0);
  si.loadProblem(ncols, nrows, start, index, value,
    collb.data(), colub.data(), obj.data(), sense, rhs, range.data());
  for (int j = 0; j < ncols; ++j)
    si.setInteger(j);
}

// -----------------------------------------------------------------------
// Set packing / partitioning / covering — canonical form
// -----------------------------------------------------------------------
//   row 0: x0 + x1 + x2       <= 1  (packing)
//   row 1: x0 + x1 + x2 + x3  = 1  (partitioning)
//   row 2: x1 + x2 + x3       >= 1  (covering)
static int testSPCCanonical()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 5, 8, 10 };
  int idx[] = { 0, 1, 0, 1, 2, 0, 1, 2, 1, 2 };
  double val[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  char sense[] = { 'L', 'E', 'G' };
  double rhs[] = { 1, 1, 1 };
  loadBinary(si, 4, 3, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("SPCCanonical", "packing", f[OFrowsPacking], 1);
  check("SPCCanonical", "partitioning", f[OFrowsPartitioning], 1);
  check("SPCCanonical", "covering", f[OFrowsCovering], 1);
  if (nErrors == e0)
    printf("PASS testSPCCanonical\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Set packing / partitioning / covering — negated (sign-flipped)
// -----------------------------------------------------------------------
//   row 0: -x0 - x1 - x2       >= -1  (packing)
//   row 1: -x0 - x1 - x2 - x3  = -1  (partitioning)
//   row 2: -x1 - x2 - x3       <= -1  (covering)
static int testSPCNegated()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 5, 8, 10 };
  int idx[] = { 0, 1, 0, 1, 2, 0, 1, 2, 1, 2 };
  double val[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
  char sense[] = { 'G', 'E', 'L' };
  double rhs[] = { -1, -1, -1 };
  loadBinary(si, 4, 3, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("SPCNegated", "packing", f[OFrowsPacking], 1);
  check("SPCNegated", "partitioning", f[OFrowsPartitioning], 1);
  check("SPCNegated", "covering", f[OFrowsCovering], 1);
  if (nErrors == e0)
    printf("PASS testSPCNegated\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Negative case: no SPC rows
// -----------------------------------------------------------------------
//   row 0: 2*x0 + 3*x1 <= 5  (knapsack)
//   row 1: x0 + x1      <= 3  (inv-knapsack, rhs != 1)
static int testNoSPC()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 4 };
  int idx[] = { 0, 1, 0, 1 };
  double val[] = { 2, 1, 3, 1 };
  char sense[] = { 'L', 'L' };
  double rhs[] = { 5, 3 };
  loadBinary(si, 2, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("NoSPC", "packing", f[OFrowsPacking], 0);
  check("NoSPC", "partitioning", f[OFrowsPartitioning], 0);
  check("NoSPC", "covering", f[OFrowsCovering], 0);
  if (nErrors == e0)
    printf("PASS testNoSPC\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Cardinality and invariant knapsack — canonical
// -----------------------------------------------------------------------
//   row 0: x0 + x1 + x2 + x3  = 3   (cardinality)
//   row 1: x0 + x1 + x2 + x3 <= 2   (invariant knapsack)
static int testCardinalityInvKnapsack()
{
  OsiClpSolverInterface si;
  // 4 cols, each appears in both rows
  CoinBigIndex start[] = { 0, 2, 4, 6, 8 };
  int idx[] = { 0, 1, 0, 1, 0, 1, 0, 1 };
  double val[] = { 1, 1, 1, 1, 1, 1, 1, 1 };
  char sense[] = { 'E', 'L' };
  double rhs[] = { 3, 2 };
  loadBinary(si, 4, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("CardInvKnap", "cardinality", f[OFrowsCardinality], 1);
  check("CardInvKnap", "invKnapsack", f[OFrowsInvKnapsack], 1);
  // These should NOT be classified as SPC
  check("CardInvKnap", "packing", f[OFrowsPacking], 0);
  check("CardInvKnap", "partitioning", f[OFrowsPartitioning], 0);
  check("CardInvKnap", "covering", f[OFrowsCovering], 0);
  if (nErrors == e0)
    printf("PASS testCardinalityInvKnapsack\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Cardinality and invariant knapsack — negated
// -----------------------------------------------------------------------
//   row 0: -x0 - x1 - x2 - x3  = -3  (cardinality, negated)
//   row 1: -x0 - x1 - x2 - x3 >= -2  (invariant knapsack, negated)
static int testCardinalityInvKnapsackNegated()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 4, 6, 8 };
  int idx[] = { 0, 1, 0, 1, 0, 1, 0, 1 };
  double val[] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  char sense[] = { 'E', 'G' };
  double rhs[] = { -3, -2 };
  loadBinary(si, 4, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("CardInvKnapNeg", "cardinality", f[OFrowsCardinality], 1);
  check("CardInvKnapNeg", "invKnapsack", f[OFrowsInvKnapsack], 1);
  if (nErrors == e0)
    printf("PASS testCardinalityInvKnapsackNegated\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Knapsack and integer knapsack
// -----------------------------------------------------------------------
// Conditions: all binary vars, different positive weights (max-min >= 0.1),
//             no negative coefs, rhs >= 1.1
//
//   row 0: 2*x0 + 3*x1 + 5*x2 <= 7   (knapsack + integer knapsack)
//   row 1: 1.5*x0 + 2.7*x1    <= 4   (knapsack, NOT integer knapsack)
static int testKnapsack()
{
  OsiClpSolverInterface si;
  // 3 cols for row 0, 2 cols for row 1 (cols 0,1 in both, col 2 in row 0 only)
  CoinBigIndex start[] = { 0, 2, 4, 5 };
  int idx[] = { 0, 1, 0, 1, 0 };
  double val[] = { 2, 1.5, 3, 2.7, 5 };
  char sense[] = { 'L', 'L' };
  double rhs[] = { 7, 4 };
  loadBinary(si, 3, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("Knapsack", "knapsack", f[OFrowsKnapsack], 2);
  check("Knapsack", "integerKnapsack", f[OFrowsIntegerKnapsack], 1);
  // Should not be SPC or inv-knapsack
  check("Knapsack", "packing", f[OFrowsPacking], 0);
  check("Knapsack", "invKnapsack", f[OFrowsInvKnapsack], 0);
  if (nErrors == e0)
    printf("PASS testKnapsack\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Knapsack negative case: equal weights → not knapsack (max-min < 0.1)
// -----------------------------------------------------------------------
//   row 0: 3*x0 + 3*x1 + 3*x2 <= 6   (all same weight → inv-knapsack path)
static int testKnapsackNegative()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 1, 2, 3 };
  int idx[] = { 0, 0, 0 };
  double val[] = { 3, 3, 3 };
  char sense[] = { 'L' };
  double rhs[] = { 6 };
  loadBinary(si, 3, 1, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("KnapsackNeg", "knapsack", f[OFrowsKnapsack], 0);
  check("KnapsackNeg", "integerKnapsack", f[OFrowsIntegerKnapsack], 0);
  if (nErrors == e0)
    printf("PASS testKnapsackNegative\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Bin packing: exactly 1 negative coef, nz >= 2, rhs >= 1.1, all binary
// -----------------------------------------------------------------------
//   row 0: 3*x0 + 2*x1 - 4*x2 <= 5   (bin packing)
static int testBinPacking()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 1, 2, 3 };
  int idx[] = { 0, 0, 0 };
  double val[] = { 3, 2, -4 };
  char sense[] = { 'L' };
  double rhs[] = { 5 };
  loadBinary(si, 3, 1, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("BinPacking", "binPacking", f[OFrowsBinPacking], 1);
  if (nErrors == e0)
    printf("PASS testBinPacking\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Bin packing negative: 2 negative coefs → not bin packing
// -----------------------------------------------------------------------
//   row 0: 3*x0 - 2*x1 - 4*x2 <= 5
static int testBinPackingNegative()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 1, 2, 3 };
  int idx[] = { 0, 0, 0 };
  double val[] = { 3, -2, -4 };
  char sense[] = { 'L' };
  double rhs[] = { 5 };
  loadBinary(si, 3, 1, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("BinPackingNeg", "binPacking", f[OFrowsBinPacking], 0);
  if (nErrors == e0)
    printf("PASS testBinPackingNegative\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Combined model: one of each type, verify no cross-contamination
// -----------------------------------------------------------------------
//   row 0: x0 + x1 + x2       <= 1   (packing)
//   row 1: x0 + x1 + x2        = 1   (partitioning)
//   row 2: x0 + x1 + x2       >= 1   (covering)
//   row 3: x0 + x1 + x2 + x3  = 3   (cardinality)
//   row 4: x0 + x1 + x2 + x3 <= 2   (inv-knapsack)
//   row 5: 2*x0 + 5*x1        <= 6   (knapsack + int-knapsack)
//   row 6: 3*x0 + 2*x1 - x2   <= 4   (bin packing)
static int testCombined()
{
  OsiClpSolverInterface si;
  const int ncols = 4;
  const int nrows = 7;

  CoinPackedMatrix mtx(false, 0, 0); // row-major
  mtx.setDimensions(0, ncols);

  {
    // row 0: x0+x1+x2 <= 1
    int c[] = { 0, 1, 2 };
    double v[] = { 1, 1, 1 };
    mtx.appendRow(3, c, v);
  }
  {
    // row 1: x0+x1+x2 = 1
    int c[] = { 0, 1, 2 };
    double v[] = { 1, 1, 1 };
    mtx.appendRow(3, c, v);
  }
  {
    // row 2: x0+x1+x2 >= 1
    int c[] = { 0, 1, 2 };
    double v[] = { 1, 1, 1 };
    mtx.appendRow(3, c, v);
  }
  {
    // row 3: x0+x1+x2+x3 = 3
    int c[] = { 0, 1, 2, 3 };
    double v[] = { 1, 1, 1, 1 };
    mtx.appendRow(4, c, v);
  }
  {
    // row 4: x0+x1+x2+x3 <= 2
    int c[] = { 0, 1, 2, 3 };
    double v[] = { 1, 1, 1, 1 };
    mtx.appendRow(4, c, v);
  }
  {
    // row 5: 2*x0 + 5*x1 <= 6
    int c[] = { 0, 1 };
    double v[] = { 2, 5 };
    mtx.appendRow(2, c, v);
  }
  {
    // row 6: 3*x0 + 2*x1 - x2 <= 4
    int c[] = { 0, 1, 2 };
    double v[] = { 3, 2, -1 };
    mtx.appendRow(3, c, v);
  }

  double collb[] = { 0, 0, 0, 0 };
  double colub[] = { 1, 1, 1, 1 };
  double obj[] = { 0, 0, 0, 0 };
  char sense[] = { 'L', 'E', 'G', 'E', 'L', 'L', 'L' };
  double rhs[] = { 1, 1, 1, 3, 2, 6, 4 };
  double range[] = { 0, 0, 0, 0, 0, 0, 0 };

  si.loadProblem(mtx, collb, colub, obj, sense, rhs, range);
  for (int j = 0; j < ncols; ++j)
    si.setInteger(j);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("Combined", "packing", f[OFrowsPacking], 1);
  check("Combined", "partitioning", f[OFrowsPartitioning], 1);
  check("Combined", "covering", f[OFrowsCovering], 1);
  check("Combined", "cardinality", f[OFrowsCardinality], 1);
  check("Combined", "invKnapsack", f[OFrowsInvKnapsack], 1);
  check("Combined", "knapsack", f[OFrowsKnapsack], 1);
  check("Combined", "integerKnapsack", f[OFrowsIntegerKnapsack], 1);
  check("Combined", "binPacking", f[OFrowsBinPacking], 1);
  if (nErrors == e0)
    printf("PASS testCombined\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Singleton: nz == 1
// -----------------------------------------------------------------------
//   row 0: 5*x0 <= 3  (singleton)
//   row 1: x0 + x1 <= 1  (NOT singleton)
static int testSingleton()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 3 };
  int idx[] = { 0, 1, 1 };
  double val[] = { 5, 1, 1 };
  char sense[] = { 'L', 'L' };
  double rhs[] = { 3, 1 };
  loadBinary(si, 2, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("Singleton", "singleton", f[OFrowsSingleton], 1);
  if (nErrors == e0)
    printf("PASS testSingleton\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Aggregation: nz == 2, equality
// -----------------------------------------------------------------------
//   row 0: 3*x0 + 2*x1 = 5  (aggregation)
//   row 1: 3*x0 + 2*x1 <= 5  (NOT aggregation — inequality)
static int testAggregation()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 4 };
  int idx[] = { 0, 1, 0, 1 };
  double val[] = { 3, 3, 2, 2 };
  char sense[] = { 'E', 'L' };
  double rhs[] = { 5, 5 };
  loadBinary(si, 2, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("Aggregation", "aggr", f[OFrowsAggr], 1);
  if (nErrors == e0)
    printf("PASS testAggregation\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Variable bound: nz == 2, exactly 1 binary variable
// -----------------------------------------------------------------------
//   row 0: x_bin + 3*y_cont <= 5  (variable bound)
//   row 1: x_bin + y_bin <= 1     (NOT var bound — both binary)
static int testVarBound()
{
  OsiClpSolverInterface si;
  // 3 cols: col0 = binary, col1 = continuous, col2 = binary
  const int ncols = 3, nrows = 2;
  CoinBigIndex start[] = { 0, 2, 3, 4 };
  int idx[] = { 0, 1, 0, 1 };
  double val[] = { 1, 1, 3, 1 };
  double collb[] = { 0, 0, 0 };
  double colub[] = { 1, 10, 1 };
  double obj[] = { 0, 0, 0 };
  char sense[] = { 'L', 'L' };
  double rhs[] = { 5, 1 };
  double range[] = { 0, 0 };
  si.loadProblem(ncols, nrows, start, idx, val,
    collb, colub, obj, sense, rhs, range);
  si.setInteger(0);
  // col1 stays continuous
  si.setInteger(2);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("VarBound", "varBnd", f[OFrowsVarBnd], 1);
  if (nErrors == e0)
    printf("PASS testVarBound\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Precedence: nz == 2, same-type vars, 1 pos + 1 neg, equal abs coefs
//   ax - ay <= b
// -----------------------------------------------------------------------
//   row 0: 3*x0 - 3*x1 <= 5  (should be precedence)
//
// NOTE: the current code checks dbl_equal(summRow.minV, summRow.maxV)
// which compares -3 == 3 → always false.  This is a known bug.
// The correct check should be dbl_equal(-summRow.minV, summRow.maxV).
// This test documents the current (broken) behavior: 0 precedence rows.
// When the bug is fixed, update expected value to 1.
static int testPrecedence()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 1, 2 };
  int idx[] = { 0, 0 };
  double val[] = { 3, -3 };
  char sense[] = { 'L' };
  double rhs[] = { 5 };
  loadBinary(si, 2, 1, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  // BUG: should be 1, but current code always yields 0
  check("Precedence", "prec", f[OFrowsPrec], 0);
  if (nErrors == e0)
    printf("PASS testPrecedence (documents known bug: prec always 0)\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Flow (binary): all binary, >= 2 pos + >= 2 neg coefs, equality
// -----------------------------------------------------------------------
//   row 0: x0 + x1 - x2 - x3 = 0  (flow binary)
//   row 1: x0 + x1 - x2 - x3 <= 0  (NOT flow — inequality)
static int testFlowBin()
{
  OsiClpSolverInterface si;
  CoinBigIndex start[] = { 0, 2, 4, 6, 8 };
  int idx[] = { 0, 1, 0, 1, 0, 1, 0, 1 };
  double val[] = { 1, 1, 1, 1, -1, -1, -1, -1 };
  char sense[] = { 'E', 'L' };
  double rhs[] = { 0, 0 };
  loadBinary(si, 4, 2, start, idx, val, sense, rhs);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("FlowBin", "flowBin", f[OFrowsFlowBin], 1);
  if (nErrors == e0)
    printf("PASS testFlowBin\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Flow (mixed): NOT all binary, >= 2 pos + >= 2 neg coefs, equality
// -----------------------------------------------------------------------
//   row 0: x0_cont + x1_cont - x2_bin - x3_bin = 0  (flow mixed)
static int testFlowMixed()
{
  OsiClpSolverInterface si;
  const int ncols = 4, nrows = 1;
  CoinBigIndex start[] = { 0, 1, 2, 3, 4 };
  int idx[] = { 0, 0, 0, 0 };
  double val[] = { 1, 1, -1, -1 };
  double collb[] = { 0, 0, 0, 0 };
  double colub[] = { 100, 100, 1, 1 };
  double obj[] = { 0, 0, 0, 0 };
  char sense[] = { 'E' };
  double rhs[] = { 0 };
  double range[] = { 0 };
  si.loadProblem(ncols, nrows, start, idx, val,
    collb, colub, obj, sense, rhs, range);
  si.setInteger(2);
  si.setInteger(3);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("FlowMixed", "flowMx", f[OFrowsFlowMx], 1);
  check("FlowMixed", "flowBin", f[OFrowsFlowBin], 0);
  if (nErrors == e0)
    printf("PASS testFlowMixed\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Mixed binary: constraint with both binary and continuous variables
// -----------------------------------------------------------------------
//   row 0: x0_bin + 3*y_cont <= 5
static int testMixedBin()
{
  OsiClpSolverInterface si;
  const int ncols = 2, nrows = 1;
  CoinBigIndex start[] = { 0, 1, 2 };
  int idx[] = { 0, 0 };
  double val[] = { 1, 3 };
  double collb[] = { 0, 0 };
  double colub[] = { 1, 100 };
  double obj[] = { 0, 0 };
  char sense[] = { 'L' };
  double rhs[] = { 5 };
  double range[] = { 0 };
  si.loadProblem(ncols, nrows, start, idx, val,
    collb, colub, obj, sense, rhs, range);
  si.setInteger(0);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("MixedBin", "mixedBin", f[OFrowsMixedBin], 1);
  if (nErrors == e0)
    printf("PASS testMixedBin\n");
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// General integer: constraint with general integer (non-binary) variables
// -----------------------------------------------------------------------
//   row 0: 2*x0_genint + 3*x1_genint <= 10
static int testGenInt()
{
  OsiClpSolverInterface si;
  const int ncols = 2, nrows = 1;
  CoinBigIndex start[] = { 0, 1, 2 };
  int idx[] = { 0, 0 };
  double val[] = { 2, 3 };
  double collb[] = { 0, 0 };
  double colub[] = { 10, 10 };
  double obj[] = { 0, 0 };
  char sense[] = { 'L' };
  double rhs[] = { 10 };
  double range[] = { 0 };
  si.loadProblem(ncols, nrows, start, idx, val,
    collb, colub, obj, sense, rhs, range);
  si.setInteger(0);
  si.setInteger(1);

  double f[OFCount];
  OsiFeatures::compute(f, &si);

  int e0 = nErrors;
  check("GenInt", "genInt", f[OFrowsGenInt], 1);
  check("GenInt", "packing", f[OFrowsPacking], 0);
  if (nErrors == e0)
    printf("PASS testGenInt\n");
  return nErrors - e0;
}

int main()
{
  testSPCCanonical();
  testSPCNegated();
  testNoSPC();
  testCardinalityInvKnapsack();
  testCardinalityInvKnapsackNegated();
  testKnapsack();
  testKnapsackNegative();
  testBinPacking();
  testBinPackingNegative();
  testSingleton();
  testAggregation();
  testVarBound();
  testPrecedence();
  testFlowBin();
  testFlowMixed();
  testMixedBin();
  testGenInt();
  testCombined();

  printf("\n%d tests, %d failures\n", nTests, nErrors);
  if (nErrors) {
    fprintf(stderr, "%d test(s) FAILED\n", nErrors);
    return 1;
  }
  printf("All OsiFeatures constraint-type tests passed.\n");
  return 0;
}
