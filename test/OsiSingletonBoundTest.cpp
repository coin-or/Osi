// Copyright (C) 2026, COIN-OR Foundation
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Tests for OsiSolverInterface::tightenBoundsFromSingletonRows.
//
// Cases covered:
//   - Positive coefficient, tightens upper bound
//   - Positive coefficient, tightens lower bound
//   - Negative coefficient (bounds flip)
//   - One-sided row (only upper or lower bound finite)
//   - Multiple singleton rows on the same variable (tightest wins)
//   - Contradictory singleton rows (infeasibility => returns -1)
//   - Near-zero coefficient (should be skipped)
//   - Non-singleton rows are ignored
//   - nFixed count: variable not previously fixed, becomes fixed
//   - Already-fixed variable does not increment nFixed
//   - No tightening when implied bounds are not tighter

#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

static int nErrors = 0;
static int nTests = 0;

static bool eq(double a, double b, double tol = 1e-8)
{
  return fabs(a - b) <= tol;
}

static void check(const char *test, const char *what, double got, double expected)
{
  ++nTests;
  if (!eq(got, expected)) {
    fprintf(stderr, "FAIL %s: %s = %.10g, expected %.10g\n", test, what, got, expected);
    ++nErrors;
  }
}

static void checkInt(const char *test, const char *what, int got, int expected)
{
  ++nTests;
  if (got != expected) {
    fprintf(stderr, "FAIL %s: %s = %d, expected %d\n", test, what, got, expected);
    ++nErrors;
  }
}

// Build a solver with one variable (col bounds colLb..colUb) and one row.
// The row has coefficient `a` on x and bounds rowLb <= a*x <= rowUb.
static OsiClpSolverInterface makeSingleVar(
  double colLb, double colUb,
  double a,
  double rowLb, double rowUb)
{
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);
  const double inf = si.getInfinity();
  CoinBigIndex start[] = { 0, 1 };
  int idx[] = { 0 };
  double val[] = { a };
  double lo[] = { colLb };
  double up[] = { colUb };
  double obj[] = { 0.0 };
  double rlo[] = { rowLb };
  double rup[] = { rowUb };
  si.loadProblem(1, 1, start, idx, val, lo, up, obj, rlo, rup);
  return si;
}

// -----------------------------------------------------------------------
// Positive coefficient — tightens upper bound
//   x in [0, 10], row: x <= 3  => upper becomes 3
// -----------------------------------------------------------------------
static int testPositiveCoefUpper()
{
  const char *name = "testPositiveCoefUpper";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1.0, -inf, 3.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  checkInt(name, "nFixed",     nFixed, 0);
  check   (name, "colLb",      si.getColLower()[0], 0.0);
  check   (name, "colUb",      si.getColUpper()[0], 3.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Positive coefficient — tightens lower bound
//   x in [0, 10], row: x >= 5  => lower becomes 5
// -----------------------------------------------------------------------
static int testPositiveCoefLower()
{
  const char *name = "testPositiveCoefLower";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1.0, 5.0, inf);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  checkInt(name, "nFixed",     nFixed, 0);
  check   (name, "colLb",      si.getColLower()[0], 5.0);
  check   (name, "colUb",      si.getColUpper()[0], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Positive coefficient, coefficient != 1
//   x in [0, 10], row: 2x <= 6  => upper becomes 3
// -----------------------------------------------------------------------
static int testPositiveCoefScale()
{
  const char *name = "testPositiveCoefScale";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 2.0, -inf, 6.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  check   (name, "colUb",      si.getColUpper()[0], 3.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Negative coefficient — bounds flip
//   x in [0, 10], row: -x <= -2  => lower becomes 2
// -----------------------------------------------------------------------
static int testNegativeCoef()
{
  const char *name = "testNegativeCoef";
  const double inf = 1e30;
  // -1*x <= -2  means x >= 2
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, -1.0, -inf, -2.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  check   (name, "colLb",      si.getColLower()[0], 2.0);
  check   (name, "colUb",      si.getColUpper()[0], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Ranged row tightens both bounds
//   x in [0, 10], row: 2 <= x <= 7  => both bounds tightened
// -----------------------------------------------------------------------
static int testRangedRowBothBounds()
{
  const char *name = "testRangedRowBothBounds";
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1.0, 2.0, 7.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  checkInt(name, "nFixed",     nFixed, 0);
  check   (name, "colLb",      si.getColLower()[0], 2.0);
  check   (name, "colUb",      si.getColUpper()[0], 7.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Multiple singleton rows on the same variable — tightest wins
//   x in [0, 10], row0: x <= 8, row1: x <= 5  => upper becomes 5
// -----------------------------------------------------------------------
static int testMultipleSingletonsSameVar()
{
  const char *name = "testMultipleSingletonsSameVar";
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);
  const double inf = si.getInfinity();

  // 1 variable, 2 rows, both singleton
  CoinBigIndex start[] = { 0, 1 };
  int idx[] = { 0 };
  double val[] = { 1.0 };
  double lo[] = { 0.0 };
  double up[] = { 10.0 };
  double obj[] = { 0.0 };
  double rlo[] = { -inf, -inf };
  double rup[] = { 8.0, 5.0 };

  // We need to load 2 rows. Use CoinPackedMatrix approach.
  CoinPackedMatrix colMat(false, 0, 0);
  colMat.setDimensions(0, 1);
  int idx0[] = { 0 };
  double v0[] = { 1.0 };
  int idx1[] = { 0 };
  double v1[] = { 1.0 };
  colMat.appendRow(1, idx0, v0);
  colMat.appendRow(1, idx1, v1);

  si.loadProblem(colMat, lo, up, obj, rlo, rup);

  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  check   (name, "colUb",      si.getColUpper()[0], 5.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Contradictory singleton rows — infeasibility detected (returns -1)
//   x in [0, 10], row0: x >= 8, row1: x <= 3
// -----------------------------------------------------------------------
static int testContradictoryRows()
{
  const char *name = "testContradictoryRows";
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);
  const double inf = si.getInfinity();

  CoinPackedMatrix colMat(false, 0, 0);
  colMat.setDimensions(0, 1);
  int idx0[] = { 0 };
  double v0[] = { 1.0 };
  colMat.appendRow(1, idx0, v0); // x >= 8
  colMat.appendRow(1, idx0, v0); // x <= 3
  double lo[] = { 0.0 };
  double up[] = { 10.0 };
  double obj[] = { 0.0 };
  double rlo[] = { 8.0, -inf };
  double rup[] = { inf,  3.0 };
  si.loadProblem(colMat, lo, up, obj, rlo, rup);

  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "return -1 (infeasible)", nTight, -1);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Near-zero coefficient is skipped — no tightening
//   x in [0, 10], row: 1e-10 * x <= 5  => no change
// -----------------------------------------------------------------------
static int testNearZeroCoef()
{
  const char *name = "testNearZeroCoef";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1e-10, -inf, 5.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 0);
  check   (name, "colUb",      si.getColUpper()[0], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Non-singleton rows are ignored
//   x0 in [0,10], x1 in [0,10], row: x0 + x1 <= 5  => no change
// -----------------------------------------------------------------------
static int testNonSingletonIgnored()
{
  const char *name = "testNonSingletonIgnored";
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);
  const double inf = si.getInfinity();

  // Column-major: col 0 in row 0 (val 1), col 1 in row 0 (val 1)
  CoinBigIndex start[] = { 0, 1, 2 }; // numcols+1 entries
  int idx[] = { 0, 0 };
  double val[] = { 1.0, 1.0 };
  double lo[] = { 0.0, 0.0 };
  double up[] = { 10.0, 10.0 };
  double obj[] = { 0.0, 0.0 };
  double rlo[] = { -inf };
  double rup[] = { 5.0 };
  si.loadProblem(2, 1, start, idx, val, lo, up, obj, rlo, rup);

  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 0);
  check   (name, "col0Ub",     si.getColUpper()[0], 10.0);
  check   (name, "col1Ub",     si.getColUpper()[1], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// nFixed: variable not previously fixed, becomes fixed
//   x in [0, 10], row: 4 <= x <= 4  => fixed at 4
// -----------------------------------------------------------------------
static int testBecomesFixed()
{
  const char *name = "testBecomesFixed";
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1.0, 4.0, 4.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  checkInt(name, "nFixed",     nFixed, 1);
  check   (name, "colLb",      si.getColLower()[0], 4.0);
  check   (name, "colUb",      si.getColUpper()[0], 4.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Already-fixed variable: tighter-or-equal implied bounds don't add to nFixed
//   x fixed at 5 (lb=ub=5), row: x <= 7 => no change at all
// -----------------------------------------------------------------------
static int testAlreadyFixed()
{
  const char *name = "testAlreadyFixed";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(5.0, 5.0, 1.0, -inf, 7.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 0);
  checkInt(name, "nFixed",     nFixed, 0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Implied bound not tighter than current — no change
//   x in [0, 10], row: x <= 15  => no tightening
// -----------------------------------------------------------------------
static int testNoTightening()
{
  const char *name = "testNoTightening";
  const double inf = 1e30;
  OsiClpSolverInterface si = makeSingleVar(0.0, 10.0, 1.0, -inf, 15.0);
  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 0);
  checkInt(name, "nFixed",     nFixed, 0);
  check   (name, "colUb",      si.getColUpper()[0], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

// -----------------------------------------------------------------------
// Mixed: one singleton row and one non-singleton row
//   x0 in [0,10], x1 in [0,10]
//   row0: x0 <= 3  (singleton → tighten x0)
//   row1: x0 + x1 <= 8  (non-singleton → ignore)
// -----------------------------------------------------------------------
static int testMixedRows()
{
  const char *name = "testMixedRows";
  OsiClpSolverInterface si;
  si.messageHandler()->setLogLevel(0);
  const double inf = si.getInfinity();

  CoinPackedMatrix colMat(false, 0, 0);
  colMat.setDimensions(0, 2);
  int idx0[] = { 0 };
  double v0[] = { 1.0 };
  colMat.appendRow(1, idx0, v0); // singleton: x0 <= 3
  int idx1[] = { 0, 1 };
  double v1[] = { 1.0, 1.0 };
  colMat.appendRow(2, idx1, v1); // non-singleton
  double lo[] = { 0.0, 0.0 };
  double up[] = { 10.0, 10.0 };
  double obj[] = { 0.0, 0.0 };
  double rlo[] = { -inf, -inf };
  double rup[] = { 3.0, 8.0 };
  si.loadProblem(colMat, lo, up, obj, rlo, rup);

  int nFixed = -1;
  int nTight = si.tightenBoundsFromSingletonRows(nFixed);

  int e0 = nErrors;
  checkInt(name, "nTightened", nTight, 1);
  check   (name, "col0Ub",     si.getColUpper()[0], 3.0);
  check   (name, "col1Ub",     si.getColUpper()[1], 10.0);
  if (nErrors == e0)
    printf("PASS %s\n", name);
  return nErrors - e0;
}

int main()
{
  testPositiveCoefUpper();
  testPositiveCoefLower();
  testPositiveCoefScale();
  testNegativeCoef();
  testRangedRowBothBounds();
  testMultipleSingletonsSameVar();
  testContradictoryRows();
  testNearZeroCoef();
  testNonSingletonIgnored();
  testBecomesFixed();
  testAlreadyFixed();
  testNoTightening();
  testMixedRows();

  printf("\n%d tests, %d failures\n", nTests, nErrors);
  if (nErrors) {
    fprintf(stderr, "%d test(s) FAILED\n", nErrors);
    return 1;
  }
  printf("All singleton-bound tightening tests passed.\n");
  return 0;
}
