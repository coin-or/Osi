// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiSolverParameters_H
#define OsiSolverParameters_H

enum OsiIntParam {
  /** The maximum number of iterations (whatever that means for the given
      solver) the solver can execute in the OsiSolverinterface::initialSolve()
      and the OsiSolverinterface::resolve() methods before terminating. */
  OsiMaxNumIteration = 0,
  /** The maximum number of iterations (whatever that means for the given
      solver) the solver can execute in the
      OsiSolverinterface::solveFromHotStart() method before terminating. */
  OsiMaxNumIterationHotStart,
  /** Just a marker, so that OsiSolverInterface can allocate a static sized
      array to store parameters. */
  OsiLastIntParam
};

enum OsiDblParam {
  /** Set Dual objective limit. This is to be used as a termination
      criteria in methods where the dual objective monotonically changes
      (e.g., dual simplex, the volume algorithm) */
  OsiDualObjectiveLimit = 0,
  /** Primal objective limit. This is to be used as a termination
      criteria in methods where the primal objective monotonically changes
      (e.g., primal simplex) */
  OsiPrimalObjectiveLimit,
  /** The maximum amount the dual constraints can be violated and still be
      considered feasible. */
  OsiDualTolerance,
  /** The maximum amount the primal constraints can be violated and still be
      considered feasible. */
  OsiPrimalTolerance,
  /** Just a marker, so that OsiSolverInterface can allocate a static sized
      array to store parameters. */
  OsiLastDblParam
};

#endif
