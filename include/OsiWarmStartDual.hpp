// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiWarmStartDual_H
#define OsiWarmStartDual_H

#include "CoinHelperFunctions.hpp"
#include "OsiWarmStart.hpp"

//#############################################################################

/** WarmStart information that is only a dual vector */

class OsiWarmStartDual : public OsiWarmStart {
private:
  OsiWarmStartDual(const OsiWarmStartDual&);
  OsiWarmStartDual& operator=(const OsiWarmStartDual&);

public:
  /// return the size of the dual vector
  int size() const { return dualSize_; }
  /// return a pointer to the array of duals
  const double * dual() const { return dualVector_; }

  /** Assign the dual vector to be the warmstart information. In this method
      the object assumes ownership of the pointer and upon return "dual" will
      be a NULL pointer. If copying is desirable use the constructor. */
  void assignDual(int size, double *& dual) {
    dualSize_ = size;
    delete[] dualVector_;
    dualVector_ = dual;
    dual = NULL;
  }

  OsiWarmStartDual() : dualSize_(0), dualVector_(NULL) {}
  OsiWarmStartDual(int size, const double * dual) :
    dualSize_(size), dualVector_(new double[size]) {
    CoinDisjointCopyN(dual, size, dualVector_);
  }
  virtual ~OsiWarmStartDual() {
    delete[] dualVector_;
  }

private:
  ///@name Private data members
  //@{
    /// the size of the dual vector
    int dualSize_;
    /// the dual vector itself
    double * dualVector_;
  //@}
};

#endif
