// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef OsiWarmStart_H
#define OsiWarmStart_H

//#############################################################################

/** Warmstart information abstract base class. <br>
    Really nothing can be generalized for warmstarting information. All we
    know that it exists. Hence the abstract base class contains only a virtual
    destructor. */

class OsiWarmStart {
public:
  virtual ~OsiWarmStart() {}
};

#endif
