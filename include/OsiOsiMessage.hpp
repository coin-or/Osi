// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiOsiMessage_H
#define OsiOsiMessage_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

/** This deals with Osi messages (as against Clp messages etc).
    OsiMessageHandler.hpp is the general part of message handling.
    All it has are enum's for the various messages.
    OsiOsiMessage.cpp has text in various languages.

    It is trivial to use the .hpp and .cpp file as a basis for
    messages for other components.
 */

#include "OsiMessageHandler.hpp"
enum OSI_Message
{
  OSI_CONTINUATION=-1,
  OSI_MPS_LINE=0,
  OSI_MPS_STATS,
  OSI_BAB_MAXITS,
  OSI_BAB_SOLUTION,
  OSI_BAB_END,
  OSI_BAB_INFEAS,
  OSI_BAB_STRONG,
  OSI_MPS_ILLEGAL,
  OSI_MPS_BADIMAGE,
  OSI_MPS_DUPOBJ,
  OSI_MPS_DUPROW,
  OSI_MPS_NOMATCHROW,
  OSI_MPS_NOMATCHCOL,
  OSI_BAB_NOINT,
  OSI_MPS_FILE,
  OSI_MPS_BADFILE1,
  OSI_MPS_BADFILE2,
  OSI_MPS_EOF,
  OSI_MPS_RETURNING,
  OSI_DUMMY_END
};

class OsiOsiMessage : public OsiMessages {

public:

  /**@name Constructors etc */
  //@{
  /** Constructor */
  OsiOsiMessage(Language language=us_en);
  //@}

};

#endif
