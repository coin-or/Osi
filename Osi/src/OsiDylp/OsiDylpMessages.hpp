#ifndef OsiDylpMessages_H
#define OsiDylpMessages_H

/*! \legal
  Copyright (C) 2004 -- 2007
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.
*/

/*
  sccs: @(#)OsiDylpMessages.hpp	1.5     09/16/04
  cvs: $Id$
*/


#include "CoinMessageHandler.hpp"

/*
  Enum used to specify ODSI messages to the message handler. There is no
  need for the order here to match the order of message definition in
  OsiDylpMessages.cpp, but all enum values must be here. ODSI_DUMMY_END must
  be last, however.
*/

typedef enum { ODSI_TEST_MSG,
	       ODSI_MPSFILEIO,
	       ODSI_UNSUPFORCEDO,
	       ODSI_IGNOREDHINT,
	       ODSI_EMPTYODWSB,
	       ODSI_NOTODWSB,
	       ODSI_ODWSBBADSIZE,
	       ODSI_ODWSBBADSTATUS,
	       ODSI_ODWSBSHORTBASIS,
	       ODSI_PRESOL_STATS,
	       ODSI_PRESOL_PASS,
	       ODSI_POSTSOL,
	       ODSI_POSTSOL_ACT,
	       ODSI_COLD,
	       ODSI_WARM,
	       ODSI_HOT,
	       ODSI_ALLDYLP,
	       ODSI_ATTACH,
	       ODSI_DETACH,
	       ODSI_NOSOLVE,
	       ODSI_FAILEDCALL,
	       ODSI_ACCESS_STALE,
	       ODSI_SHORTSTATS,
	       ODSI_CONFUSION,
	       ODSI_DUMMY_END } OsiDylpMessageID_enum ;

#endif /* OsiDylpMessages_H */
