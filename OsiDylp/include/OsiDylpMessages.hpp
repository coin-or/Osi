#ifdef COIN_USE_DYLP
#ifndef OsiDylpMessages_H
#define OsiDylpMessages_H

/*! \legal
  Copyright (C) 2004.
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.
*/

/*
  sccs: @(#)OsiDylpMessages.hpp	1.5     09/16/04
  cvs: $Id$
*/


#include "CoinMessageHandler.hpp"

typedef enum { ODSI_TEST_MSG,
	       ODSI_MPSFILEIO,
	       ODSI_UNSUPFORCEDO,
	       ODSI_IGNORED,
	       ODSI_EMPTYODWSB,
	       ODSI_NOTODWSB,
	       ODSI_ODWSBBADSIZE,
	       ODSI_ODWSBBADSTATUS,
	       ODSI_PRESOL_STATS,
	       ODSI_PRESOL_PASS,
	       ODSI_POSTSOL,
	       ODSI_POSTSOL_ACT,
	       ODSI_LPRESULT,
	       ODSI_DUMMY_END } OsiDylpMessageID_enum ;

#endif /* OsiDylpMessages_H */
#endif /* COIN_USE_DYLP */
