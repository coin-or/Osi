/*! \legal
  Copyright (C) 2004 -- 2007
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.

  This file is a portion of the COIN/OSI interface for dylp.
*/

#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpMessages.hpp"

namespace {
  char sccsid[] UNUSED = "@(#)OsiDylpMessages.cpp	1.6	09/25/04" ;
  char cvsid[] UNUSED = "$Id$" ;
}

/* Cut name lengths for readability. */

#define ODSI OsiDylpSolverInterface
#define OSI OsiSolverInterface

/*
  The general strategy for setting up the COIN messaging facility is as
  follows:

    * Create a CoinMessages object, sized to accommodate the number of messages
      you plan to define. The initial size can be zero, but incremental growth
      is not particularly efficient.

    * Define your messages. There are three pieces of information that need
      to be supplied for each message: an external ID number, a level, and a
      format string. There is a notion of severity wired into the external ID
      number; see CoinMessageHandler for details. Level controls when the
      message is printed; CoinMessageHandler associates higher numbers with
      increased verbosity.  The format string can be null, if you want to
      build the whole thing on the fly.

      Severity, as of March, 2004:
	  0 - 2999	informational
       3000 - 5999	warnings
       6000 - 8999	non-fatal errors
       9000 -	 	fatal errors

    * Load the messages into the CoinMessage object. This involves repeated
      creation of a CoinOneMessage object, which is then loaded into the
      CoinMessages object using CoinMessages::addMessage. You specify an
      internal ID number at the time you add the message to the CoinMessages
      object.
    
      One approach to this is to define a derived class, XXXMessages, whose
      sole reason for existence is a constructor that performs the load. The
      internal ID numbers can be defined with an enum, which makes it a bit
      easier to specify a particular message and makes your code a bit more
      resilient against the inevitable insertion/deletion of messages as the
      code evolves.

    * Create a CoinMessageHandler object.

  To print a message, start construction of the message with

    CoinMessageHandler::message(<internal ID number>,<CoinMessages object>)
  
  Then use the various overloaded << operators to supply the necessary
  parameters. Finish with CoinMessageHandler::finish() or '<< CoinMessageEol'
  (where CoinMessageEol is part of the CoinMessageMarker enum).
*/

/*
  Message definitions.

  The precise form isn't important here, so long as the method that loads them
  into a CoinMessages object can come up with values for external number,
  detail level, and format string. The use of an enum to provide a local ID
  for each message is mainly useful with the internationalisation feature. It
  makes it easy to slap the same ID on alternate versions of a message.
*/

namespace { /* unnamed:file */

typedef struct { OsiDylpMessageID_enum inID ;
		 int exID ;
		 int lvl ;
		 const char *fmt ; } MsgDefn ;

static MsgDefn us_en_defns[] = {
  // informational (0 -- 2999)
  { ODSI_TEST_MSG, 1, 2, "This is the us_en test message, eh." },
  { ODSI_MPSFILEIO, 10, 5, "MPS file %s %s with %d errors." },
  { ODSI_COLD, 50, 3, "dylp cold start%? (%s)%? result %s, z = %g, iters = %d." },
  { ODSI_WARM, 51, 3, "dylp warm start result %s, z = %g, iters = %d." },
  { ODSI_HOT, 52, 3, "dylp hot start result %s, z = %g, iters = %d." },
  { ODSI_ALLDYLP, 53, 4, "dylp %s start odsi object %#x." },
  { ODSI_DETACH, 54, 5, "dylp detach from %#x." },
  { ODSI_ATTACH, 55, 5, "dylp attach in %s by %#x." },
  { ODSI_SHORTSTATS, 56, 4,
    "Problem %s: tot %.4f pre %.4f lp1 %d %.4f post %.4f lp2 %d %.4f" },
  { ODSI_PRESOL_STATS, 100, 3,
    "%s %d constraints, %d variables, %d coefficients." }, 
  { ODSI_PRESOL_PASS, 101, 6,
    "Presolve pass %d: dropped %d constraints (%.2f), %d variables (%.2f)." },
  { ODSI_POSTSOL, 200, 3, "Postsolve %s."},
  { ODSI_POSTSOL_ACT, 201, 6, "Applying postsolve transform %s."},
  // warning (3000 -- 5999)
  { ODSI_IGNOREDHINT, 3001, 2, "Ignored unsupported hint; %s." },
  { ODSI_ODWSBSHORTBASIS, 3100, 1,
    "[%s]: basis has only %d variables for %d constraints." },
  { ODSI_NOSOLVE, 3200, 1, "Impossible to call dylp; %s." },
  // Non-fatal errors (6000 -- 8999)
  { ODSI_UNSUPFORCEDO, 6001, 1, "Attempt to force unsupported hint; %s." },
  { ODSI_ACCESS_STALE, 6050, 1,
    "(%s) request to return a value from a stale solution."},
  { ODSI_EMPTYODWSB, 6101, 1, "Empty warm start basis object." },
  { ODSI_NOTODWSB, 6102, 1,
    "The warm start basis object is not a %sWarmStartBasis object." },
  { ODSI_ODWSBBADSIZE, 6103, 1,
    "Basis size %d x %d does not match constraint system size %d x %d." },
  { ODSI_ODWSBBADSTATUS, 6104, 1,
    "Flipping %s (%d) from %s to %s; lack of finite bound." },
  // Fatal errors (9000 and up)
  { ODSI_CONFUSION, 9001, 1,
    "Internal confusion, line %d." },
  { ODSI_DUMMY_END, 999999, 0, "" }
} ;

static MsgDefn uk_en_defns[] = {
  { ODSI_TEST_MSG, 1, 2, "Blimey, this must be the uk_en test message" },
  { ODSI_DUMMY_END, 999999, 0, "" }
} ;


/*
  We seem to need a dummy CoinMessages object to prevent the compiler from
  complaining that CoinMessages::Language is unintialised.
*/
  const CoinMessages dummy(0) ;
/*
  The author is Canadian, eh. But we'll go with us_en anyways.
*/
  const CoinMessages::Language default_language = CoinMessages::us_en ;

} /* End unnamed:file */

/*!
  This function constructs a CoinMessages object filled with a default set of
  messages, overlaid with whatever is available for the specified language.
  It is used to establish the initial set of messages, and is also called
  whenever the language is changed. The latter, because there's no way of
  guaranteeing that the sets of messages for alternate languages will all
  replace the same messages. This approach guarantees that the set of
  messages is always composed of the default language and the messages for
  one alternate language.
*/

void ODSI::setOsiDylpMessages (CoinMessages::Language local_language)

  
{ CoinMessages odsiMessages(sizeof(us_en_defns)/sizeof(MsgDefn)) ;

  odsiMessages.setLanguage(local_language) ;
  strcpy(odsiMessages.source_,"dylp");

/*
  Yes, this is gloriously redundant, but it's set up in anticipation of
  future extensions.
*/
  MsgDefn *msgdefn ;
  switch (default_language)
  { case CoinMessages::us_en:
    { msgdefn = us_en_defns ;
      break ; }
    default:
    { msgdefn = us_en_defns ;
      break ; } }
/*
  Open a loop to create and load the messages.
*/
  while (msgdefn->inID != ODSI_DUMMY_END)
  { CoinOneMessage coinmsg(msgdefn->exID,msgdefn->lvl,msgdefn->fmt) ;
    odsiMessages.addMessage(msgdefn->inID,coinmsg) ;
    msgdefn++ ; }
/*
  Now, if the local language differs from the default language, load any
  overrides. There's a trivial uk_en message set, solely for the purpose of
  testing language change routines.
*/
  if (local_language != default_language)
  { switch (local_language)
    { case CoinMessages::us_en:
      { msgdefn = us_en_defns ;
	break; }
      case CoinMessages::uk_en:
      { msgdefn = uk_en_defns ;
	break; }
      default:
      { msgdefn = us_en_defns ;
	break; } }

    while (msgdefn->inID != ODSI_DUMMY_END)
    { odsiMessages.replaceMessage(msgdefn->inID,msgdefn->fmt) ;
      msgdefn++ ; } }
/*
  Each CoinOneMessage has a fixed-length array to hold the message; by default
  this is 400 chars. Convert to `compressed' CoinOneMessage objects where the
  array is only as large as necessary. Any attempt to replace a message, or the
  message text, will automatically trigger a decompress operation before doing
  the replacement, but the messages will *not* be automatically recompressed.
*/
  odsiMessages.toCompact() ;

  messages_ = odsiMessages ;
  return ; }
