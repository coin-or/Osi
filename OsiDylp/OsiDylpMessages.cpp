#ifdef COIN_USE_DYLP

/*! \legal
  Copyright (C) 2003.
  Lou Hafer, International Business Machines Corporation and others. All
  Rights Reserved.
*/

namespace {
  char sccsid[] = "%W%	%G%" ;
  char cvsid[] = "$Id$" ;
}

#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpMessages.hpp"

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
  { ODSI_TEST_MSG, 0001, 2, "This is the us_en test message, eh." },
  { ODSI_MPSFILEIO, 0010, 2, "MPS file %s %s with %d errors." },
  { ODSI_UNSUPFORCEDO, 6001, 1, "Attempt to force unsupported hint; %s." },
  { ODSI_IGNORED, 3001, 2, "Ignored unsupported hint; %s." },
  { ODSI_DUMMY_END, 999999, 0, "" }
} ;

static MsgDefn uk_en_defns[] = {
  { ODSI_TEST_MSG, 0001, 2, "Blimey, this must be the uk_en test message" },
  { ODSI_DUMMY_END, 999999, 0, "" }
} ;

/*
  We seem to need a dummy CoinMessages object to prevent the compiler from
  complaining that CoinMessages::Language is unintialised.
*/
  const CoinMessages dummy(0) ;
/*
  The author is Canadian, eh? But we'll go with us_en anyways.
*/
  const CoinMessages::Language default_language = CoinMessages::us_en ;

} /* End unnamed:file */

/*
  This function constructs a CoinMessage object filled with a default
  set of messages, overlaid with whatever is available for the specified
  language.
*/

CoinMessages ODSI::setOsiDylpMessages (CoinMessages::Language local_language)

  
{ CoinMessages *msgs =
	new CoinMessages(sizeof(us_en_defns)/sizeof(MsgDefn)) ;

  msgs->setLanguage(local_language) ;
  strcpy(msgs->source_,"dylp");

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
    msgs->addMessage(msgdefn->inID,coinmsg) ;
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
    { msgs->replaceMessage(msgdefn->inID,msgdefn->fmt) ;
      msgdefn++ ; } }

  return (*msgs) ; }


#endif /* COIN_USE_DYLP */
