// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiOsiMessage.hpp"

typedef struct {
  OSI_Message internalNumber;
  int externalNumber; // or continuation
  char detail;
  const char * message;
} Osi_message;
static Osi_message us_english[]=
{
  {OSI_MPS_LINE,1,1,"At line %d %s"},
  {OSI_MPS_STATS,2,1,"Problem %s has %d rows, %d columns and %d elements"},
  {OSI_BAB_MAXITS,3,1,"Exiting on maximum iterations"},
  {OSI_BAB_SOLUTION,4,1,"Integer solution of %g found after %d iterations and %d nodes"},
  {OSI_BAB_END,5,1,"Search took %d iterations and %d nodes"},
  {OSI_BAB_INFEAS,6,1,"The LP relaxation is infeasible or too expensive"},
  {OSI_BAB_STRONG,7,3,"Strong branching on %d (%d), down %g up %g value %g"},
  {OSI_MPS_ILLEGAL,3001,0,"Illegal value for %s of %g"},
  {OSI_MPS_BADIMAGE,3002,0,"Bad image at line %d < %s >"},
  {OSI_MPS_DUPOBJ,3003,0,"Duplicate objective at line %d < %s >"},
  {OSI_MPS_DUPROW,3004,0,"Duplicate row %s at line %d < %s >"},
  {OSI_MPS_NOMATCHROW,3005,0,"No match for row %s at line %d < %s >"},
  {OSI_MPS_NOMATCHCOL,3006,0,"No match for column %s at line %d < %s >"},
  {OSI_BAB_NOINT,3007,0,"No integer variables - nothing to do"},
  {OSI_MPS_FILE,6001,0,"Unable to open mps input file %s"},
  {OSI_MPS_BADFILE1,6002,0,"Unknown image %s at line 1 of file %s"},
  {OSI_MPS_BADFILE2,6003,0,"Consider the possibility of a compressed file\
 which zlib is unable to read"},
  {OSI_MPS_EOF,6004,0,"EOF on file %s"},
  {OSI_MPS_RETURNING,6005,0,"Returning as too many errors"},
  {OSI_SOLVER_MPS,8,1,"%s read with %d errors"},
  {OSI_DUMMY_END,999999,0,""}
};
// **** aiutami!
static Osi_message italian[]=
{
  {OSI_MPS_LINE,1,1,"al numero %d %s"},
  {OSI_MPS_STATS,2,1,"matrice %s ha %d file, %d colonne and %d elementi (diverso da zero)"},
  {OSI_DUMMY_END,999999,0,""}
};
/* Constructor */
OsiOsiMessage::OsiOsiMessage(Language language) :
  OsiMessages(sizeof(us_english)/sizeof(Osi_message))
{
  language_=language;
  strcpy(source_,"Osi");
  Osi_message * message = us_english;

  while (message->internalNumber!=OSI_DUMMY_END) {
    OsiOneMessage oneMessage(message->externalNumber,message->detail,
		message->message);
    addMessage(message->internalNumber,oneMessage);
    message ++;
  }

  // now override any language ones

  switch (language) {
  case it:
    message = italian;
    break;

  default:
    message=NULL;
    break;
  }

  // replace if any found
  if (message) {
    while (message->internalNumber!=OSI_DUMMY_END) {
      replaceMessage(message->internalNumber,message->message);
      message ++;
    }
  }
}
