// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "OsiMessageHandler.hpp"
#include <cassert>
#include <map>
using std::min;
using std::max;

/* Default constructor. */
OsiOneMessage::OsiOneMessage()
{
  externalNumber_=-1;
  int i;
  for (i=0;i<MAXMESSAGE;i++) 
    message_[i]=NULL;
  severity_='I';
  detail_=128;
}
/* Destructor */
OsiOneMessage::~OsiOneMessage()
{
  int i;
  for (i=0;i<MAXMESSAGE;i++) 
    free(message_[i]);
}
/* The copy constructor */
OsiOneMessage::OsiOneMessage(const OsiOneMessage & rhs)
{
  externalNumber_=rhs.externalNumber_;
  int i;
  for (i=0;i<MAXMESSAGE;i++) {
    if (rhs.message_[i])
      message_[i]=strdup(rhs.message_[i]);
    else
      message_[i]=NULL;
  }
  severity_=rhs.severity_;
  detail_=rhs.detail_;
}
/* assignment operator. */
OsiOneMessage& 
OsiOneMessage::operator=(const OsiOneMessage & rhs)
{
  if (this != &rhs) {
    externalNumber_=rhs.externalNumber_;
    int i;
    for (i=0;i<MAXMESSAGE;i++) {
      free(message_[i]);
      if (rhs.message_[i])
	message_[i]=strdup(rhs.message_[i]);
      else
	message_[i]=NULL;
    }
    severity_=rhs.severity_;
    detail_=rhs.detail_;
  }
  return *this;
}
/* Normal constructor */
OsiOneMessage::OsiOneMessage(int externalNumber, char detail,
		int numberParts, const char ** message)
{
  externalNumber_=externalNumber;
  int i;
  for (i=0;i<MAXMESSAGE;i++) 
    message_[i]=NULL;

  for (i=0;i<numberParts;i++) 
    message_[i]=strdup(message[i]);
  if (externalNumber<3000)
    severity_='I';
  else if (externalNumber<6000)
    severity_='W';
  else if (externalNumber<9000)
    severity_='E';
  else
    severity_='S';
  detail_=detail;
}
// Replaces messages (i.e. a different language)
void 
OsiOneMessage::replaceMessage(int numberParts, const char ** message)
{
  int i;
  for (i=0;i<MAXMESSAGE;i++) 
    message_[i]=NULL;

  for (i=0;i<numberParts;i++) 
    message_[i]=strdup(message[i]);
}


/* Constructor with number of messages. */
OsiMessages::OsiMessages(int numberMessages)
{
  numberMessages_=numberMessages;
  language_=us_en;
  strcpy(source_,"Unk");
  if (numberMessages_) {
    message_ = new OsiOneMessage * [numberMessages_];
    int i;
    for (i=0;i<numberMessages_;i++) 
      message_[i]=NULL;
  } else {
    message_=NULL;
  }
}
/* Destructor */
OsiMessages::~OsiMessages()
{
  int i;
  for (i=0;i<numberMessages_;i++) 
    delete message_[i];
  delete [] message_;
}
/* The copy constructor */
OsiMessages::OsiMessages(const OsiMessages & rhs)
{
  numberMessages_=rhs.numberMessages_;
  language_=rhs.language_;
  strcpy(source_,rhs.source_);
  if (numberMessages_) {
    message_ = new OsiOneMessage * [numberMessages_];
    int i;
    for (i=0;i<numberMessages_;i++) 
      if (rhs.message_[i])
	message_[i]=new OsiOneMessage(*(rhs.message_[i]));
      else
	message_[i] = NULL;
  } else {
    message_=NULL;
  }
}
/* assignment operator. */
OsiMessages& 
OsiMessages::operator=(const OsiMessages & rhs)
{
  if (this != &rhs) {
    language_=rhs.language_;
    strcpy(source_,rhs.source_);
    int i;
    for (i=0;i<numberMessages_;i++)
	delete message_[i];
    delete [] message_;
    numberMessages_=rhs.numberMessages_;
    if (numberMessages_) {
      message_ = new OsiOneMessage * [numberMessages_];
      int i;
      for (i=0;i<numberMessages_;i++) 
	if (rhs.message_[i])
	  message_[i]=new OsiOneMessage(*(rhs.message_[i]));
	else
	  message_[i] = NULL;
    } else {
      message_=NULL;
    }
  }
  return *this;
}
// Puts message in correct place
void 
OsiMessages::addMessage(int messageNumber, const OsiOneMessage & message)
{
  if (messageNumber>=numberMessages_) {
    // should not happen but allow for it
    OsiOneMessage ** temp = new OsiOneMessage * [messageNumber+1];
    int i;
    for (i=0;i<numberMessages_;i++) 
      temp[i] = message_[i];
    for (;i<=messageNumber;i++) 
      temp[i] = NULL;
    delete [] message_;
    message_ = temp;
  }
  delete message_[messageNumber];
  message_[messageNumber]=new OsiOneMessage(message);
}
// Replaces messages (i.e. a different language)
void 
OsiMessages::replaceMessage(int messageNumber, 
			    int numberParts, const char ** message)
{
  assert(messageNumber<numberMessages_);
  message_[messageNumber]->replaceMessage(numberParts,message);
}

// Print message, return 0 normally
int 
OsiMessageHandler::print() 
{
  if (messageOut_>messageBuffer_) {
    *messageOut_=0;
    //take off trailing spaces and commas
    messageOut_--;
    while (messageOut_>=messageBuffer_) {
      if (*messageOut_==' '||*messageOut_==',') {
	*messageOut_=0;
	messageOut_--;
      } else {
	break;
      } /* endif */
    } /* endwhile */
    fprintf(fp_,"%s\n",messageBuffer_);
    if (currentMessage_.severity_=='S') {
      fprintf(fp_,"Stopping due to previous errors.\n");
      //Should do walkback
      abort();
    } 
  }
  return 0;
}
/* Amount of print out:
   0 - none
   1 - minimal
   2 - normal low
   3 - normal high
   4 - verbose
   above that 8,16,32 etc just for selective debug and are for
   printf messages in code
*/
void 
OsiMessageHandler::setLogLevel(int value)
{
  if (value>=0)
    logLevel_=value;
}
void 
OsiMessageHandler::setPrefix(bool value)
{
  if (value)
    prefix_ = 255;
  else
    prefix_ =0;
}
// Constructor
OsiMessageHandler::OsiMessageHandler() :
  logLevel_(1),
  prefix_(255),
  currentMessage_(),
  part_(0),
  format_(NULL),
  numberDoubleFields_(0),
  numberIntFields_(0),
  numberCharFields_(0),
  numberStringFields_(0),
  skip_(true),
  highestNumber_(-1),
  fp_(stdout)
{
  messageBuffer_[0]='\0';
  messageOut_ = messageBuffer_;
  source_="Unk";
}
// Constructor
OsiMessageHandler::OsiMessageHandler(FILE * fp) :
  logLevel_(1),
  prefix_(255),
  currentMessage_(),
  part_(0),
  format_(NULL),
  numberDoubleFields_(0),
  numberIntFields_(0),
  numberCharFields_(0),
  numberStringFields_(0),
  skip_(true),
  highestNumber_(-1),
  fp_(fp)
{
  messageBuffer_[0]='\0';
  messageOut_ = messageBuffer_;
  source_="Unk";
}
/* Destructor */
OsiMessageHandler::~OsiMessageHandler()
{
}
/* The copy constructor */
OsiMessageHandler::OsiMessageHandler(const OsiMessageHandler& rhs)
{
  logLevel_=rhs.logLevel_;
  prefix_ = rhs.prefix_;
  currentMessage_=rhs.currentMessage_;
  part_=rhs.part_;
  int i;
  numberDoubleFields_ = rhs.numberDoubleFields_;
  for (i=0;i<numberDoubleFields_;i++) 
    doubleValue_[i]=rhs.doubleValue_[i];
  numberIntFields_ = rhs.numberIntFields_;
  for (i=0;i<numberIntFields_;i++) 
    intValue_[i]=rhs.intValue_[i];
  numberCharFields_ = rhs.numberCharFields_;
  for (i=0;i<numberCharFields_;i++) 
    charValue_[i]=rhs.charValue_[i];
  numberStringFields_ = rhs.numberStringFields_;
  for (i=0;i<numberStringFields_;i++) 
    stringValue_[i]=rhs.stringValue_[i];
  int offset = rhs.format_ - rhs.currentMessage_.part(part_);
  format_ = currentMessage_.part(part_)+offset;
  strcpy(messageBuffer_,rhs.messageBuffer_);
  offset = rhs.messageOut_-rhs.messageBuffer_;
  messageOut_= messageBuffer_+offset;
  skip_= rhs.skip_;
  highestNumber_= rhs.highestNumber_;
  fp_ = rhs.fp_;
  source_ = rhs.source_;
}
/* assignment operator. */
OsiMessageHandler & 
OsiMessageHandler::operator=(const OsiMessageHandler& rhs)
{
  if (this != &rhs) {
    logLevel_=rhs.logLevel_;
    prefix_ = rhs.prefix_;
    currentMessage_=rhs.currentMessage_;
    part_=rhs.part_;
    int i;
    numberDoubleFields_ = rhs.numberDoubleFields_;
    for (i=0;i<numberDoubleFields_;i++) 
      doubleValue_[i]=rhs.doubleValue_[i];
    numberIntFields_ = rhs.numberIntFields_;
    for (i=0;i<numberIntFields_;i++) 
      intValue_[i]=rhs.intValue_[i];
    numberCharFields_ = rhs.numberCharFields_;
    for (i=0;i<numberCharFields_;i++) 
      charValue_[i]=rhs.charValue_[i];
    numberStringFields_ = rhs.numberStringFields_;
    for (i=0;i<numberStringFields_;i++) 
      stringValue_[i]=rhs.stringValue_[i];
    int offset = rhs.format_ - rhs.currentMessage_.part(part_);
    format_ = currentMessage_.part(part_)+offset;
    strcpy(messageBuffer_,rhs.messageBuffer_);
    offset = rhs.messageOut_-rhs.messageBuffer_;
    messageOut_= messageBuffer_+offset;
    skip_= rhs.skip_;
    highestNumber_= rhs.highestNumber_;
    fp_ = rhs.fp_;
    source_ = rhs.source_;
  }
  return *this;
}
// Clone
OsiMessageHandler * 
OsiMessageHandler::clone() const
{
  return new OsiMessageHandler(*this);
}
// Start a message
OsiMessageHandler & 
OsiMessageHandler::message(int messageNumber,
			   const OsiMessages &normalMessage)
{
  if (messageOut_!=messageBuffer_) {
    // put out last message
    print();
  }
  numberDoubleFields_=0;
  numberIntFields_=0;
  numberCharFields_=0;
  numberStringFields_=0;
  part_=0;
  currentMessage_= *(normalMessage.message_[messageNumber]);
  source_ = normalMessage.source_;
  format_ = currentMessage_.message_[0];
  messageBuffer_[0]='\0';
  messageOut_=messageBuffer_;
  highestNumber_ = max(highestNumber_,currentMessage_.externalNumber_);
  // do we print
  int detail = currentMessage_.detail_;
  skip_=true;
  if (detail>=8) {
    // bit setting - debug
    if ((detail&logLevel_)!=0)
      skip_ = false;
  } else if (logLevel_>=detail) {
    skip_ = false;
  }
  if (!skip_) {
    if (prefix_) {
      sprintf(messageOut_,"%s%4.4d%c ",source_.c_str(),
	      currentMessage_.externalNumber_,
	      currentMessage_.severity_);
      messageOut_ += strlen(messageOut_);
    }
    format_ = nextPerCent(format_,true);
  }
  return *this;
}
/* The following is to help existing codes interface
   Starts message giving number and complete text
*/
OsiMessageHandler & 
OsiMessageHandler::message(int externalNumber,const char * source,
			   const char * msg, char severity)
{
  if (messageOut_!=messageBuffer_) {
    // put out last message
    print();
  }
  numberDoubleFields_=0;
  numberIntFields_=0;
  numberCharFields_=0;
  numberStringFields_=0;
  // mark so will not update buffer
  part_=-1;
  currentMessage_= OsiOneMessage();
  currentMessage_.setExternalNumber(externalNumber);
  source_ = source;
  format_ = NULL;
  highestNumber_ = max(highestNumber_,externalNumber);
  // If we get here we always print
  if (prefix_) {
    sprintf(messageOut_,"%s%4.4d%c ",source_.c_str(),
	    externalNumber,
	    severity);
  }
  strcat(messageBuffer_,msg);
  messageOut_=messageBuffer_+strlen(messageBuffer_);
  return *this;
}
/* Stop (and print) 
*/
int 
OsiMessageHandler::finish()
{
  if (messageOut_!=messageBuffer_) {
    // put out last message
    print();
  }
  numberDoubleFields_=0;
  numberIntFields_=0;
  numberCharFields_=0;
  numberStringFields_=0;
  part_=0;
  format_ = NULL;
  messageBuffer_[0]='\0';
  messageOut_=messageBuffer_;
  skip_=true;
  return 0;
}
/* Gets position of next field in format
   if initial then copies until first % */
char * 
OsiMessageHandler::nextPerCent(char * start , const bool initial)
{
  if (start) {
    bool foundNext=false;
    while (!foundNext) {
      char * nextPerCent = strchr(start,'%');
      if (nextPerCent) {
        if (initial) {
          int numberToCopy=nextPerCent-start;
          strncpy(messageOut_,start,numberToCopy);
          messageOut_+=numberToCopy;
        } /* endif */
        start=nextPerCent;
        if (start[1]!='%') {
          foundNext=true;
        } else {
          start+=2;
          if (initial) {
            *messageOut_='%';
            messageOut_++;
          } /* endif */
        } /* endif */
      } else {
        if (initial) {
          strcpy(messageOut_,start);
          messageOut_+=strlen(messageOut_);
        } /* endif */
        start=0;
        foundNext=true;
      } /* endif */
    } /* endwhile */
  } /* endif */
  return start;
}
// Adds into message
OsiMessageHandler & 
OsiMessageHandler::operator<< (int intvalue)
{
  if (skip_)
    return *this; // not doing this message
  intValue_[numberIntFields_++] = intvalue;
  if (part_>=0) {
    if (format_) {
      //format is at % (but may be changed to null)
      *format_='%';
      char * next = nextPerCent(format_+1);
      if (next) 
	*next=0;
      // could check
      sprintf(messageOut_,format_,intvalue);
      format_=next;
    } else {
      sprintf(messageOut_," %d",intvalue);
    } /* endif */
    messageOut_+=strlen(messageOut_);
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (double doublevalue)
{
  if (skip_)
    return *this; // not doing this message
  doubleValue_[numberDoubleFields_++] = doublevalue;
  if (part_>=0) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      if (next) 
	*next=0;
      // could check
      sprintf(messageOut_,format_,doublevalue);
      format_=next;
    } else {
      sprintf(messageOut_," %g",doublevalue);
    } /* endif */
    messageOut_+=strlen(messageOut_);
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (std::string stringvalue)
{
  if (skip_)
    return *this; // not doing this message
  stringValue_[numberStringFields_++] = stringvalue;
  if (part_>=0) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      if (next) 
	*next=0;
      // could check
      sprintf(messageOut_,format_,stringvalue.c_str());
      format_=next;
    } else {
      sprintf(messageOut_," %s",stringvalue.c_str());
    } /* endif */
    messageOut_+=strlen(messageOut_);
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (char charvalue)
{
  if (skip_)
    return *this; // not doing this message
  intValue_[numberCharFields_++] = charvalue;
  if (part_>=0) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      if (next) 
	*next=0;
      // could check
      sprintf(messageOut_,format_,charvalue);
      format_=next;
    } else {
      sprintf(messageOut_," %c",charvalue);
    } /* endif */
    messageOut_+=strlen(messageOut_);
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (const char *stringvalue)
{
  if (skip_)
    return *this; // not doing this message
  stringValue_[numberStringFields_++] = stringvalue;
  if (part_>=0) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      if (next) 
	*next=0;
      // could check
      sprintf(messageOut_,format_,stringvalue);
      format_=next;
    } else {
      sprintf(messageOut_," %s",stringvalue);
    } /* endif */
    messageOut_+=strlen(messageOut_);
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (OsiMessageMarker marker)
{
  if (!skip_) {
    switch (marker) {
    case OsiMessageEol:
      finish();
      break;
    case OsiMessageNewline:
      strcat(messageOut_,"\n");
      messageOut_++;
      break;
    }
  } else {
    // skipping - tidy up
    format_ = NULL;
  }
  return *this;
}
// Positions to part n
OsiMessageHandler & 
OsiMessageHandler::messageSection(int part)
{
  if (skip_)
    return *this; // not doing this message
  part_=part;
  format_ = currentMessage_.message_[part];
  format_ = nextPerCent(format_,true);
  return *this;
}
// returns current 
OsiMessageHandler & 
OsiMessageHandler::message()
{
  return * this;
}
