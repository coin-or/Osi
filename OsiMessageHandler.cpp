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
  message_=NULL;
  severity_='I';
  detail_=-1;
}
/* Destructor */
OsiOneMessage::~OsiOneMessage()
{
  free(message_);
}
/* The copy constructor */
OsiOneMessage::OsiOneMessage(const OsiOneMessage & rhs)
{
  externalNumber_=rhs.externalNumber_;
  if (rhs.message_)
    message_=strdup(rhs.message_);
  else
    message_=NULL;
  severity_=rhs.severity_;
  detail_=rhs.detail_;
}
/* assignment operator. */
OsiOneMessage& 
OsiOneMessage::operator=(const OsiOneMessage & rhs)
{
  if (this != &rhs) {
    externalNumber_=rhs.externalNumber_;
    free(message_);
    if (rhs.message_)
      message_=strdup(rhs.message_);
    else
      message_=NULL;
    severity_=rhs.severity_;
    detail_=rhs.detail_;
  }
  return *this;
}
/* Normal constructor */
OsiOneMessage::OsiOneMessage(int externalNumber, char detail,
		const char * message)
{
  externalNumber_=externalNumber;
  message_=strdup(message);
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
OsiOneMessage::replaceMessage( const char * message)
{
  free(message_);
  message_=strdup(message);
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
			     const char * message)
{
  assert(messageNumber<numberMessages_);
  message_[messageNumber]->replaceMessage(message);
}
// Changes detail level for one message
void 
OsiMessages::setDetailMessage(int newLevel, int messageNumber)
{
  int i;
  for (i=0;i<numberMessages_;i++) {
    if (message_[i]->externalNumber()==messageNumber) {
      message_[i]->setDetail(newLevel);
      break;
    }
  }
}
// Changes detail level for several messages
void 
OsiMessages::setDetailMessages(int newLevel, int numberMessages,
			       int * messageNumbers)
{
  int i;
  if (numberMessages<3) {
    // do one by one
    int j;
    for (j=0;j<numberMessages;j++) {
      int messageNumber = messageNumbers[j];
      for (i=0;i<numberMessages_;i++) {
	if (message_[i]->externalNumber()==messageNumber) {
	  message_[i]->setDetail(newLevel);
	  break;
	}
      }
    }
  } else {
    // do backward mapping
    int backward[10000];
    for (i=0;i<10000;i++) 
      backward[i]=-1;
    for (i=0;i<numberMessages_;i++) 
      backward[message_[i]->externalNumber()]=i;
    for (i=0;i<numberMessages;i++) {
      int iback = backward[messageNumbers[i]];
      if (iback>=0)
	message_[iback]->setDetail(newLevel);
    }
  }
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
      } 
    } 
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
  internalNumber_(0),
  format_(NULL),
  numberDoubleFields_(0),
  numberIntFields_(0),
  numberCharFields_(0),
  numberStringFields_(0),
  printStatus_(0),
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
  internalNumber_(0),
  format_(NULL),
  numberDoubleFields_(0),
  numberIntFields_(0),
  numberCharFields_(0),
  numberStringFields_(0),
  printStatus_(0),
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
  internalNumber_=rhs.internalNumber_;
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
  int offset = rhs.format_ - rhs.currentMessage_.message();
  format_ = currentMessage_.message()+offset;
  strcpy(messageBuffer_,rhs.messageBuffer_);
  offset = rhs.messageOut_-rhs.messageBuffer_;
  messageOut_= messageBuffer_+offset;
  printStatus_= rhs.printStatus_;
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
    internalNumber_=rhs.internalNumber_;
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
    int offset = rhs.format_ - rhs.currentMessage_.message();
    format_ = currentMessage_.message()+offset;
    strcpy(messageBuffer_,rhs.messageBuffer_);
    offset = rhs.messageOut_-rhs.messageBuffer_;
    messageOut_= messageBuffer_+offset;
    printStatus_= rhs.printStatus_;
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
  internalNumber_=messageNumber;
  currentMessage_= *(normalMessage.message_[messageNumber]);
  source_ = normalMessage.source_;
  format_ = currentMessage_.message_;
  messageBuffer_[0]='\0';
  messageOut_=messageBuffer_;
  highestNumber_ = max(highestNumber_,currentMessage_.externalNumber_);
  // do we print
  int detail = currentMessage_.detail_;
  printStatus_=0;
  if (detail>=8) {
    // bit setting - debug
    if ((detail&logLevel_)==0)
      printStatus_ = 3;
  } else if (logLevel_<detail) {
    printStatus_ = 3;
  }
  if (!printStatus_) {
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
  internalNumber_=externalNumber;
  currentMessage_= OsiOneMessage();
  currentMessage_.setExternalNumber(externalNumber);
  source_ = source;
  // mark so will not update buffer
  printStatus_=2;
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
  /* Allows for skipping printing of part of message,
      but putting in data */
  OsiMessageHandler & 
OsiMessageHandler::printing(bool onOff)
{
  // has no effect if skipping or whole message in
  if (printStatus_<2) {
    assert(format_[1]=='?');
    if (onOff)
      printStatus_=0;
    else
      printStatus_=1;
    format_ = nextPerCent(format_+2,true);
  }
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
  internalNumber_=-1;
  format_ = NULL;
  messageBuffer_[0]='\0';
  messageOut_=messageBuffer_;
  printStatus_=true;
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
	// %? is skipped over as it is just a separator
	if (nextPerCent[1]!='?') {
	  if (initial&&!printStatus_) {
	    int numberToCopy=nextPerCent-start;
	    strncpy(messageOut_,start,numberToCopy);
	    messageOut_+=numberToCopy;
	  } 
	  start=nextPerCent;
	  if (start[1]!='%') {
	    foundNext=true;
	    if (!initial) 
	      *start='\0'; //zap
	  } else {
	    start+=2;
	    if (initial) {
	      *messageOut_='%';
	      messageOut_++;
	    } 
	  }
	} else {
	  foundNext=true;
	  // skip to % and zap
	  start=nextPerCent;
	  *start='\0'; 
	}
      } else {
        if (initial&&!printStatus_) {
          strcpy(messageOut_,start);
          messageOut_+=strlen(messageOut_);
        } 
        start=0;
        foundNext=true;
      } 
    } 
  } 
  return start;
}
// Adds into message
OsiMessageHandler & 
OsiMessageHandler::operator<< (int intvalue)
{
  if (printStatus_==3)
    return *this; // not doing this message
  intValue_[numberIntFields_++] = intvalue;
  if (printStatus_<2) {
    if (format_) {
      //format is at % (but may be changed to null)
      *format_='%';
      char * next = nextPerCent(format_+1);
      // could check
      if (!printStatus_) {
	sprintf(messageOut_,format_,intvalue);
	messageOut_+=strlen(messageOut_);
      }
      format_=next;
    } else {
      sprintf(messageOut_," %d",intvalue);
      messageOut_+=strlen(messageOut_);
    } 
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (double doublevalue)
{
  if (printStatus_==3)
    return *this; // not doing this message
  doubleValue_[numberDoubleFields_++] = doublevalue;
  if (printStatus_<2) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      // could check
      if (!printStatus_) {
	sprintf(messageOut_,format_,doublevalue);
	messageOut_+=strlen(messageOut_);
      }
      format_=next;
    } else {
      sprintf(messageOut_," %g",doublevalue);
      messageOut_+=strlen(messageOut_);
    } 
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (std::string stringvalue)
{
  if (printStatus_==3)
    return *this; // not doing this message
  stringValue_[numberStringFields_++] = stringvalue;
  if (printStatus_<2) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      // could check
      if (!printStatus_) {
	sprintf(messageOut_,format_,stringvalue.c_str());
	messageOut_+=strlen(messageOut_);
      }
      format_=next;
    } else {
      sprintf(messageOut_," %s",stringvalue.c_str());
      messageOut_+=strlen(messageOut_);
    } 
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (char charvalue)
{
  if (printStatus_==3)
    return *this; // not doing this message
  intValue_[numberCharFields_++] = charvalue;
  if (printStatus_<2) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      // could check
      if (!printStatus_) {
	sprintf(messageOut_,format_,charvalue);
	messageOut_+=strlen(messageOut_);
      }
      format_=next;
    } else {
      sprintf(messageOut_," %c",charvalue);
      messageOut_+=strlen(messageOut_);
    } 
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (const char *stringvalue)
{
  if (printStatus_==3)
    return *this; // not doing this message
  stringValue_[numberStringFields_++] = stringvalue;
  if (printStatus_<2) {
    if (format_) {
      //format is at % (but changed to 0)
      *format_='%';
      char * next = nextPerCent(format_+1);
      // could check
      if (!printStatus_) {
	sprintf(messageOut_,format_,stringvalue);
	messageOut_+=strlen(messageOut_);
      }
      format_=next;
    } else {
      sprintf(messageOut_," %s",stringvalue);
      messageOut_+=strlen(messageOut_);
    } 
  }
  return *this;
}
OsiMessageHandler & 
OsiMessageHandler::operator<< (OsiMessageMarker marker)
{
  if (printStatus_!=3) {
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

// returns current 
OsiMessageHandler & 
OsiMessageHandler::message()
{
  return * this;
}
