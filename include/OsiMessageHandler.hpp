// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef OsiMessageHandler_H
#define OsiMessageHandler_H

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <iostream>
#include <stdio.h>
#include <string>

/** This is a first attempt at a message handler.

 I am strongly in favo(u)r of language support for an open
 source project, but I have tried to make it as lightweight as
 possible in the sense that only a subset of messages need to be
 defined - the rest default to US English.  There will be different
 sets of messages for each component - so at present there is a
 Clp component and an Osi component.

 The default handler at present just prints to stdout or to a FILE pointer

 Because messages are only used in a controlled environment and have no
 impact on code and are tested by other tests I have include tests such
 as language and derivation in other unit tests.

 I need help with iso character codes as I don't understand how
 it all works.

*/

/** Class for one massaged message.
    At this stage it is in the correct language.
    A message may have several parts of which all except the
    first may be omitted.
 */

class OsiOneMessage {

public:
  /**@name Constructors etc */
  //@{
  /** Default constructor. */
  OsiOneMessage();
  /** Normal constructor */
  OsiOneMessage(int externalNumber, char detail,
		int numberParts, const char ** message);
  /** Destructor */
  ~OsiOneMessage();
  /** The copy constructor */
  OsiOneMessage(const OsiOneMessage&);
  /** assignment operator. */
  OsiOneMessage& operator=(const OsiOneMessage&);
  //@}

  /**@name Useful stuff */
  //@{
  /// Replaces messages (i.e. a different language)
  void replaceMessage(int numberParts, const char ** message);
  /// Returns part of message
  char * part(int iPart) const
  {return message_[iPart];};
  //@}

   /**@name gets and sets methods */
   //@{
  /// number to print out (also determines severity)
    inline int externalNumber() const
  {return externalNumber_;};
    inline void setExternalNumber(int number) 
  {externalNumber_=number;};
  /// Severity
  inline char severity() const
  {return severity_;};
  //@}

  /**@name member data */
  //@{
  // maximum number of parts to a message
#define MAXMESSAGE 4
  /// number to print out (also determines severity)
    int externalNumber_;
  /// Will only print if detail matches
    char detail_;
  /// Severity
    char severity_;
  /// Messages (in correct language)
    char * message_[MAXMESSAGE];
  //@}
};

/** Class for massaged messages.
    By now knows language and source
 */

class OsiMessages {

public:
  /** These are the languages that are supported.  At present only
      us_en is serious and the rest are for testing.  I can't
      really see uk_en being used seriously!
  */
  enum Language {
    us_en = 0,
    uk_en,
    it
  };
  /**@name Constructors etc */
  //@{
  /** Constructor with number of messages. */
  OsiMessages(int numberMessages=0);
  /** Destructor */
  ~OsiMessages();
  /** The copy constructor */
  OsiMessages(const OsiMessages&);
  /** assignment operator. */
  OsiMessages& operator=(const OsiMessages&);
  //@}

  /**@name Useful stuff */
  //@{
  /// Puts message in correct place
  void addMessage(int messageNumber, const OsiOneMessage & message);
  /// Replaces messages (i.e. a different language)
  void replaceMessage(int messageNumber, 
		int numberParts, const char ** message);
  /** Language.  Need to think about iso codes */
  inline Language language() const
  {return language_;};
  void setLanguage(Language language)
  {language_ = language;};
  //@}

  /**@name member data */
  //@{
  /// Message number
  int numberMessages_;
  /// Language 
  Language language_;
  /// Source e.g. Clp
  char source_[4];
  /// Messages
  OsiOneMessage ** message_;
  //@}
};
// for convenience eol
enum OsiMessageMarker {
  OsiMessageEol = 0,
  OsiMessageNewline = 1
};

/** Base class for message handling

    The default behavior is just to print. Inherit to change that
    This is a first attempt and it shows.
    messages <3000 are informational
             <6000 warnings
	     <9000 errors - but will try and carry on or return code
	     >=9000 stops
    Where there are derived classes I have started message numbers
    at 1001

    Messages can be printed with or without a prefix e.g. Clp0102I.
    A prefix makes the messages look less nimble but is very useful
    for "grep" etc.
*/

class OsiMessageHandler  {
  
public:
   /**@name Virtual methods that the derived classes may provide */
   //@{
  /** Print message, return 0 normally.
      At present this is the only virtual method
  */
   virtual int print() ;
   //@}

  /**@name Constructors etc */
  //@{
  /// Constructor
  OsiMessageHandler();
  /// Constructor to put to file pointer (won't be closed)
  OsiMessageHandler(FILE *fp);
  /** Destructor */
  virtual ~OsiMessageHandler();
  /** The copy constructor */
  OsiMessageHandler(const OsiMessageHandler&);
  /** assignment operator. */
  OsiMessageHandler& operator=(const OsiMessageHandler&);
  /// Clone
  virtual OsiMessageHandler * clone() const;
  //@}
   /**@name gets and sets methods */
   //@{
  /** Amount of print out:
      0 - none
      1 - minimal
      2 - normal low
      3 - normal high
      4 - verbose
      above that 8,16,32 etc just for selective debug and are for
      printf messages in code
  */
  inline int logLevel() const
          { return logLevel_;};
  void setLogLevel(int value);
  /// Switch on or off prefix
  void setPrefix(bool yesNo);
  /// values in message
  inline double doubleValue(int position) const
  { return doubleValue_[position];};
  inline int numberDoubleFields() const
  {return numberDoubleFields_;};
  inline int intValue(int position) const
  { return intValue_[position];};
  inline int numberIntFields() const
  {return numberIntFields_;};
  inline char charValue(int position) const
  { return charValue_[position];};
  inline int numberCharFields() const
  {return numberCharFields_;};
  inline std::string stringValue(int position) const
  { return stringValue_[position];};
  inline int numberStringFields() const
  {return numberStringFields_;};
  /// Current message
  inline OsiOneMessage  currentMessage() const
  {return currentMessage_;};
  /// Source of current message
  inline std::string currentSource() const
  {return source_;};
  /// Output buffer
  inline const char * messageBuffer() const
  {return messageBuffer_;};
  /// Highest message number (indicates any errors)
  inline int highestNumber() const
  {return highestNumber_;};
  //@}
  
  /**@name Actions to create a message  */
  //@{
  /// Start a message
  OsiMessageHandler & message(int messageNumber,
			      const OsiMessages & messages);
  /// returns current (just to make more uniform in look)
  OsiMessageHandler & message();
  /** Gets position of next field in format
      if initial then copies until first % */
  char * nextPerCent(char * start , const bool initial=false);
  /// Adds into message
  OsiMessageHandler & operator<< (int intvalue);
  OsiMessageHandler & operator<< (double doublevalue);
  OsiMessageHandler & operator<< (std::string stringvalue);
  OsiMessageHandler & operator<< (char charvalue);
  OsiMessageHandler & operator<< (const char *stringvalue);
  OsiMessageHandler & operator<< (OsiMessageMarker);
  /** Stop (and print) - or use eol
  */
  int finish();
  /// Positions to part n
  OsiMessageHandler & messageSection(int part);

  /** The following is to help existing codes interface
      Starts message giving number and complete text
  */
  OsiMessageHandler & message(int externalNumber,const char * header,
			      const char * msg,char severity);
  
  //@}
  
private:
  /**@name Private member data */
  //@{
  /// values in message
  double doubleValue_[10];
  int intValue_[10];
  char charValue_[10];
  std::string stringValue_[10];
  /// Log level
  int logLevel_;
  /// Whether we want prefix (may get more subtle so is int)
  int prefix_;
  /// Current message
  OsiOneMessage  currentMessage_;
  /// Part number
  int part_;
  /// Format string for message (remainder)
  char * format_;
  /// Number fields filled in,  0 in constructor then incremented
  int numberDoubleFields_;
  int numberIntFields_;
  int numberCharFields_;
  int numberStringFields_;
  /// Output buffer
  char messageBuffer_[1000];
  /// Position in output buffer
  char * messageOut_;
  /// Current source of message
  std::string source_;
  /// If we are skipping this message
  bool skip_;
  /// Highest message number (indicates any errors)
  int highestNumber_;
  /// File pointer
  FILE * fp_;
   //@}
};

#endif
