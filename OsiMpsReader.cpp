// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiMpsReader.hpp"
const bool True = true;
const bool False = false;

#include <math.h>

/// The following lengths are in decreasing order (for 64 bit etc)
/// Large enough to contain element index
typedef int OSIElementIndex;

/// Large enough to contain column index
typedef int OSIColumnIndex;

/// Large enough to contain row index (or basis)
typedef int OSIRowIndex;

#define MAX_FIELD_LENGTH 100
#define MAX_CARD_LENGTH 5*MAX_FIELD_LENGTH+80

enum OSISectionType { OSI_NO_SECTION, OSI_NAME_SECTION, OSI_ROW_SECTION,
  OSI_COLUMN_SECTION,
  OSI_RHS_SECTION, OSI_RANGE_SECTION, OSI_BOUND_SECTION,
  OSI_ENDATA_SECTION, OSI_EOF_SECTION, OSI_UNKNOWN_SECTION
};

enum OSIMpsType { OSI_N_ROW, OSI_E_ROW, OSI_L_ROW, OSI_G_ROW,
  OSI_BLANK_COLUMN, OSI_S1_COLUMN, OSI_S2_COLUMN, OSI_S3_COLUMN,
  OSI_INTORG, OSI_INTEND, OSI_SOSEND, OSI_UNSET_BOUND,
  OSI_UP_BOUND, OSI_FX_BOUND, OSI_LO_BOUND, OSI_FR_BOUND,
  OSI_MI_BOUND, OSI_PL_BOUND, OSI_BV_BOUND, OSI_UI_BOUND,
  OSI_SC_BOUND, OSI_UNKNOWN_MPS_TYPE
};

/// Very simple code for reading MPS data
class OSIMpsio {

public:

  /**@name Constructor and destructor */
  //@{
  /// Constructor expects file to be open - reads down to (and reads) NAME card
  OSIMpsio ( FILE * fp );
  /// Destructor
  ~OSIMpsio (  );
  //@}


  /**@name card stuff */
  //@{
  /// Gets next field and returns section type e.g. OSI_COLUMN_SECTION
  OSISectionType nextField (  );
  /// Returns current section type
  inline OSISectionType whichSection (  ) const {
    return section_;
  };
  /// Only for first field on card otherwise BLANK_COLUMN
  /// e.g. OSI_E_ROW
  inline OSIMpsType mpsType (  ) const {
    return mpsType_;
  };
  /// Cleans card - taking out trailing blanks
  void cleanCard();
  /// Returns row name of current field
  inline const char *rowName (  ) const {
    return rowName_;
  };
  /// Returns column name of current field
  inline const char *columnName (  ) const {
    return columnName_;
  };
  /// Returns value in current field
  inline double value (  ) const {
    return value_;
  };
  /// Whole card (for printing)
  inline const char *card (  ) const {
    return card_;
  };
  /// Returns card number
  inline OSIElementIndex cardNumber (  ) const {
    return cardNumber_;
  };
  //@}

////////////////// data //////////////////
private:

  /**@name data */
  //@{
  /// Current card image
  char card_[MAX_CARD_LENGTH];
  /// Current position within card image
  char *position_;
  /// End of card
  char *eol_;
  /// Current OSIMpsType
  OSIMpsType mpsType_;
  /// Current row name
  char rowName_[MAX_FIELD_LENGTH];
  /// Current column name
  char columnName_[MAX_FIELD_LENGTH];
  /// Current value
  double value_;
  /// File pointer
  FILE *fp_;
  /// Which section we think we are in
  OSISectionType section_;
  /// Card number
  OSIElementIndex cardNumber_;
  /// Whether free format.  Just for blank RHS etc
  bool freeFormat_;
  /// If all names <= 8 characters then allow embedded blanks
  bool eightChar_;
  //@}
};
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

//#############################################################################
// sections
const static char *section[] = {
  "", "NAME", "ROW", "COLUMN", "RHS", "RANGE", "BOUND", "ENDATA", " "
};

// what is allowed in each section - must line up with OSISectionType
const static OSIMpsType startType[] = {
  OSI_UNKNOWN_MPS_TYPE, OSI_UNKNOWN_MPS_TYPE,
  OSI_N_ROW, OSI_BLANK_COLUMN,
  OSI_BLANK_COLUMN, OSI_BLANK_COLUMN,
  OSI_UP_BOUND, OSI_UNKNOWN_MPS_TYPE,
  OSI_UNKNOWN_MPS_TYPE, OSI_UNKNOWN_MPS_TYPE
};
const static OSIMpsType endType[] = {
  OSI_UNKNOWN_MPS_TYPE, OSI_UNKNOWN_MPS_TYPE,
  OSI_BLANK_COLUMN, OSI_UNSET_BOUND,
  OSI_S1_COLUMN, OSI_S1_COLUMN,
  OSI_UNKNOWN_MPS_TYPE, OSI_UNKNOWN_MPS_TYPE,
  OSI_UNKNOWN_MPS_TYPE, OSI_UNKNOWN_MPS_TYPE
};
const static int allowedLength[] = {
  0, 0,
  1, 2,
  0, 0,
  2, 0,
  0, 0
};

// names of types
const static char *mpsTypes[] = {
  "N", "E", "L", "G",
  "  ", "S1", "S2", "S3", "  ", "  ", "  ",
  "  ", "UP", "FX", "LO", "FR", "MI", "PL", "BV", "UI", "SC"
};

void OSIMpsio::cleanCard()
{
  cardNumber_++;
  char * lastNonBlank = card_-1;
  char * image = card_;
  while ( *image != '\0' ) {
    if ( *image != '\t' && *image < ' ' ) {
      break;
    } else if ( *image != '\t' && *image != ' ') {
      lastNonBlank = image;
    }
    image++;
  }
  *(lastNonBlank+1)='\0';
}

char *
nextBlankOr ( char *image )
{
  bool blank = false;

  while ( !blank ) {
    if ( *image == ' ' || *image == '\t' ) {
      break;
    }
    if ( *image == '\0' )
      return NULL;
    image++;
  }
  return image;
}

//  OSIMpsio.  Constructor
OSIMpsio::OSIMpsio (  FILE * fp )
{
  memset ( card_, 0, MAX_CARD_LENGTH );
  position_ = card_;
  eol_ = card_;
  mpsType_ = OSI_UNKNOWN_MPS_TYPE;
  memset ( rowName_, 0, MAX_FIELD_LENGTH );
  memset ( columnName_, 0, MAX_FIELD_LENGTH );
  value_ = 0.0;
  fp_ = fp;
  section_ = OSI_EOF_SECTION;
  cardNumber_ = 0;
  freeFormat_ = false;
  eightChar_ = true;
  bool found = false;

  while ( !found ) {
    // need new image
    char *getit = fgets ( card_, MAX_CARD_LENGTH, fp_ );

    if ( !getit ) {
      section_ = OSI_EOF_SECTION;
      break;
    }
    // strip off newline etc
    cleanCard();
    if ( !strncmp ( card_, "NAME", 4 ) ) {
      section_ = OSI_NAME_SECTION;
      char *next = card_ + 4;

      {
	std::cout<<"At line "<< cardNumber_ <<" "<< card_<<std::endl;
      }
      while ( next != eol_ ) {
	if ( *next == ' ' || *next == '\t' ) {
	  next++;
	} else {
	  break;
	}
      }
      if ( next != eol_ ) {
	char *nextBlank = nextBlankOr ( next );
	char save;

	if ( nextBlank ) {
	  save = *nextBlank;
	  *nextBlank = '\0';
	  strcpy ( columnName_, next );
	  *nextBlank = save;
	  if ( strstr ( nextBlank, "FREE" ) )
	    freeFormat_ = true;
	} else {
	  strcpy ( columnName_, next );
	}
      } else {
	strcpy ( columnName_, "no_name" );
      }
      break;
    } else if ( card_[0] != '*' && card_[0] != '#' ) {
      section_ = OSI_UNKNOWN_SECTION;
      break;
    }
  }
}

//  ~OSIMpsio.  Destructor
OSIMpsio::~OSIMpsio (  )
{
}

void
strcpyAndCompress ( char *to, const char *from )
{
  int n = strlen ( from );
  int i;
  int nto = 0;

  for ( i = 0; i < n; i++ ) {
    if ( from[i] != ' ' ) {
      to[nto++] = from[i];
    }
  }
  if ( !nto )
    to[nto++] = ' ';
  to[nto] = '\0';
}

//  nextField
OSISectionType
OSIMpsio::nextField (  )
{
  mpsType_ = OSI_BLANK_COLUMN;
  // find next non blank character
  char *next = position_;

  while ( next != eol_ ) {
    if ( *next == ' ' || *next == '\t' ) {
      next++;
    } else {
      break;
    }
  }
  bool gotCard;

  if ( next == eol_ ) {
    gotCard = false;
  } else {
    gotCard = true;
  }
  while ( !gotCard ) {
    // need new image
    char *getit = fgets ( card_, MAX_CARD_LENGTH, fp_ );

    if ( !getit ) {
      return OSI_EOF_SECTION;
    }
    // strip off newline
    cleanCard();
    if ( card_[0] == ' ' ) {
      // not a section or comment
      position_ = card_;
      eol_ = card_ + strlen ( card_ );
      // get mps type and column name
      // scan to first non blank
      next = card_;
      while ( next != eol_ ) {
	if ( *next == ' ' || *next == '\t' ) {
	  next++;
	} else {
	  break;
	}
      }
      if ( next != eol_ ) {
	char *nextBlank = nextBlankOr ( next );
	int nchar;

	if ( nextBlank ) {
	  nchar = nextBlank - next;
	} else {
	  nchar = -1;
	}
	mpsType_ = OSI_BLANK_COLUMN;
	// special coding if RHS or RANGE, not free format and blanks
	if ( ( section_ != OSI_RHS_SECTION && section_ != OSI_RANGE_SECTION )
	     || freeFormat_ || strncmp ( card_ + 4, "        ", 8 ) ) {
	  // if columns section only look for first field if MARKER
	  if ( section_ == OSI_COLUMN_SECTION
	       && !strstr ( next, "'MARKER'" ) ) nchar = -1;
	  if ( nchar == allowedLength[section_] ) {
	    //could be a type
	    int i;

	    for ( i = startType[section_]; i < endType[section_]; i++ ) {
	      if ( !strncmp ( next, mpsTypes[i], nchar ) ) {
		mpsType_ = ( OSIMpsType ) i;
		break;
	      }
	    }
	    if ( mpsType_ != OSI_BLANK_COLUMN ) {
	      //we know all we need so we can skip over
	      next = nextBlank;
	      while ( next != eol_ ) {
		if ( *next == ' ' || *next == '\t' ) {
		  next++;
		} else {
		  break;
		}
	      }
	      if ( next == eol_ ) {
		// error
		position_ = eol_;
		mpsType_ = OSI_UNKNOWN_MPS_TYPE;
	      } else {
		nextBlank = nextBlankOr ( next );
	      }
	    }
	  }
	  if ( mpsType_ != OSI_UNKNOWN_MPS_TYPE ) {
	    // special coding if BOUND, not free format and blanks
	    if ( section_ != OSI_BOUND_SECTION ||
		 freeFormat_ || strncmp ( card_ + 4, "        ", 8 ) ) {
	      char save = '?';

	      if ( !freeFormat_ && eightChar_ && next == card_ + 4 ) {
		if ( eol_ - next >= 8 ) {
		  if ( *( next + 8 ) != ' ' && *( next + 8 ) != '\0' ) {
		    eightChar_ = false;
		  } else {
		    nextBlank = next + 8;
		  }
		  save = *nextBlank;
		  *nextBlank = '\0';
		} else {
		  nextBlank = NULL;
		}
	      } else {
		if ( nextBlank ) {
		  save = *nextBlank;
		  *nextBlank = '\0';
		}
	      }
	      strcpyAndCompress ( columnName_, next );
	      if ( nextBlank ) {
		*nextBlank = save;
		// on to next
		next = nextBlank;
	      } else {
		next = eol_;
	      }
	    } else {
	      // blank bounds name
	      strcpy ( columnName_, "        " );
	    }
	    while ( next != eol_ ) {
	      if ( *next == ' ' || *next == '\t' ) {
		next++;
	      } else {
		break;
	      }
	    }
	    if ( next == eol_ ) {
	      // error unless row section
	      position_ = eol_;
	      value_ = -1.0e100;
	      if ( section_ != OSI_ROW_SECTION )
		mpsType_ = OSI_UNKNOWN_MPS_TYPE;
	    } else {
	      nextBlank = nextBlankOr ( next );
	    }
	    if ( section_ != OSI_ROW_SECTION ) {
	      char save = '?';

	      if ( !freeFormat_ && eightChar_ && next == card_ + 14 ) {
		if ( eol_ - next >= 8 ) {
		  if ( *( next + 8 ) != ' ' && *( next + 8 ) != '\0' ) {
		    eightChar_ = false;
		  } else {
		    nextBlank = next + 8;
		  }
		  save = *nextBlank;
		  *nextBlank = '\0';
		} else {
		  nextBlank = NULL;
		}
	      } else {
		if ( nextBlank ) {
		  save = *nextBlank;
		  *nextBlank = '\0';
		}
	      }
	      strcpyAndCompress ( rowName_, next );
	      if ( nextBlank ) {
		*nextBlank = save;
		// on to next
		next = nextBlank;
	      } else {
		next = eol_;
	      }
	      while ( next != eol_ ) {
		if ( *next == ' ' || *next == '\t' ) {
		  next++;
		} else {
		  break;
		}
	      }
	      // special coding for markers
	      if ( section_ == OSI_COLUMN_SECTION &&
		   !strncmp ( rowName_, "'MARKER'", 8 ) && next != eol_ ) {
		if ( !strncmp ( next, "'INTORG'", 8 ) ) {
		  mpsType_ = OSI_INTORG;
		} else if ( !strncmp ( next, "'INTEND'", 8 ) ) {
		  mpsType_ = OSI_INTEND;
		} else if ( !strncmp ( next, "'SOSORG'", 8 ) ) {
		  if ( mpsType_ == OSI_BLANK_COLUMN )
		    mpsType_ = OSI_S1_COLUMN;
		} else if ( !strncmp ( next, "'SOSEND'", 8 ) ) {
		  mpsType_ = OSI_SOSEND;
		} else {
		  mpsType_ = OSI_UNKNOWN_MPS_TYPE;
		}
		position_ = eol_;
		return section_;
	      }
	      if ( next == eol_ ) {
		// error unless bounds
		position_ = eol_;
		value_ = -1.0e100;
		if ( section_ != OSI_BOUND_SECTION )
		  mpsType_ = OSI_UNKNOWN_MPS_TYPE;
	      } else {
		nextBlank = nextBlankOr ( next );
		if ( nextBlank ) {
		  save = *nextBlank;
		  *nextBlank = '\0';
		}
		value_ = -1.0e100;
		sscanf ( next, "%lg", &value_ );
		//may have to allow some extra formats
		assert ( value_ > -1.0e50 );
		if ( nextBlank ) {
		  *nextBlank = save;
		  position_ = nextBlank;
		} else {
		  position_ = eol_;
		}
	      }
	    }
	  }
	} else {
	  //blank name in RHS or RANGE
	  strcpy ( columnName_, "        " );
	  char save = '?';

	  if ( !freeFormat_ && eightChar_ && next == card_ + 14 ) {
	    if ( eol_ - next >= 8 ) {
	      if ( *( next + 8 ) != ' ' && *( next + 8 ) != '\0' ) {
		eightChar_ = false;
	      } else {
		nextBlank = next + 8;
	      }
	      save = *nextBlank;
	      *nextBlank = '\0';
	    } else {
	      nextBlank = NULL;
	    }
	  } else {
	    if ( nextBlank ) {
	      save = *nextBlank;
	      *nextBlank = '\0';
	    }
	  }
	  strcpyAndCompress ( rowName_, next );
	  if ( nextBlank ) {
	    *nextBlank = save;
	    // on to next
	    next = nextBlank;
	  } else {
	    next = eol_;
	  }
	  while ( next != eol_ ) {
	    if ( *next == ' ' || *next == '\t' ) {
	      next++;
	    } else {
	      break;
	    }
	  }
	  if ( next == eol_ ) {
	    // error 
	    position_ = eol_;
	    value_ = -1.0e100;
	    mpsType_ = OSI_UNKNOWN_MPS_TYPE;
	  } else {
	    nextBlank = nextBlankOr ( next );
	    value_ = -1.0e100;
	    if ( nextBlank ) {
	      save = *nextBlank;
	      *nextBlank = '\0';
	    }
	    sscanf ( next, "%lg", &value_ );
	    //may have to allow some extra formats
	    assert ( value_ > -1.0e50 );
	    if ( nextBlank ) {
	      *nextBlank = save;
	      position_ = nextBlank;
	    } else {
	      position_ = eol_;
	    }
	  }
	}
      }
      return section_;
    } else if ( card_[0] != '*' ) {
      // not a comment
      int i;

      {
	std::cout<<"At line "<< cardNumber_ <<" "<< card_<<std::endl;
      }
      for ( i = OSI_ROW_SECTION; i < OSI_UNKNOWN_SECTION; i++ ) {
	if ( !strncmp ( card_, section[i], strlen ( section[i] ) ) ) {
	  break;
	}
      }
      position_ = card_;
      eol_ = card_;
      section_ = ( OSISectionType ) i;
      return section_;
    } else {
      // comment
    }
  }
  // we only get here for second field (we even allow more???)
  {
    char save = '?';
    char *nextBlank = nextBlankOr ( next );

    if ( !freeFormat_ && eightChar_ && next == card_ + 39 ) {
      if ( eol_ - next >= 8 ) {
	if ( *( next + 8 ) != ' ' && *( next + 8 ) != '\0' ) {
	  eightChar_ = false;
	} else {
	  nextBlank = next + 8;
	}
	save = *nextBlank;
	*nextBlank = '\0';
      } else {
	nextBlank = NULL;
      }
    } else {
      if ( nextBlank ) {
	save = *nextBlank;
	*nextBlank = '\0';
      }
    }
    strcpyAndCompress ( rowName_, next );
    // on to next
    if ( nextBlank ) {
      *nextBlank = save;
      next = nextBlank;
    } else {
      next = eol_;
    }
    while ( next != eol_ ) {
      if ( *next == ' ' || *next == '\t' ) {
	next++;
      } else {
	break;
      }
    }
    if ( next == eol_ ) {
      // error
      position_ = eol_;
      mpsType_ = OSI_UNKNOWN_MPS_TYPE;
    } else {
      nextBlank = nextBlankOr ( next );
    }
    if ( nextBlank ) {
      save = *nextBlank;
      *nextBlank = '\0';
    }
    value_ = -1.0e100;
    sscanf ( next, "%lg", &value_ );
    //may have to allow some extra formats
    assert ( value_ > -1.0e50 );
    if ( nextBlank ) {
      *nextBlank = save;
      position_ = nextBlank;
    } else {
      position_ = eol_;
    }
  }
  return section_;
}
static int
hash ( const char *name, int maxsiz, int length )
{
  static int mmult[] = {
    262139, 259459, 256889, 254291, 251701, 249133, 246709, 244247,
    241667, 239179, 236609, 233983, 231289, 228859, 226357, 223829,
    221281, 218849, 216319, 213721, 211093, 208673, 206263, 203773,
    201233, 198637, 196159, 193603, 191161, 188701, 186149, 183761,
    181303, 178873, 176389, 173897, 171469, 169049, 166471, 163871,
    161387, 158941, 156437, 153949, 151531, 149159, 146749, 144299,
    141709, 139369, 136889, 134591, 132169, 129641, 127343, 124853,
    122477, 120163, 117757, 115361, 112979, 110567, 108179, 105727,
    103387, 101021, 98639, 96179, 93911, 91583, 89317, 86939, 84521,
    82183, 79939, 77587, 75307, 72959, 70793, 68447, 66103
  };
  int n = 0;
  int j;

  for ( j = 0; j < length; ++j ) {
    int iname = name[j];

    n += mmult[j] * iname;
  }
  return ( abs ( n ) % maxsiz );	/* integer abs */
}

//  startHash.  Creates hash list for names
void
OsiMpsReader::startHash ( char **names, const OSIColumnIndex number , int section )
{
  names_[section] = names;
  numberHash_[section] = number;
  startHash(section);
}
void
OsiMpsReader::startHash ( int section ) const
{
  char ** names = names_[section];
  OSIColumnIndex number = numberHash_[section];
  OSIColumnIndex i;
  OSIColumnIndex maxhash = 4 * number;
  OSIColumnIndex ipos, iput;

  //hash_=(OsiHashLink *) malloc(maxhash*sizeof(OsiHashLink));
  hash_[section] = new OsiHashLink[maxhash];
  
  OsiHashLink * hashThis = hash_[section];

  for ( i = 0; i < maxhash; i++ ) {
    hashThis[i].index = -1;
    hashThis[i].next = -1;
  }

  /*
   * Initialize the hash table.  Only the index of the first name that
   * hashes to a value is entered in the table; subsequent names that
   * collide with it are not entered.
   */
  for ( i = 0; i < number; ++i ) {
    char *thisName = names[i];
    int length = strlen ( thisName );

    ipos = hash ( thisName, maxhash, length );
    if ( hashThis[ipos].index == -1 ) {
      hashThis[ipos].index = i;
    }
  }

  /*
   * Now take care of the names that collided in the preceding loop,
   * by finding some other entry in the table for them.
   * Since there are as many entries in the table as there are names,
   * there must be room for them.
   */
  iput = -1;
  for ( i = 0; i < number; ++i ) {
    char *thisName = names[i];
    int length = strlen ( thisName );

    ipos = hash ( thisName, maxhash, length );

    while ( 1 ) {
      OSIColumnIndex j1 = hashThis[ipos].index;

      if ( j1 == i )
	break;
      else {
	char *thisName2 = names[j1];

	if ( strcmp ( thisName, thisName2 ) == 0 ) {
	  printf ( "** duplicate name %s\n", names[i] );
	  break;
	} else {
	  OSIColumnIndex k = hashThis[ipos].next;

	  if ( k == -1 ) {
	    while ( 1 ) {
	      ++iput;
	      if ( iput > number ) {
		printf ( "** too many names\n" );
		break;
	      }
	      if ( hashThis[iput].index == -1 ) {
		break;
	      }
	    }
	    hashThis[ipos].next = iput;
	    hashThis[iput].index = i;
	    break;
	  } else {
	    ipos = k;
	    /* nothing worked - try it again */
	  }
	}
      }
    }
  }
}

//  stopHash.  Deletes hash storage
void
OsiMpsReader::stopHash ( int section )
{
  delete [] hash_[section];
  hash_[section] = NULL;
}

//  findHash.  -1 not found
OSIColumnIndex
OsiMpsReader::findHash ( const char *name , int section ) const
{
  OSIColumnIndex found = -1;

#if 0
  OSIColumnIndex i;

  for ( i = 0; i < numberHash_; i++ ) {
    if ( !strcmp ( names_[section][i], name ) ) {
      found = i;
      break;
    }
  }
#else
  char ** names = names_[section];
  OsiHashLink * hashThis = hash_[section];
  OSIColumnIndex maxhash = 4 * numberHash_[section];
  OSIColumnIndex ipos;

  /* default if we don't find anything */
  if ( !maxhash )
    return -1;
  int length = strlen ( name );

  ipos = hash ( name, maxhash, length );
  while ( 1 ) {
    OSIColumnIndex j1 = hashThis[ipos].index;

    if ( j1 >= 0 ) {
      char *thisName2 = names[j1];

      if ( strcmp ( name, thisName2 ) != 0 ) {
	OSIColumnIndex k = hashThis[ipos].next;

	if ( k != -1 )
	  ipos = k;
	else
	  break;
      } else {
	found = j1;
	break;
      }
    } else {
      found = -1;
      break;
    }
  }
#endif
  return found;
}

//------------------------------------------------------------------
// Get value for infinity
//------------------------------------------------------------------
double OsiMpsReader::getInfinity() const
{
  return infinity_;
}
//------------------------------------------------------------------
// Set value for infinity
//------------------------------------------------------------------
void OsiMpsReader::setInfinity(double value) 
{
  if ( value >= 1.020 ) {
    infinity_ = value;
  } else {
    std::cout << "Illegal value for infinity of " << value << endl;
  }

}
// Set file name
void OsiMpsReader::setFileName(const char * name)
{
  free(fileName_);
  fileName_=strdup(name);
}
// Get file name
const char * OsiMpsReader::getFileName() const
{
  return fileName_;
}
// Test if current file exists
const bool OsiMpsReader::fileExists() const
{
  // I am opening it to make sure not odd
  FILE *fp;
  if (strcmp(fileName_,"stdin")) {
    fp = fopen ( fileName_, "r" );
  } else {
    fp = stdin;
  }
  if (!fp) {
    return false;
  } else {
    fclose(fp);
    return true;
  }
}
#define MAX_INTEGER 1000000
// Sets default upper bound for integer variables
void OsiMpsReader::setDefaultBound(int value)
{
  if ( value >= 1 && value <=MAX_INTEGER ) {
    defaultBound_ = value;
  } else {
    std::cout << "Illegal default integer bound of " << value << endl;
  }
}
// gets default upper bound for integer variables
int OsiMpsReader::getDefaultBound() const
{
  return defaultBound_;
}
//------------------------------------------------------------------
// Read mps files
//------------------------------------------------------------------
int OsiMpsReader::readMps(const char * filename,  const char * extension)
{
  // clean up later as pointless using strings
  std::string f(filename);
  std::string fullname;
  free(fileName_);
  if (f!="stdin"&&f!="-") {
    std::string e(extension);
    if (f.find('.')<=0) 
      fullname = f + "." + e;
    else
      fullname = f;
    fileName_=strdup(fullname.c_str());    
  } else {
    fileName_=strdup("stdin");    
  }
  return readMps();
}
int OsiMpsReader::readMps()
{
  FILE *fp;
  if (strcmp(fileName_,"stdin")) {
    fp = fopen ( fileName_, "r" );
  } else {
    fp = stdin;
  }
  if (!fp) {
    std::cout << "Unable to open file " << fileName_ << endl;
    return -1;
  }
  OSIRowIndex numberRows;
  OSIColumnIndex numberColumns;
  OSIElementIndex numberElements;
  char name[100];
  bool ifmps;
  OSIMpsio mpsfile ( fp );

  if ( mpsfile.whichSection (  ) == OSI_NAME_SECTION ) {
    ifmps = true;
    strcpy ( name, mpsfile.columnName (  ) );
  } else if ( mpsfile.whichSection (  ) == OSI_UNKNOWN_SECTION ) {
    std::cout << "Unknown image " << mpsfile.
      card (  ) << " at line 1 of file" << fileName_ << endl;
    return -2;
  } else if ( mpsfile.whichSection (  ) != OSI_EOF_SECTION ) {
    strcpy ( name, mpsfile.card (  ) );
    ifmps = false;
  } else {
    std::cout << "EOF on file" << fileName_ << endl;
    return -3;
  }
  double *rowlower;
  double *rowupper;
  double *collower;
  double *colupper;
  double *objective;
  OSIElementIndex *start;
  OSIRowIndex *row;
  double *element;
  double objectiveOffset = 0.0;

  int numberErrors = 0;
  int i;

  if ( ifmps ) {
    // mps file - always read in free format
    bool gotNrow = false;

    //get ROWS
    assert ( mpsfile.nextField (  ) == OSI_ROW_SECTION );
    //use malloc etc as I don't know how to do realloc in C++
    numberRows = 0;
    numberColumns = 0;
    numberElements = 0;
    OSIRowIndex maxRows = 1000;
    char objectiveName[200];
    OSIMpsType *rowType =

      ( OSIMpsType * ) malloc ( maxRows * sizeof ( OSIMpsType ) );
    char **rowName = ( char ** ) malloc ( maxRows * sizeof ( char * ) );

    // for discarded free rows
    OSIRowIndex maxFreeRows = 100;
    OSIRowIndex numberOtherFreeRows = 0;
    char **freeRowName =

      ( char ** ) malloc ( maxFreeRows * sizeof ( char * ) );
    while ( mpsfile.nextField (  ) == OSI_ROW_SECTION ) {
      switch ( mpsfile.mpsType (  ) ) {
      case OSI_N_ROW:
	if ( !gotNrow ) {
	  gotNrow = true;
	  strcpy ( objectiveName, mpsfile.columnName (  ) );
	} else {
	  // add to discard list
	  if ( numberOtherFreeRows == maxFreeRows ) {
	    maxFreeRows = ( 3 * maxFreeRows ) / 2 + 100;
	    freeRowName =
	      ( char ** ) realloc ( freeRowName,

				    maxFreeRows * sizeof ( char * ) );
	  }
	  freeRowName[numberOtherFreeRows] =
	    strdup ( mpsfile.columnName (  ) );
	  numberOtherFreeRows++;
	}
	break;
      case OSI_E_ROW:
      case OSI_L_ROW:
      case OSI_G_ROW:
	if ( numberRows == maxRows ) {
	  maxRows = ( 3 * maxRows ) / 2 + 1000;
	  rowType =
	    ( OSIMpsType * ) realloc ( rowType,
				       maxRows * sizeof ( OSIMpsType ) );
	  rowName =

	    ( char ** ) realloc ( rowName, maxRows * sizeof ( char * ) );
	}
	rowType[numberRows] = mpsfile.mpsType (  );
	rowName[numberRows] = strdup ( mpsfile.columnName (  ) );
	numberRows++;
	break;
      default:
	numberErrors++;
	if ( numberErrors < 100 ) {
	  std::cout << "Bad image at card " << mpsfile.
	    cardNumber (  ) << " " << mpsfile.card (  ) << endl;
	} else if (numberErrors > 100000) {
	  std::cout << "Returning as too many errors"<< endl;
	  return numberErrors;
	}
      }
    }
    assert ( mpsfile.whichSection (  ) == OSI_COLUMN_SECTION );
    assert ( gotNrow );
    rowType =
      ( OSIMpsType * ) realloc ( rowType,
				 numberRows * sizeof ( OSIMpsType ) );
    // put objective and other free rows at end
    rowName =
      ( char ** ) realloc ( rowName,
			    ( numberRows + 1 +

			      numberOtherFreeRows ) * sizeof ( char * ) );
    rowName[numberRows] = objectiveName;
    memcpy ( rowName + numberRows + 1, freeRowName,
	     numberOtherFreeRows * sizeof ( char * ) );

    startHash ( rowName, numberRows + 1 + numberOtherFreeRows , 0 );
    OSIColumnIndex maxColumns = 1000 + numberRows / 5;
    OSIElementIndex maxElements = 5000 + numberRows / 2;
    OSIMpsType *columnType = ( OSIMpsType * )
      malloc ( maxColumns * sizeof ( OSIMpsType ) );
    char **columnName = ( char ** ) malloc ( maxColumns * sizeof ( char * ) );

    objective = ( double * ) malloc ( maxColumns * sizeof ( double ) );
    start = ( OSIElementIndex * )
      malloc ( ( maxColumns + 1 ) * sizeof ( OSIElementIndex ) );
    row = ( OSIRowIndex * )
      malloc ( maxElements * sizeof ( OSIRowIndex ) );
    element =
      ( double * ) malloc ( maxElements * sizeof ( double ) );
    // for duplicates
    OSIElementIndex *rowUsed = new OSIElementIndex[numberRows];

    for (i=0;i<numberRows;i++) {
      rowUsed[i]=-1;
    }
    bool objUsed = false;

    numberElements = 0;
    char lastColumn[200];

    memset ( lastColumn, '\0', 200 );
    OSIColumnIndex column = -1;
    bool inIntegerSet = false;
    OSIColumnIndex numberIntegers = 0;
    const double tinyElement = 1.0e-14;

    while ( mpsfile.nextField (  ) == OSI_COLUMN_SECTION ) {
      switch ( mpsfile.mpsType (  ) ) {
      case OSI_BLANK_COLUMN:
	if ( strcmp ( lastColumn, mpsfile.columnName (  ) ) ) {
	  // new column

	  // reset old column and take out tiny
	  if ( numberColumns ) {
	    objUsed = false;
	    OSIElementIndex i;
	    OSIElementIndex k = start[column];

	    for ( i = k; i < numberElements; i++ ) {
	      OSIRowIndex irow = row[i];

	      if ( fabs ( element[i] ) > tinyElement ) {
		element[k++] = element[i];
	      }
	      rowUsed[irow] = -1;
	    }
	    numberElements = k;
	  }
	  column = numberColumns;
	  if ( numberColumns == maxColumns ) {
	    maxColumns = ( 3 * maxColumns ) / 2 + 1000;
	    columnType = ( OSIMpsType * )
	      realloc ( columnType, maxColumns * sizeof ( OSIMpsType ) );
	    columnName = ( char ** )
	      realloc ( columnName, maxColumns * sizeof ( char * ) );

	    objective = ( double * )
	      realloc ( objective, maxColumns * sizeof ( double ) );
	    start = ( OSIElementIndex * )
	      realloc ( start,
			( maxColumns + 1 ) * sizeof ( OSIElementIndex ) );
	  }
	  if ( !inIntegerSet ) {
	    columnType[column] = OSI_UNSET_BOUND;
	  } else {
	    columnType[column] = OSI_INTORG;
	    numberIntegers++;
	  }
	  columnName[column] = strdup ( mpsfile.columnName (  ) );
	  strcpy ( lastColumn, mpsfile.columnName (  ) );
	  objective[column] = 0.0;
	  start[column] = numberElements;
	  numberColumns++;
	}
	if ( fabs ( mpsfile.value (  ) ) > tinyElement ) {
	  if ( numberElements == maxElements ) {
	    maxElements = ( 3 * maxElements ) / 2 + 1000;
	    row = ( OSIRowIndex * )
	      realloc ( row, maxElements * sizeof ( OSIRowIndex ) );
	    element = ( double * )
	      realloc ( element, maxElements * sizeof ( double ) );
	  }
	  // get row number
	  OSIRowIndex irow = findHash ( mpsfile.rowName (  ) , 0 );

	  if ( irow >= 0 ) {
	    double value = mpsfile.value (  );

	    // check for duplicates
	    if ( irow == numberRows ) {
	      // objective
	      if ( objUsed ) {
		numberErrors++;
		if ( numberErrors < 100 ) {
		  std::cout << "Duplicate objective at card " <<
		    mpsfile.cardNumber (  )
		    << " " << mpsfile.card (  ) << endl;
		} else if (numberErrors > 100000) {
		  std::cout << "Returning as too many errors"<< endl;
		  return numberErrors;
		}
	      } else {
		objUsed = true;
	      }
	      value += objective[column];
	      if ( fabs ( value ) <= tinyElement )
		value = 0.0;
	      objective[column] = value;
	    } else if ( irow < numberRows ) {
	      // other free rows will just be discarded so won't get here
	      if ( rowUsed[irow] >= 0 ) {
		element[rowUsed[irow]] += value;
		numberErrors++;
		if ( numberErrors < 100 ) {
		  std::cout << "Duplicate row at card " << mpsfile.
		    cardNumber (  ) << " " << mpsfile.
		    card (  ) << " " << mpsfile.rowName (  ) << endl;
		} else if (numberErrors > 100000) {
		  std::cout << "Returning as too many errors"<< endl;
		  return numberErrors;
		}
	      } else {
		row[numberElements] = irow;
		element[numberElements] = value;
		rowUsed[irow] = numberElements;
		numberElements++;
	      }
	    }
	  } else {
	    numberErrors++;
	    if ( numberErrors < 100 ) {
	      std::cout << "No match for row at card " << mpsfile.
		cardNumber (  ) << " " << mpsfile.card (  ) << " " << mpsfile.
		rowName (  ) << endl;
	    } else if (numberErrors > 100000) {
	      std::cout << "Returning as too many errors"<< endl;
	      return numberErrors;
	    }
	  }
	}
	break;
      case OSI_INTORG:
	inIntegerSet = true;
	break;
      case OSI_INTEND:
	inIntegerSet = false;
	break;
      case OSI_S1_COLUMN:
      case OSI_S2_COLUMN:
      case OSI_S3_COLUMN:
      case OSI_SOSEND:
	std::cout << "** code sos etc later" << endl;
	abort (  );
	break;
      default:
	numberErrors++;
	if ( numberErrors < 100 ) {
	  std::cout << "Bad image at card " << mpsfile.
	    cardNumber (  ) << " " << mpsfile.card (  ) << endl;
	} else if (numberErrors > 100000) {
	  std::cout << "Returning as too many errors"<< endl;
	  return numberErrors;
	}
      }
    }
    start[numberColumns] = numberElements;
    delete[]rowUsed;
    assert ( mpsfile.whichSection (  ) == OSI_RHS_SECTION );
    columnType =
      ( OSIMpsType * ) realloc ( columnType,
				 numberColumns * sizeof ( OSIMpsType ) );
    columnName =

      ( char ** ) realloc ( columnName, numberColumns * sizeof ( char * ) );
    objective = ( double * )
      realloc ( objective, numberColumns * sizeof ( double ) );
    start = ( OSIElementIndex * )
      realloc ( start, ( numberColumns + 1 ) * sizeof ( OSIElementIndex ) );
    row = ( OSIRowIndex * )
      realloc ( row, numberElements * sizeof ( OSIRowIndex ) );
    element = ( double * )
      realloc ( element, numberElements * sizeof ( double ) );

    rowlower = ( double * ) malloc ( numberRows * sizeof ( double ) );
    rowupper = ( double * ) malloc ( numberRows * sizeof ( double ) );
    for (i=0;i<numberRows;i++) {
      rowlower[i]=-infinity_;
      rowupper[i]=infinity_;
    }
    objUsed = false;
    memset ( lastColumn, '\0', 200 );
    bool gotRhs = false;

    // need coding for blank rhs
    while ( mpsfile.nextField (  ) == OSI_RHS_SECTION ) {
      OSIRowIndex irow;

      switch ( mpsfile.mpsType (  ) ) {
      case OSI_BLANK_COLUMN:
	if ( strcmp ( lastColumn, mpsfile.columnName (  ) ) ) {

	  // skip rest if got a rhs
	  if ( gotRhs ) {
	    while ( mpsfile.nextField (  ) == OSI_RHS_SECTION ) {
	    }
	  } else {
	    gotRhs = true;
	    strcpy ( lastColumn, mpsfile.columnName (  ) );
	  }
	}
	// get row number
	irow = findHash ( mpsfile.rowName (  ) , 0 );
	if ( irow >= 0 ) {
	  double value = mpsfile.value (  );

	  // check for duplicates
	  if ( irow == numberRows ) {
	    // objective
	    if ( objUsed ) {
	      numberErrors++;
	      if ( numberErrors < 100 ) {
		std::cout << "Duplicate objective at card " <<
		  mpsfile.cardNumber (  )
		  << " " << mpsfile.card (  ) << endl;
	      } else if (numberErrors > 100000) {
		std::cout << "Returning as too many errors"<< endl;
		return numberErrors;
	      }
	    } else {
	      objUsed = true;
	    }
	    objectiveOffset += value;
	  } else if ( irow < numberRows ) {
	    if ( rowlower[irow] != -infinity_ ) {
	      numberErrors++;
	      if ( numberErrors < 100 ) {
		std::cout << "Duplicate row at card " << mpsfile.
		  cardNumber (  ) << " " << mpsfile.
		  card (  ) << " " << mpsfile.rowName (  ) << endl;
	      } else if (numberErrors > 100000) {
		std::cout << "Returning as too many errors"<< endl;
		return numberErrors;
	      }
	    } else {
	      rowlower[irow] = value;
	    }
	  }
	} else {
	  numberErrors++;
	  if ( numberErrors < 100 ) {
	    std::cout << "No match for row at card " << mpsfile.
	      cardNumber (  ) << " " << mpsfile.card (  ) << " " << mpsfile.
	      rowName (  ) << endl;
	  } else if (numberErrors > 100000) {
	    std::cout << "Returning as too many errors"<< endl;
	    return numberErrors;
	  }
	}
	break;
      default:
	numberErrors++;
	if ( numberErrors < 100 ) {
	  std::cout << "Bad image at card " << mpsfile.
	    cardNumber (  ) << " " << mpsfile.card (  ) << endl;
	} else if (numberErrors > 100000) {
	  std::cout << "Returning as too many errors"<< endl;
	  return numberErrors;
	}
      }
    }
    if ( mpsfile.whichSection (  ) == OSI_RANGE_SECTION ) {
      memset ( lastColumn, '\0', 200 );
      bool gotRange = false;
      OSIRowIndex irow;

      // need coding for blank range
      while ( mpsfile.nextField (  ) == OSI_RANGE_SECTION ) {
	switch ( mpsfile.mpsType (  ) ) {
	case OSI_BLANK_COLUMN:
	  if ( strcmp ( lastColumn, mpsfile.columnName (  ) ) ) {

	    // skip rest if got a range
	    if ( gotRange ) {
	      while ( mpsfile.nextField (  ) == OSI_RANGE_SECTION ) {
	      }
	    } else {
	      gotRange = true;
	      strcpy ( lastColumn, mpsfile.columnName (  ) );
	    }
	  }
	  // get row number
	  irow = findHash ( mpsfile.rowName (  ) , 0 );
	  if ( irow >= 0 ) {
	    double value = mpsfile.value (  );

	    // check for duplicates
	    if ( irow == numberRows ) {
	      // objective
	      numberErrors++;
	      if ( numberErrors < 100 ) {
		std::cout << "Duplicate objective at card " <<
		  mpsfile.cardNumber (  )
		  << " " << mpsfile.card (  ) << endl;
	      } else if (numberErrors > 100000) {
		std::cout << "Returning as too many errors"<< endl;
		return numberErrors;
	      }
	    } else {
	      if ( rowupper[irow] != infinity_ ) {
		numberErrors++;
		if ( numberErrors < 100 ) {
		  std::cout << "Duplicate row at card " << mpsfile.
		    cardNumber (  ) << " " << mpsfile.
		    card (  ) << " " << mpsfile.rowName (  ) << endl;
		} else if (numberErrors > 100000) {
		  std::cout << "Returning as too many errors"<< endl;
		  return numberErrors;
		}
	      } else {
		rowupper[irow] = value;
	      }
	    }
	  } else {
	    numberErrors++;
	    if ( numberErrors < 100 ) {
	      std::cout << "No match for row at card " << mpsfile.
		cardNumber (  ) << " " << mpsfile.card (  ) << " " << mpsfile.
		rowName (  ) << endl;
	    } else if (numberErrors > 100000) {
	      std::cout << "Returning as too many errors"<< endl;
	      return numberErrors;
	    }
	  }
	  break;
	default:
	  numberErrors++;
	  if ( numberErrors < 100 ) {
	    std::cout << "Bad image at card " << mpsfile.
	      cardNumber (  ) << " " << mpsfile.card (  ) << endl;
	  } else if (numberErrors > 100000) {
	    std::cout << "Returning as too many errors"<< endl;
	    return numberErrors;
	  }
	}
      }
    }
    stopHash ( 0 );
    // massage ranges
    {
      OSIRowIndex irow;

      for ( irow = 0; irow < numberRows; irow++ ) {
	double lo = rowlower[irow];
	double up = rowupper[irow];
	double up2 = rowupper[irow];	//range

	switch ( rowType[irow] ) {
	case OSI_E_ROW:
	  if ( lo == -infinity_ )
	    lo = 0.0;
	  if ( up == infinity_ ) {
	    up = lo;
	  } else if ( up > 0.0 ) {
	    up += lo;
	  } else {
	    up = lo;
	    lo += up2;
	  }
	  break;
	case OSI_L_ROW:
	  if ( lo == -infinity_ ) {
	    up = 0.0;
	  } else {
	    up = lo;
	    lo = -infinity_;
	  }
	  if ( up2 != infinity_ ) {
	    lo = up - fabs ( up2 );
	  }
	  break;
	case OSI_G_ROW:
	  if ( lo == -infinity_ ) {
	    lo = 0.0;
	    up = infinity_;
	  } else {
	    up = infinity_;
	  }
	  if ( up2 != infinity_ ) {
	    up = lo + fabs ( up2 );
	  }
	  break;
	}
	rowlower[irow] = lo;
	rowupper[irow] = up;
      }
    }
    free ( rowType );
#if 0
    // get rid of row names !
    {
      OSIRowIndex irow;

      for ( irow = 0; irow < numberRows; irow++ ) {
	free ( rowName[irow] );
      }
      for ( irow = 0; irow < numberOtherFreeRows; irow++ ) {
	free ( freeRowName[irow] );
      }
      free ( rowName );
      free ( freeRowName );
    }
#endif
    // default bounds
    collower = ( double * ) malloc ( numberColumns * sizeof ( double ) );
    colupper = ( double * ) malloc ( numberColumns * sizeof ( double ) );
    for (i=0;i<numberColumns;i++) {
      collower[i]=0.0;
      colupper[i]=infinity_;
    }
    // set up integer region just in case
    integerType_ = new char[numberColumns];

    for ( column = 0; column < numberColumns; column++ ) {
      if ( columnType[column] == OSI_INTORG ) {
	columnType[column] = OSI_UNSET_BOUND;
	integerType_[column] = 1;
      } else {
	integerType_[column] = 0;
      }
    }
    if ( mpsfile.whichSection (  ) == OSI_BOUND_SECTION ) {
      memset ( lastColumn, '\0', 200 );
      bool gotBound = false;

      startHash ( columnName, numberColumns , 1 );
      while ( mpsfile.nextField (  ) == OSI_BOUND_SECTION ) {
	if ( strcmp ( lastColumn, mpsfile.columnName (  ) ) ) {

	  // skip rest if got a bound
	  if ( gotBound ) {
	    while ( mpsfile.nextField (  ) == OSI_BOUND_SECTION ) {
	    }
	  } else {
	    gotBound = true;;
	    strcpy ( lastColumn, mpsfile.columnName (  ) );
	  }
	}
	// get column number
	OSIColumnIndex icolumn = findHash ( mpsfile.rowName (  ) , 1 );

	if ( icolumn >= 0 ) {
	  double value = mpsfile.value (  );
	  bool ifError = false;

	  switch ( mpsfile.mpsType (  ) ) {
	  case OSI_UP_BOUND:
	    if ( value == -1.0e100 )
	      ifError = true;
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	      if ( value < 0.0 ) {
		collower[icolumn] = -infinity_;
	      }
	    } else if ( columnType[icolumn] == OSI_LO_BOUND ) {
	      if ( value < collower[icolumn] ) {
		ifError = true;
	      } else if ( value < collower[icolumn] + tinyElement ) {
		value = collower[icolumn];
	      }
	    } else if ( columnType[icolumn] == OSI_MI_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    colupper[icolumn] = value;
	    columnType[icolumn] = OSI_UP_BOUND;
	    break;
	  case OSI_LO_BOUND:
	    if ( value == -1.0e100 )
	      ifError = true;
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else if ( columnType[icolumn] == OSI_UP_BOUND ||
			columnType[icolumn] == OSI_UI_BOUND ) {
	      if ( value > colupper[icolumn] ) {
		ifError = true;
	      } else if ( value > colupper[icolumn] - tinyElement ) {
		value = colupper[icolumn];
	      }
	    } else {
	      ifError = true;
	    }
	    collower[icolumn] = value;
	    columnType[icolumn] = OSI_LO_BOUND;
	    break;
	  case OSI_FX_BOUND:
	    if ( value == -1.0e100 )
	      ifError = true;
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    collower[icolumn] = value;
	    colupper[icolumn] = value;
	    columnType[icolumn] = OSI_FX_BOUND;
	    break;
	  case OSI_FR_BOUND:
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    collower[icolumn] = -infinity_;
	    colupper[icolumn] = infinity_;
	    columnType[icolumn] = OSI_FR_BOUND;
	    break;
	  case OSI_MI_BOUND:
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	      colupper[icolumn] = 0.0;
	    } else if ( columnType[icolumn] == OSI_UP_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    collower[icolumn] = -infinity_;
	    columnType[icolumn] = OSI_MI_BOUND;
	    break;
	  case OSI_PL_BOUND:
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    columnType[icolumn] = OSI_PL_BOUND;
	    break;
	  case OSI_UI_BOUND:
	    if ( value == -1.0e100 )
	      ifError = true;
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else if ( columnType[icolumn] == OSI_LO_BOUND ) {
	      if ( value < collower[icolumn] ) {
		ifError = true;
	      } else if ( value < collower[icolumn] + tinyElement ) {
		value = collower[icolumn];
	      }
	    } else {
	      ifError = true;
	    }
	    colupper[icolumn] = value;
	    columnType[icolumn] = OSI_UI_BOUND;
	    if ( !integerType_[icolumn] ) {
	      numberIntegers++;
	      integerType_[icolumn] = 1;
	    }
	    break;
	  case OSI_BV_BOUND:
	    if ( columnType[icolumn] == OSI_UNSET_BOUND ) {
	    } else {
	      ifError = true;
	    }
	    collower[icolumn] = 0.0;
	    colupper[icolumn] = 1.0;
	    columnType[icolumn] = OSI_BV_BOUND;
	    if ( !integerType_[icolumn] ) {
	      numberIntegers++;
	      integerType_[icolumn] = 1;
	    }
	    break;
	  default:
	    ifError = true;
	    break;
	  }
	  if ( ifError ) {
	    numberErrors++;
	    if ( numberErrors < 100 ) {
	      std::cout << "Bad image at card " << mpsfile.
		cardNumber (  ) << " " << mpsfile.card (  ) << endl;
	    } else if (numberErrors > 100000) {
	      std::cout << "Returning as too many errors"<< endl;
	      return numberErrors;
	    }
	  }
	} else {
	  numberErrors++;
	  if ( numberErrors < 100 ) {
	    std::cout << "No match for column at card " << mpsfile.
	      cardNumber (  ) << " " << mpsfile.card (  ) << " " << mpsfile.
	      rowName (  ) << endl;
	  } else if (numberErrors > 100000) {
	    std::cout << "Returning as too many errors"<< endl;
	    return numberErrors;
	  }
	}
      }
      stopHash ( 1 );
    }
    free ( columnType );
#if 0
    // get rid of column names !
    {
      OSIColumnIndex icolumn;

      for ( icolumn = 0; icolumn < numberColumns; icolumn++ ) {
	free ( columnName[icolumn] );
      }
      free ( columnName );
    }
#endif
    // clean up integers
    if ( !numberIntegers ) {
      delete[]integerType_;
      integerType_ = NULL;
    } else {
      OSIColumnIndex icolumn;

      for ( icolumn = 0; icolumn < numberColumns; icolumn++ ) {
	if ( integerType_[icolumn] ) {
	  assert ( collower[icolumn] >= -MAX_INTEGER );
	  // if 0 infinity make 0-1 ???
	  if ( columnType[icolumn] == OSI_UNSET_BOUND ) 
	    colupper[icolumn] = defaultBound_;
	  if ( colupper[icolumn] > MAX_INTEGER ) 
	    colupper[icolumn] = MAX_INTEGER;
	}
      }
    }
    assert ( mpsfile.whichSection (  ) == OSI_ENDATA_SECTION );
  } else {
    fscanf ( fp, "%d %d %d\n", &numberRows, &numberColumns, &numberElements );
    OSIColumnIndex i;

    rowlower = ( double * ) malloc ( numberRows * sizeof ( double ) );
    rowupper = ( double * ) malloc ( numberRows * sizeof ( double ) );
    for ( i = 0; i < numberRows; i++ ) {
      int j;

      fscanf ( fp, "%d %lg %lg\n", &j, &rowlower[i], &rowupper[i] );
      assert ( i == j );
    }
    collower = ( double * ) malloc ( numberColumns * sizeof ( double ) );
    colupper = ( double * ) malloc ( numberColumns * sizeof ( double ) );
    objective= ( double * ) malloc ( numberColumns * sizeof ( double ) );
    start = ( OSIElementIndex *) malloc ((numberColumns + 1) *
					sizeof (OSIElementIndex) );
    row = ( OSIRowIndex * ) malloc (numberElements * sizeof (OSIRowIndex));
    element = ( double * ) malloc (numberElements * sizeof (double) );

    start[0] = 0;
    numberElements = 0;
    for ( i = 0; i < numberColumns; i++ ) {
      int j;
      int n;

      fscanf ( fp, "%d %d %lg %lg %lg\n", &j, &n, &collower[i], &colupper[i],
	       &objective[i] );
      assert ( i == j );
      for ( j = 0; j < n; j++ ) {
	fscanf ( fp, "       %d %lg\n", &row[numberElements],
		 &element[numberElements] );
	numberElements++;
      }
      start[i + 1] = numberElements;
    }
  }
  // construct packed matrix
  matrixByColumn_ = 
    new OsiPackedMatrix(true,
			numberRows,numberColumns,numberElements,
			element,row,start,NULL);
  free ( row );
  free ( start );
  free ( element );

  rowlower_ = rowlower;
  rowupper_ = rowupper;
  collower_ = collower;
  colupper_ = colupper;
  objective_ = objective;
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  numberElements_ = numberElements;
  std::cout<<"Problem has " << numberRows_ << " rows, " << numberColumns_
	   << " columns and " << numberElements_ << " elements" <<endl;
  fclose ( fp );
  return numberErrors;
}

//------------------------------------------------------------------
// Get number of rows, columns and elements
//------------------------------------------------------------------
int OsiMpsReader::getNumCols() const
{
  return numberColumns_;
}
int OsiMpsReader::getNumRows() const
{
  return numberRows_;
}
int OsiMpsReader::getNumElements() const
{
  return numberElements_;
}

//------------------------------------------------------------------
// Get pointer to column lower and upper bounds.
// Set element of column lower & upper bounds.
//------------------------------------------------------------------  
const double * OsiMpsReader::getColLower() const
{
  return collower_;
}
const double * OsiMpsReader::getColUpper() const
{
  return colupper_;
}

//------------------------------------------------------------------
// Get pointer to row lower and upper bounds.
//------------------------------------------------------------------  
const double * OsiMpsReader::getRowLower() const
{
  return rowlower_;
}
const double * OsiMpsReader::getRowUpper() const
{
  return rowupper_;
}
 
/** A quick inlined function to convert from lb/ub stryle constraint
    definition to sense/rhs/range style */
inline void
OsiMpsReader::convertBoundToSense(const double lower, const double upper,
					char& sense, double& right,
					double& range) const
{
  range = 0.0;
  if (lower > -infinity_) {
    if (upper < infinity_) {
      right = upper;
      if (upper==lower) {
        sense = 'E';
      } else {
        sense = 'R';
        range = upper - lower;
      }
    } else {
      sense = 'G';
      right = lower;
    }
  } else {
    if (upper < infinity_) {
      sense = 'L';
      right = upper;
    } else {
      sense = 'N';
      right = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
/** A quick inlined function to convert from sense/rhs/range stryle constraint
    definition to lb/ub style */
inline void
OsiMpsReader::convertSenseToBound(const char sense, const double right,
					const double range,
					double& lower, double& upper) const
{
  switch (sense) {
  case 'E':
    lower = upper = right;
    break;
  case 'L':
    lower = -infinity_;
    upper = right;
    break;
  case 'G':
    lower = right;
    upper = infinity_;
    break;
  case 'R':
    lower = right - range;
    upper = right;
    break;
  case 'N':
    lower = -infinity_;
    upper = infinity_;
    break;
  }
}
//------------------------------------------------------------------
// Get sense of row constraints.
//------------------------------------------------------------------ 
const char * OsiMpsReader::getRowSense() const
{
  if ( rowsense_==NULL ) {

    int nr=numberRows_;
    rowsense_ = (char *) malloc(nr*sizeof(char));


    double dum1,dum2;
    int i;
    for ( i=0; i<nr; i++ ) {
      convertBoundToSense(rowlower_[i],rowupper_[i],rowsense_[i],dum1,dum2);
    }
  }
  return rowsense_;
}

//------------------------------------------------------------------
// Get the rhs of rows.
//------------------------------------------------------------------ 
const double * OsiMpsReader::getRightHandSide() const
{
  if ( rhs_==NULL ) {

    int nr=numberRows_;
    rhs_ = (double *) malloc(nr*sizeof(double));


    char dum1;
    double dum2;
    int i;
    for ( i=0; i<nr; i++ ) {
      convertBoundToSense(rowlower_[i],rowupper_[i],dum1,rhs_[i],dum2);
    }
  }
  return rhs_;
}

//------------------------------------------------------------------
// Get the range of rows.
// Length of returned vector is getNumRows();
//------------------------------------------------------------------ 
const double * OsiMpsReader::getRowRange() const
{
  if ( rowrange_==NULL ) {

    int nr=numberRows_;
    rowrange_ = (double *) malloc(nr*sizeof(double));
    std::fill(rowrange_,rowrange_+nr,0.0);

    char dum1;
    double dum2;
    int i;
    for ( i=0; i<nr; i++ ) {
      convertBoundToSense(rowlower_[i],rowupper_[i],dum1,dum2,rowrange_[i]);
    }
  }
  return rowrange_;
}

const double * OsiMpsReader::getObjCoefficients() const
{
  return objective_;
}
 
//------------------------------------------------------------------
// Create a row copy of the matrix ...
//------------------------------------------------------------------
const OsiPackedMatrix * OsiMpsReader::getMatrixByRow() const
{
  if ( matrixByRow_ == NULL && matrixByColumn_) {
    matrixByRow_ = new OsiPackedMatrix(*matrixByColumn_);
    matrixByRow_->reverseOrdering();
  }
  return matrixByRow_;
}

//------------------------------------------------------------------
// Create a column copy of the matrix ...
//------------------------------------------------------------------
const OsiPackedMatrix * OsiMpsReader::getMatrixByCol() const
{
  return matrixByColumn_;
}



//------------------------------------------------------------------
// Return true if column is a continuous, binary, ...
//------------------------------------------------------------------
bool OsiMpsReader::isContinuous(int columnNumber) const
{
  const char * intType = integerType_;
  if ( intType==NULL ) return true;
  assert (columnNumber>=0 && columnNumber < numberColumns_);
  if ( intType[columnNumber]==0 ) return true;
  return false;
}

/* Return true if column is integer.
   Note: This function returns true if the the column
   is binary or a general integer.
*/
bool OsiMpsReader::isInteger(int columnNumber) const
{
  const char * intType = integerType_;
  if ( intType==NULL ) return false;
  assert (columnNumber>=0 && columnNumber < numberColumns_);
  if ( intType[columnNumber]==0 ) return true;
  return false;
}
// if integer
const char * OsiMpsReader::integer() const
{
  return integerType_;
}
// names - returns NULL if out of range
const char * OsiMpsReader::rowName(int index) const
{
  if (index>=0&&index<numberRows_) {
    return names_[0][index];
  } else {
    return NULL;
  }
}
const char * OsiMpsReader::columnName(int index) const
{
  if (index>=0&&index<numberColumns_) {
    return names_[1][index];
  } else {
    return NULL;
  }
}
// names - returns -1 if name not found
int OsiMpsReader::rowIndex(const char * name) const
{
  if (!hash_[0]) {
    if (numberRows_) {
      startHash(0);
    } else {
      return -1;
    }
  }
  return findHash ( name , 0 );
}
    int OsiMpsReader::columnIndex(const char * name) const
{
  if (!hash_[1]) {
    if (numberColumns_) {
      startHash(1);
    } else {
      return -1;
    }
  }
  return findHash ( name , 1 );
}

// Release all row information (lower, upper)
void OsiMpsReader::releaseRowInformation()
{
  free(rowlower_);
  free(rowupper_);
  rowlower_=NULL;
  rowupper_=NULL;
}
// Release all column information (lower, upper, objective)
void OsiMpsReader::releaseColumnInformation()
{
  free(collower_);
  free(colupper_);
  free(objective_);
  collower_=NULL;
  colupper_=NULL;
  objective_=NULL;
}
// Release integer information
void OsiMpsReader::releaseIntegerInformation()
{
  free(integerType_);
  integerType_=NULL;
}
// Release row names
void OsiMpsReader::releaseRowNames()
{
  releaseRedundantInformation();
  int i;
  for (i=0;i<numberHash_[0];i++) {
    free(names_[0][i]);
  }
  free(names_[0]);
  names_[0]=NULL;
  numberHash_[0]=0;
}
// Release column names
void OsiMpsReader::releaseColumnNames()
{
  releaseRedundantInformation();
  int i;
  for (i=0;i<numberHash_[1];i++) {
    free(names_[1][i]);
  }
  free(names_[1]);
  names_[1]=NULL;
  numberHash_[1]=0;
}
// Release matrix information
void OsiMpsReader::releaseMatrixInformation()
{
  releaseRedundantInformation();
  delete matrixByColumn_;
  matrixByColumn_=NULL;
}
  


//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
OsiMpsReader::OsiMpsReader ()
:
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
matrixByRow_(NULL),
matrixByColumn_(NULL),
rowlower_(NULL),
rowupper_(NULL),
collower_(NULL),
colupper_(NULL),
objective_(NULL),
integerType_(NULL),
fileName_(strdup("stdin")),
numberRows_(0),
numberColumns_(0),
numberElements_(0),
defaultBound_(1),
infinity_(DBL_MAX)
{
  numberHash_[0]=0;
  hash_[0]=NULL;
  names_[0]=NULL;
  numberHash_[1]=0;
  hash_[1]=NULL;
  names_[1]=NULL;
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
OsiMpsReader::OsiMpsReader (
                  const OsiMpsReader & source)
:
rowsense_(NULL),
rhs_(NULL),
rowrange_(NULL),
matrixByRow_(NULL),
matrixByColumn_(NULL),
rowlower_(NULL),
rowupper_(NULL),
collower_(NULL),
colupper_(NULL),
objective_(NULL),
integerType_(NULL),
fileName_(strdup("stdin")),
numberRows_(0),
numberColumns_(0),
numberElements_(0),
defaultBound_(1),
infinity_(DBL_MAX)
{
  numberHash_[0]=0;
  hash_[0]=NULL;
  names_[0]=NULL;
  numberHash_[1]=0;
  hash_[1]=NULL;
  names_[1]=NULL;
  if ( source.rowlower_ !=NULL || source.collower_ != NULL) {
    gutsOfCopy(source);
    // OK and proper to leave rowsense_, rhs_, and
    // rowrange_ (also row copy and hash) to NULL.  They will be constructed
    // if they are required.
  }
}

void OsiMpsReader::gutsOfCopy(const OsiMpsReader & rhs)
{
  if (rhs.matrixByColumn_)
    matrixByColumn_=new OsiPackedMatrix(*(rhs.matrixByColumn_));
  numberElements_=rhs.numberElements_;
  numberRows_=rhs.numberRows_;
  numberColumns_=rhs.numberColumns_;
  if (rhs.rowlower_) {
    rowlower_ = (double *) malloc(numberRows_*sizeof(double));
    rowupper_ = (double *) malloc(numberRows_*sizeof(double));
    memcpy(rowlower_,rhs.rowlower_,numberRows_*sizeof(double));
    memcpy(rowupper_,rhs.rowupper_,numberRows_*sizeof(double));
  }
  if (rhs.collower_) {
    collower_ = (double *) malloc(numberColumns_*sizeof(double));
    colupper_ = (double *) malloc(numberColumns_*sizeof(double));
    objective_ = (double *) malloc(numberColumns_*sizeof(double));
    memcpy(collower_,rhs.collower_,numberColumns_*sizeof(double));
    memcpy(colupper_,rhs.colupper_,numberColumns_*sizeof(double));
    memcpy(objective_,rhs.objective_,numberColumns_*sizeof(double));
  }
  if (rhs.integerType_) {
    integerType_ = (char *) malloc (numberColumns_*sizeof(char));
    memcpy(integerType_,rhs.integerType_,numberColumns_*sizeof(char));
  }
  fileName_ = strdup(rhs.fileName_);
  numberHash_[0]=rhs.numberHash_[0];
  numberHash_[1]=rhs.numberHash_[1];
  defaultBound_=rhs.defaultBound_;
  infinity_=rhs.infinity_;
  int section;
  for (section=0;section<2;section++) {
    if (numberHash_[section]) {
      char ** names2 = rhs.names_[section];
      names_[section] = (char **) malloc(numberHash_[section]*
					 sizeof(char *));
      char ** names = names_[section];
      int i;
      for (i=0;i<numberHash_[section];i++) {
	names[i]=strdup(names2[i]);
      }
    }
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
OsiMpsReader::~OsiMpsReader ()
{
  gutsOfDestructor();
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
OsiMpsReader &
OsiMpsReader::operator=(
                   const OsiMpsReader& rhs)
{
  if (this != &rhs) {    
    gutsOfDestructor();
    if ( rhs.rowlower_ !=NULL || rhs.collower_ != NULL) {
      gutsOfCopy(rhs);
    }
  }
  return *this;
}



//-------------------------------------------------------------------
void OsiMpsReader::gutsOfDestructor()
{  
  freeAll();
}


void OsiMpsReader::freeAll()
{  
  releaseRedundantInformation();
  delete matrixByRow_;
  delete matrixByColumn_;
  matrixByRow_=NULL;
  matrixByColumn_=NULL;
  free(rowlower_);
  free(rowupper_);
  free(collower_);
  free(colupper_);
  free(objective_);
  free(integerType_);
  free(names_[0]);
  free(names_[1]);
  rowlower_=NULL;
  rowupper_=NULL;
  collower_=NULL;
  colupper_=NULL;
  objective_=NULL;
  integerType_=NULL;
  names_[0]=NULL;
  names_[1]=NULL;
}

/* Release all information which can be re-calculated e.g. rowsense
    also any row copies OR hash tables for names */
void OsiMpsReader::releaseRedundantInformation()
{  
  free( rowsense_);
  free( rhs_);
  free( rowrange_);
  rowsense_=NULL;
  rhs_=NULL;
  rowrange_=NULL;
  free (hash_[0]);
  free (hash_[1]);
  hash_[0]=0;
  hash_[1]=0;
  delete matrixByRow_;
  matrixByRow_=NULL;
}

//#############################################################################


