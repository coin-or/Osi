// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiWarmStartBasis.hpp"
#include <cmath>
#include <iostream>

//#############################################################################

void 
OsiWarmStartBasis::setSize(int ns, int na) {
  delete[] structuralStatus_;
  delete[] artificialStatus_;
  structuralStatus_ = new char[(ns + 3) / 4];
  artificialStatus_ = new char[(na + 3) / 4];
  // CoinFillN was used here - I am sure memset will be cleaner for char
  memset (structuralStatus_, 0, ((ns + 3) / 4) * sizeof(char));
  memset (artificialStatus_, 0, ((na + 3) / 4) * sizeof(char));
  numArtificial_ = na;
  numStructural_ = ns;
}

void 
OsiWarmStartBasis::assignBasisStatus(int ns, int na, char*& sStat, 
				     char*& aStat) {
  delete[] structuralStatus_;
  delete[] artificialStatus_;
  numStructural_ = ns;
  numArtificial_ = na;
  structuralStatus_ = sStat;
  artificialStatus_ = aStat;
  sStat = 0;
  aStat = 0;
}
OsiWarmStartBasis::OsiWarmStartBasis(int ns, int na, 
				     const char* sStat, const char* aStat) :
  numStructural_(ns), numArtificial_(na),
  structuralStatus_(NULL), artificialStatus_(NULL) {
  structuralStatus_ = new char[(ns + 3) / 4];
  artificialStatus_ = new char[(na + 3) / 4];
  memcpy (structuralStatus_, sStat, ((ns + 3) / 4) * sizeof(char));
  memcpy (artificialStatus_, aStat, ((na + 3) / 4) * sizeof(char));
}

OsiWarmStartBasis::OsiWarmStartBasis(const OsiWarmStartBasis& ws) :
  numStructural_(ws.numStructural_), numArtificial_(ws.numArtificial_),
  structuralStatus_(NULL), artificialStatus_(NULL) {
  structuralStatus_ = new char[(numStructural_ + 3) / 4];
  artificialStatus_ = new char[(numArtificial_ + 3) / 4];
  memcpy (structuralStatus_, ws.structuralStatus_, 
	  ((numStructural_ + 3) / 4) * sizeof(char));
  memcpy (artificialStatus_, ws.artificialStatus_, 
	  ((numArtificial_ + 3) / 4) * sizeof(char));
}

OsiWarmStartBasis& 
OsiWarmStartBasis::operator=(const OsiWarmStartBasis& rhs)
{
  if (this != &rhs) {
    numStructural_=rhs.numStructural_;
    numArtificial_=rhs.numArtificial_;
    delete [] structuralStatus_;
    delete [] artificialStatus_;
    structuralStatus_ = new char[(numStructural_ + 3) / 4];
    artificialStatus_ = new char[(numArtificial_ + 3) / 4];
    memcpy (structuralStatus_, rhs.structuralStatus_, 
	    ((numStructural_ + 3) / 4) * sizeof(char));
    memcpy (artificialStatus_, rhs.artificialStatus_, 
	    ((numArtificial_ + 3) / 4) * sizeof(char));
  }
  return *this;
}

// Resizes 
void 
OsiWarmStartBasis::resize (int newNumberRows, int newNumberColumns)
{
  char * array;
  int i , nCharNew, nCharOld;
  if (newNumberRows!=numArtificial_) {
    nCharOld  = (numArtificial_+3)>>2;
    nCharNew  = (newNumberRows+3)>>2;
    array = new char[nCharNew];
    // zap all for clarity and zerofault etc
    memset(array,0,nCharNew*sizeof(char));
    memcpy(array,artificialStatus_,(nCharOld>nCharNew)?nCharNew:nCharOld);
    delete [] artificialStatus_;
    artificialStatus_ = array;
    for (i=numArtificial_;i<newNumberRows;i++) 
      setArtifStatus(i, basic);
    numArtificial_ = newNumberRows;
  }
  if (newNumberColumns!=numStructural_) {
    nCharOld  = (numStructural_+3)>>2;
    nCharNew  = (newNumberColumns+3)>>2;
    array = new char[nCharNew];
    // zap all for clarity and zerofault etc
    memset(array,0,nCharNew*sizeof(char));
    memcpy(array,structuralStatus_,(nCharOld>nCharNew)?nCharNew:nCharOld);
    delete [] structuralStatus_;
    structuralStatus_ = array;
    for (i=numStructural_;i<newNumberColumns;i++) 
      setStructStatus(i, atLowerBound);
    numStructural_ = newNumberColumns;
  }
}
// Deletes rows
void 
OsiWarmStartBasis::deleteRows(int number, const int * which)
{
  int i ;
  char * deleted = new char[numArtificial_];
  int numberDeleted=0;
  memset(deleted,0,numArtificial_*sizeof(char));
  for (i=0;i<number;i++) {
    int j = which[i];
    if (j>=0&&j<numArtificial_&&!deleted[j]) {
      numberDeleted++;
      deleted[j]=1;
    }
  }
  int nCharNew  = (numArtificial_-numberDeleted+3)>>2;
  char * array = new char[nCharNew];
  int put=0;
  int numberNotBasic=0;
  for (i=0;i<numArtificial_;i++) {
    Status status = getArtifStatus(i);
    if (!deleted[i]) {
      setArtifStatus(put, status);
      put++;
    } else if (status!=OsiWarmStartBasis::basic) {
      numberNotBasic++;
    }
  }
  memcpy(array,artificialStatus_,nCharNew*sizeof(char));
  delete [] artificialStatus_;
  artificialStatus_ = array;
  delete [] deleted;
  numArtificial_ -= numberDeleted;
#ifdef OSI_DEBUG
  if (numberNotBasic)
    std::cout<<numberNotBasic<<" non basic artificials deleted"<<std::endl;
#endif
}
// Deletes columns
void 
OsiWarmStartBasis::deleteColumns(int number, const int * which)
{
  int i ;
  char * deleted = new char[numStructural_];
  int numberDeleted=0;
  memset(deleted,0,numStructural_*sizeof(char));
  for (i=0;i<number;i++) {
    int j = which[i];
    if (j>=0&&j<numStructural_&&!deleted[j]) {
      numberDeleted++;
      deleted[j]=1;
    }
  }
  int nCharNew  = (numStructural_-numberDeleted+3)>>2;
  char * array = new char[nCharNew];
  int put=0;
  int numberBasic=0;
  for (i=0;i<numStructural_;i++) {
    Status status = getStructStatus(i);
    if (!deleted[i]) {
      setStructStatus(put, status);
      put++;
    } else if (status==OsiWarmStartBasis::basic) {
      numberBasic++;
    }
  }
  memcpy(array,structuralStatus_,nCharNew*sizeof(char));
  delete [] structuralStatus_;
  structuralStatus_ = array;
  delete [] deleted;
  numStructural_ -= numberDeleted;
#ifdef OSI_DEBUG
  if (numberBasic)
    std::cout<<numberBasic<<" basic structurals deleted"<<std::endl;
#endif
}
