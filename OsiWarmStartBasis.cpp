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
  // Round all so arrays multiple of 4
  int nint = (ns+15) >> 4;
  structuralStatus_ = new char[4*nint];
  // CoinFillN was used here - I am sure memset will be cleaner for char
  memset (structuralStatus_, 0, (4*nint) * sizeof(char));
  nint = (na+15) >> 4;
  artificialStatus_ = new char[4*nint];
  memset (artificialStatus_, 0, (4*nint) * sizeof(char));
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
  // Round all so arrays multiple of 4
  int nint = (ns+15) >> 4;
  structuralStatus_ = new char[4*nint];
  structuralStatus_[4*nint-3]=0;
  structuralStatus_[4*nint-2]=0;
  structuralStatus_[4*nint-1]=0;
  memcpy (structuralStatus_, sStat, ((ns + 3) / 4) * sizeof(char));
  nint = (na+15) >> 4;
  artificialStatus_ = new char[4*nint];
  artificialStatus_[4*nint-3]=0;
  artificialStatus_[4*nint-2]=0;
  artificialStatus_[4*nint-1]=0;
  memcpy (artificialStatus_, aStat, ((na + 3) / 4) * sizeof(char));
}

OsiWarmStartBasis::OsiWarmStartBasis(const OsiWarmStartBasis& ws) :
  numStructural_(ws.numStructural_), numArtificial_(ws.numArtificial_),
  structuralStatus_(NULL), artificialStatus_(NULL) {
  // Round all so arrays multiple of 4
  int nint = (numStructural_+15) >> 4;
  structuralStatus_ = new char[4*nint];
  memcpy (structuralStatus_, ws.structuralStatus_, 
	  (4*nint) * sizeof(char));
  nint = (numArtificial_+15) >> 4;
  artificialStatus_ = new char[4*nint];
  memcpy (artificialStatus_, ws.artificialStatus_, 
	  (4*nint) * sizeof(char));
}

OsiWarmStartBasis& 
OsiWarmStartBasis::operator=(const OsiWarmStartBasis& rhs)
{
  if (this != &rhs) {
    numStructural_=rhs.numStructural_;
    numArtificial_=rhs.numArtificial_;
    delete [] structuralStatus_;
    delete [] artificialStatus_;
    // Round all so arrays multiple of 4
    int nint = (numStructural_+15) >> 4;
    structuralStatus_ = new char[4*nint];
    memcpy (structuralStatus_, rhs.structuralStatus_, 
	    (4*nint) * sizeof(char));
    nint = (numArtificial_+15) >> 4;
    artificialStatus_ = new char[4*nint];
    memcpy (artificialStatus_, rhs.artificialStatus_, 
	  (4*nint) * sizeof(char));
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
    nCharOld  = 4*((numArtificial_+15)>>4);
    nCharNew  = 4*((newNumberRows+15)>>4);
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
    nCharOld  = 4*((numStructural_+15)>>4);
    nCharNew  = 4*((newNumberColumns+15)>>4);
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
  int nCharNew  = 4*((numArtificial_-numberDeleted+15)>>4);
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
  int nCharNew  = 4*((numStructural_-numberDeleted+3)>>4);
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
// Prints in readable format (for debug)
void 
OsiWarmStartBasis::print() const
{
  std::cout<<"Basis "<<this<<" has "<<numArtificial_<<" rows and "
	   <<numStructural_<<" columns"<<std::endl;
  std::cout<<"Rows:"<<std::endl;
  int i;
  char type[]={'F','B','U','L'};

  for (i=0;i<numArtificial_;i++) 
    std::cout<<type[getArtifStatus(i)];
  std::cout<<std::endl;
  std::cout<<"Columns:"<<std::endl;

  for (i=0;i<numStructural_;i++) 
    std::cout<<type[getStructStatus(i)];
  std::cout<<std::endl;
}
