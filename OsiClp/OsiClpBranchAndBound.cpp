// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>

#include "OsiClpSolverInterface.hpp"
#include "ClpSimplex.hpp"
#include "OsiOsiMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiCuts.hpp"
#ifdef CUTS
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#endif
#undef OSI_DEBUG
// This is just like OsiRowCut but has a count of number nodes using it
class OsiCountRowCut : public OsiRowCut {

public:
  
  /**@name Constructors and destructors */
  //@{
    /// Copy constructor 
    OsiCountRowCut ( const OsiCountRowCut &);  
    /// Copy constructor from OsiRowCut
    OsiCountRowCut ( const OsiRowCut &);  

    /// Default Constructor 
    OsiCountRowCut ();
  
    /// Destructor 
    virtual ~OsiCountRowCut ();
  //@}
  /// Increment number of references
  void increment()
  {numberPointingToThis_++;};

  /// Decrement number of references and return number left
  int decrement()
  {numberPointingToThis_--;return numberPointingToThis_;};

  /// Number of other OsiNodeInfos pointing to this
  int numberPointingToThis_;
};
// Default Constructor 
OsiCountRowCut::OsiCountRowCut ()
  :
  OsiRowCut(),
  numberPointingToThis_(0)
{
}

// Copy Constructor 
OsiCountRowCut::OsiCountRowCut (const OsiCountRowCut & rhs)
  : OsiRowCut(rhs),
  numberPointingToThis_(0)
{
}
   
// Copy Constructor 
OsiCountRowCut::OsiCountRowCut (const OsiRowCut & rhs)
  : OsiRowCut(rhs),
  numberPointingToThis_(0)
{
}
   
OsiCountRowCut::~OsiCountRowCut()
{
}


/* Node information base class. 
   This may be a complete description of status of model, or it could
   just be a difference between this and parent. */

class OsiNodeInfo {
public:
  // Default Constructor 
  OsiNodeInfo ();

  // Constructor given parent
  OsiNodeInfo (OsiNodeInfo * parent);
   
  virtual ~OsiNodeInfo();
  /// returns NULL or previous one to apply (if mode ==0 )
  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model, int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const = 0;
  /// Clone
  virtual OsiNodeInfo * clone() const = 0;

  /// Increment number of references
  void increment()
  {numberPointingToThis_++;};

  /// Decrement number of references and return number left
  int decrement()
  {numberPointingToThis_--;return numberPointingToThis_;};

  /// Number of other OsiNodeInfos pointing to this
  int numberPointingToThis_;

  /// Parent of this
  OsiNodeInfo * parent() const
  {return parent_;};
  OsiNodeInfo * parent_;

  // Add cuts (increments counts)
  void addCuts(int numberCuts,OsiCountRowCut ** cuts);
  void addCuts(OsiCuts & cuts);
  /** Delete cuts (decrements counts)
      Slow unless cuts in same order as saved
  */
  void deleteCuts(int numberToDelete,OsiCountRowCut ** cuts);
  void deleteCuts(int numberToDelete,int * which);

  /// Apply cuts
  void applyCuts(OsiSolverInterface * model) const;

  /// Number of row cuts (this node)
  int numberCuts() const
  {return numberCuts_;};
  int numberCuts_;

  /// Array of pointers to cuts
  OsiCountRowCut ** cuts() const
  {return cuts_;};
  OsiCountRowCut ** cuts_;
  
};

// Default Constructor 
OsiNodeInfo::OsiNodeInfo ()
  :
  numberPointingToThis_(0),
  parent_(NULL),
  numberCuts_(0),
  cuts_(NULL)
{
}

// Constructor given parent
OsiNodeInfo::OsiNodeInfo (OsiNodeInfo * parent)
  :
  numberPointingToThis_(0),
  parent_(parent),
  numberCuts_(0),
  cuts_(NULL)
{
  if (parent_) 
    parent_->increment();
}

OsiNodeInfo::~OsiNodeInfo()
{
  assert(!numberPointingToThis_);
  int i;
  for (i=0;i<numberCuts_;i++) {
    int number = cuts_[i]->decrement();
    if (!number)
      delete cuts_[i];
  }
  delete [] cuts_;
  if (parent_) {
    int numberLinks = parent_->decrement();
    if (!numberLinks)
      delete parent_;
  }
}

// Add cuts
void
OsiNodeInfo::addCuts(int numberCuts, OsiCountRowCut ** cuts)
{
  if (numberCuts) {
    int i;
    if (!numberCuts_) {
      cuts_ = new OsiCountRowCut * [numberCuts];
    } else {
      OsiCountRowCut ** temp = new OsiCountRowCut * [numberCuts+numberCuts_];
      memcpy(temp,cuts_,numberCuts_*sizeof(OsiCountRowCut *));
      delete [] cuts_;
      cuts_ = temp;
    }
    for (i=0;i<numberCuts;i++) {
      cuts[i]->increment();
      cuts_[numberCuts_++] = cuts[i];
    }
  }
}
// Add cuts
void
OsiNodeInfo::addCuts(OsiCuts & cuts)
{
  int numberCuts = cuts.sizeRowCuts();
  if (numberCuts) {
    int i;
    if (!numberCuts_) {
      cuts_ = new OsiCountRowCut * [numberCuts];
    } else {
      OsiCountRowCut ** temp = new OsiCountRowCut * [numberCuts+numberCuts_];
      memcpy(temp,cuts_,numberCuts_*sizeof(OsiCountRowCut *));
    }
    for (i=0;i<numberCuts;i++) {
      OsiCountRowCut * thisCut = new OsiCountRowCut(*cuts.rowCutPtr(i));
      thisCut->increment();
      cuts_[numberCuts_++] = thisCut;
    }
  }
}

// delete cuts
void
OsiNodeInfo::deleteCuts(int numberToDelete, OsiCountRowCut ** cuts)
{
  int i;
  int j;
  int last=-1;
  for (i=0;i<numberToDelete;i++) {
    OsiCountRowCut * next = cuts[i];
    for (j=last+1;j<numberCuts_;j++) {
      if (next==cuts_[j])
	break;
    }
    if (j==numberCuts_) {
      // start from beginning
      for (j=0;j<last;j++) {
	if (next==cuts_[j])
	  break;
      }
      assert(j<last);
    }
    last=j;
    int number = cuts_[j]->decrement();
    if (!number)
      delete cuts_[j];
    cuts_[j]=NULL;
  }
  j=0;
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i])
      cuts_[j++]=cuts_[i];
  }
  numberCuts_ = j;
}

// delete cuts
void
OsiNodeInfo::deleteCuts(int numberToDelete, int * which)
{
  int i;
  for (i=0;i<numberToDelete;i++) {
    int iCut=which[i];
    int number = cuts_[iCut]->decrement();
    if (!number)
      delete cuts_[iCut];
    cuts_[iCut]=NULL;
  }
  int j=0;
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i])
      cuts_[j++]=cuts_[i];
  }
  numberCuts_ = j;
}

// Apply cuts
void 
OsiNodeInfo::applyCuts(OsiSolverInterface * model) const
{
  if (numberCuts_) {
    // we can treat as OsiRowCut
    OsiRowCut * cuts = new OsiRowCut[numberCuts_];
    int i;
    for (i=0;i<numberCuts_;i++) {
      cuts[i]= OsiRowCut(*cuts_[i]);
    }
    model->applyRowCuts(numberCuts_,cuts);
    delete [] cuts;
  }
}

/* Node information complete class. 
   This is a complete description of status of model and it includes 
   all cuts so far*/

class OsiFullNodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  OsiFullNodeInfo ();

  /** Constructor from current state (and list of integers).
      adds in all cuts from chain starting at lastNode
  */
  OsiFullNodeInfo (OsiSolverInterface &model,
	   int numberIntegers, int * integer,
		   OsiNodeInfo * lastNode=NULL);
  
  // Copy constructor 
  OsiFullNodeInfo ( const OsiFullNodeInfo &);
   
  // Destructor 
  ~OsiFullNodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Full basis 
  OsiWarmStartBasis basis_;
  int numberIntegers_;
  // Bounds stored in full (for integers)
  int * lower_;
  int * upper_;
};

OsiFullNodeInfo::OsiFullNodeInfo() :
  OsiNodeInfo(),
  basis_(OsiWarmStartBasis()),
  numberIntegers_(0),
  lower_(NULL),
  upper_(NULL)
{
}
OsiFullNodeInfo::OsiFullNodeInfo(OsiSolverInterface & model,
				 int numberIntegers, int * integer,
				 OsiNodeInfo * lastNode) :
  OsiNodeInfo()
{
  const OsiWarmStartBasis* ws =
    dynamic_cast<const OsiWarmStartBasis*>(model.getWarmStart());
  
  assert (ws!=NULL); // make sure not volume
  basis_ = OsiWarmStartBasis(*ws);
  delete ws;
  numberIntegers_=numberIntegers;
  lower_ = new int [numberIntegers_];
  upper_ = new int [numberIntegers_];
  const double * lower = model.getColLower();
  const double * upper = model.getColUpper();
  int i;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    lower_[i]=(int)lower[iColumn];
    upper_[i]=(int)upper[iColumn];
    /*printf("%d lo %g, up %g, value %g\n",iColumn,
      lower[iColumn],upper[iColumn],model.getColSolution()[iColumn]);*/
  }
  // add in all previous cuts
  while(lastNode) {
    addCuts(lastNode->numberCuts(),lastNode->cuts());
    lastNode = lastNode->parent();
  }
}

OsiFullNodeInfo::OsiFullNodeInfo(const OsiFullNodeInfo & rhs) :
  OsiNodeInfo(rhs)
{  
  basis_=rhs.basis_;
  numberIntegers_=rhs.numberIntegers_;
  lower_=NULL;
  upper_=NULL;
  if (rhs.lower_!=NULL) {
    lower_ = new int [numberIntegers_];
    upper_ = new int [numberIntegers_];
    assert (upper_!=NULL);
    memcpy(lower_,rhs.lower_,numberIntegers_*sizeof(int));
    memcpy(upper_,rhs.upper_,numberIntegers_*sizeof(int));
  }
}

OsiNodeInfo * 
OsiFullNodeInfo::clone() const
{ 
  return (new OsiFullNodeInfo(*this));
}

OsiFullNodeInfo::~OsiFullNodeInfo ()
{
  delete [] lower_;
  delete [] upper_;
}

OsiNodeInfo * 
OsiFullNodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<numberIntegers_;i++) {
      int iColumn=integers[i];
      model->setColBounds( iColumn,lower_[i],upper_[i]);
    }
    // move basis
    basis=basis_;
    // apply cuts
    applyCuts(model);
  }
  assert(!parent_);
  return NULL;
}
/* If going virtual - may as well specialize */

/* Node information class when == 4 words change in basis and
   single bound
 */

class Osi4NodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  Osi4NodeInfo ();

  // Constructor from current state 
  Osi4NodeInfo (OsiNodeInfo * parent,
		int variable,
		int boundChange, 
		const unsigned int * whichBasis,
		const unsigned int * changeBasis);
  
  // Copy constructor 
  Osi4NodeInfo ( const Osi4NodeInfo &);
   
  // Destructor 
  ~Osi4NodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Which parts of basis (top bit set for rows)
  unsigned int which_[4];
  /// Values of basis
  unsigned int basis_[4];
  /// Which variable (top bit if upper bound changing)
  int variable_;
  // New bound
  int newBound_;
};

Osi4NodeInfo::Osi4NodeInfo() :
  OsiNodeInfo()
{
  memset(which_,0,4*sizeof(unsigned int));
  memset(basis_,0,4*sizeof(unsigned int));
  variable_=-1;
  newBound_=-1;
}
// Constructor from current state 
Osi4NodeInfo::Osi4NodeInfo (OsiNodeInfo * parent,
			    int variable,
			    int boundChange, 
			    const unsigned int * whichBasis,
			    const unsigned int * changeBasis)
 : OsiNodeInfo(parent)
{
  int i;
  for (i=0;i<4;i++) {
    which_[i]=whichBasis[i];
    basis_[i]=changeBasis[i];
  }
  variable_ = variable;
  newBound_=boundChange;
}

Osi4NodeInfo::Osi4NodeInfo(const Osi4NodeInfo & rhs) :
  OsiNodeInfo(rhs.parent_)
{  
  memcpy(which_,rhs.which_,4*sizeof(unsigned int));
  memcpy(basis_,rhs.basis_,4*sizeof(unsigned int));
  variable_=rhs.variable_;
  newBound_=rhs.newBound_;
}

OsiNodeInfo * 
Osi4NodeInfo::clone() const
{ 
  return (new Osi4NodeInfo(*this));
}


Osi4NodeInfo::~Osi4NodeInfo ()
{
}

OsiNodeInfo * 
Osi4NodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<4;i++) {
      unsigned int which = which_[i];
      unsigned int basisWord = basis_[i];
      unsigned int * update;
      if ((which&0x80000000)==0) {
	// structural
	update = (unsigned int *) basis.getStructuralStatus();
      } else {
	// row
	update = (unsigned int *) basis.getArtificialStatus();
	which &= 0x7fffffff;
      }
      update[which]=basisWord;
    }
    if ((variable_&0x80000000)==0) {
      // lower bound changing
      model->setColLower(variable_,newBound_);
    } else {
      // upper bound changing
      model->setColUpper(variable_&0x7fffffff,newBound_);
    }
    // apply cuts
    applyCuts(model);
    return NULL;
  } else {
    return parent_;
  }
}

/* Node information class when == 3 words change in basis and
   single bound
 */

class Osi3NodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  Osi3NodeInfo ();

  // Constructor from current state 
  Osi3NodeInfo (OsiNodeInfo * parent,
		int variable,
		int boundChange, 
		const unsigned int * whichBasis,
		const unsigned int * changeBasis);
  
  // Copy constructor 
  Osi3NodeInfo ( const Osi3NodeInfo &);
   
  // Destructor 
  ~Osi3NodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Which parts of basis (top bit set for rows)
  unsigned int which_[3];
  /// Values of basis
  unsigned int basis_[3];
  /// Which variable (top bit if upper bound changing)
  int variable_;
  // New bound
  int newBound_;
};

Osi3NodeInfo::Osi3NodeInfo() :
  OsiNodeInfo()
{
  memset(which_,0,3*sizeof(unsigned int));
  memset(basis_,0,3*sizeof(unsigned int));
  variable_=-1;
  newBound_=-1;
}
// Constructor from current state 
Osi3NodeInfo::Osi3NodeInfo (OsiNodeInfo * parent,
			    int variable,
			    int boundChange, 
			    const unsigned int * whichBasis,
			    const unsigned int * changeBasis)
 : OsiNodeInfo(parent)
{
  int i;
  for (i=0;i<3;i++) {
    which_[i]=whichBasis[i];
    basis_[i]=changeBasis[i];
  }
  variable_ = variable;
  newBound_=boundChange;
}

Osi3NodeInfo::Osi3NodeInfo(const Osi3NodeInfo & rhs) :
  OsiNodeInfo(rhs.parent_)
{  
  memcpy(which_,rhs.which_,3*sizeof(unsigned int));
  memcpy(basis_,rhs.basis_,3*sizeof(unsigned int));
  variable_=rhs.variable_;
  newBound_=rhs.newBound_;
}

OsiNodeInfo * 
Osi3NodeInfo::clone() const
{ 
  return (new Osi3NodeInfo(*this));
}


Osi3NodeInfo::~Osi3NodeInfo ()
{
}

OsiNodeInfo * 
Osi3NodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<3;i++) {
      unsigned int which = which_[i];
      unsigned int basisWord = basis_[i];
      unsigned int * update;
      if ((which&0x80000000)==0) {
	// structural
	update = (unsigned int *) basis.getStructuralStatus();
      } else {
	// row
	update = (unsigned int *) basis.getArtificialStatus();
	which &= 0x7fffffff;
      }
      update[which]=basisWord;
    }
    if ((variable_&0x80000000)==0) {
      // lower bound changing
      model->setColLower(variable_,newBound_);
    } else {
      // upper bound changing
      model->setColUpper(variable_&0x7fffffff,newBound_);
    }
    // apply cuts
    applyCuts(model);
    return NULL;
  } else {
    return parent_;
  }
}

/* Node information class when == 2 words change in basis and
   single bound
 */

class Osi2NodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  Osi2NodeInfo ();

  // Constructor from current state 
  Osi2NodeInfo (OsiNodeInfo * parent,
		int variable,
		int boundChange, 
		const unsigned int * whichBasis,
		const unsigned int * changeBasis);
  
  // Copy constructor 
  Osi2NodeInfo ( const Osi2NodeInfo &);
   
  // Destructor 
  ~Osi2NodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Which parts of basis (top bit set for rows)
  unsigned int which_[2];
  /// Values of basis
  unsigned int basis_[2];
  /// Which variable (top bit if upper bound changing)
  int variable_;
  // New bound
  int newBound_;
};

Osi2NodeInfo::Osi2NodeInfo() :
  OsiNodeInfo()
{
  memset(which_,0,2*sizeof(unsigned int));
  memset(basis_,0,2*sizeof(unsigned int));
  variable_=-1;
  newBound_=-1;
}
// Constructor from current state 
Osi2NodeInfo::Osi2NodeInfo (OsiNodeInfo * parent,
			    int variable,
			    int boundChange, 
			    const unsigned int * whichBasis,
			    const unsigned int * changeBasis)
 : OsiNodeInfo(parent)
{
  int i;
  for (i=0;i<2;i++) {
    which_[i]=whichBasis[i];
    basis_[i]=changeBasis[i];
  }
  variable_ = variable;
  newBound_=boundChange;
}

Osi2NodeInfo::Osi2NodeInfo(const Osi2NodeInfo & rhs) :
  OsiNodeInfo(rhs.parent_)
{  
  memcpy(which_,rhs.which_,2*sizeof(unsigned int));
  memcpy(basis_,rhs.basis_,2*sizeof(unsigned int));
  variable_=rhs.variable_;
  newBound_=rhs.newBound_;
}

OsiNodeInfo * 
Osi2NodeInfo::clone() const
{ 
  return (new Osi2NodeInfo(*this));
}


Osi2NodeInfo::~Osi2NodeInfo ()
{
}

OsiNodeInfo * 
Osi2NodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<2;i++) {
      unsigned int which = which_[i];
      unsigned int basisWord = basis_[i];
      unsigned int * update;
      if ((which&0x80000000)==0) {
	// structural
	update = (unsigned int *) basis.getStructuralStatus();
      } else {
	// row
	update = (unsigned int *) basis.getArtificialStatus();
	which &= 0x7fffffff;
      }
      update[which]=basisWord;
    }
    if ((variable_&0x80000000)==0) {
      // lower bound changing
      model->setColLower(variable_,newBound_);
    } else {
      // upper bound changing
      model->setColUpper(variable_&0x7fffffff,newBound_);
    }
    // apply cuts
    applyCuts(model);
    return NULL;
  } else {
    return parent_;
  }
}

/* Node information class when == 1 words change in basis and
   single bound
 */

class Osi1NodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  Osi1NodeInfo ();

  // Constructor from current state 
  Osi1NodeInfo (OsiNodeInfo * parent,
		int variable,
		int boundChange, 
		const unsigned int * whichBasis,
		const unsigned int * changeBasis);
  
  // Copy constructor 
  Osi1NodeInfo ( const Osi1NodeInfo &);
   
  // Destructor 
  ~Osi1NodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Which parts of basis (top bit set for rows)
  unsigned int which_[1];
  /// Values of basis
  unsigned int basis_[1];
  /// Which variable (top bit if upper bound changing)
  int variable_;
  // New bound
  int newBound_;
};

Osi1NodeInfo::Osi1NodeInfo() :
  OsiNodeInfo()
{
  memset(which_,0,1*sizeof(unsigned int));
  memset(basis_,0,1*sizeof(unsigned int));
  variable_=-1;
  newBound_=-1;
}
// Constructor from current state 
Osi1NodeInfo::Osi1NodeInfo (OsiNodeInfo * parent,
			    int variable,
			    int boundChange, 
			    const unsigned int * whichBasis,
			    const unsigned int * changeBasis)
 : OsiNodeInfo(parent)
{
  int i;
  for (i=0;i<1;i++) {
    which_[i]=whichBasis[i];
    basis_[i]=changeBasis[i];
  }
  variable_ = variable;
  newBound_=boundChange;
}

Osi1NodeInfo::Osi1NodeInfo(const Osi1NodeInfo & rhs) :
  OsiNodeInfo(rhs.parent_)
{  
  memcpy(which_,rhs.which_,1*sizeof(unsigned int));
  memcpy(basis_,rhs.basis_,1*sizeof(unsigned int));
  variable_=rhs.variable_;
  newBound_=rhs.newBound_;
}

OsiNodeInfo * 
Osi1NodeInfo::clone() const
{ 
  return (new Osi1NodeInfo(*this));
}


Osi1NodeInfo::~Osi1NodeInfo ()
{
}

OsiNodeInfo * 
Osi1NodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<1;i++) {
      unsigned int which = which_[i];
      unsigned int basisWord = basis_[i];
      unsigned int * update;
      if ((which&0x80000000)==0) {
	// structural
	update = (unsigned int *) basis.getStructuralStatus();
      } else {
	// row
	update = (unsigned int *) basis.getArtificialStatus();
	which &= 0x7fffffff;
      }
      update[which]=basisWord;
    }
    if ((variable_&0x80000000)==0) {
      // lower bound changing
      model->setColLower(variable_,newBound_);
    } else {
      // upper bound changing
      model->setColUpper(variable_&0x7fffffff,newBound_);
    }
    // apply cuts
    applyCuts(model);
    return NULL;
  } else {
    return parent_;
  }
}

/* Node information class when words change in basis and
    bounds change
 */

class OsiPartialNodeInfo : public OsiNodeInfo {
public:

  virtual OsiNodeInfo * applyToModel(OsiSolverInterface * model,int mode,
				     const int * integers,
				     OsiWarmStartBasis & basis) const;
  // Default Constructor 
  OsiPartialNodeInfo ();

  // Constructor from current state 
  OsiPartialNodeInfo (OsiNodeInfo * parent,
		int numberChangedBounds,const int * variables,
		const int * boundChanges, int numberChanged,
		const unsigned int * whichBasis,
		const unsigned int * changeBasis);
  
  // Copy constructor 
  OsiPartialNodeInfo ( const OsiPartialNodeInfo &);
   
  // Destructor 
  ~OsiPartialNodeInfo ();
  
  /// Clone
  virtual OsiNodeInfo * clone() const;
  // Public data
  // Which parts of basis (top bit set for rows)
  unsigned int * which_;
  /// Values of basis
  unsigned int * basis_;
  /// Number of basis changes
  int numberChanged_;
  /// Which variable (top bit if upper bound changing)
  int * variables_;
  // New bound
  int * newBounds_;
  /// Number of bound changes
  int numberChangedBounds_;
};

OsiPartialNodeInfo::OsiPartialNodeInfo() :
  OsiNodeInfo(),
  which_(NULL),
  basis_(NULL),
  numberChanged_(0),
  variables_(NULL),
  newBounds_(NULL),
  numberChangedBounds_(0)
{
}
// Constructor from current state 
OsiPartialNodeInfo::OsiPartialNodeInfo (OsiNodeInfo * parent,
			    int numberChangedBounds,const int * variables,
			    const int * boundChanges, int numberChanged,
			    const unsigned int * whichBasis,
			    const unsigned int * changeBasis)
 : OsiNodeInfo(parent)
{
  numberChanged_ = numberChanged;
  which_ = new unsigned int [numberChanged_];
  basis_ = new unsigned int [numberChanged_];
  int i;
  for (i=0;i<numberChanged_;i++) {
    which_[i]=whichBasis[i];
    basis_[i]=changeBasis[i];
  }

  numberChangedBounds_ = numberChangedBounds;
  variables_ = new int [numberChangedBounds_];
  newBounds_ = new int [numberChangedBounds_];

  for (i=0;i<numberChangedBounds_;i++) {
    variables_[i]=variables[i];
    newBounds_[i]=boundChanges[i];
  }
}

OsiPartialNodeInfo::OsiPartialNodeInfo(const OsiPartialNodeInfo & rhs) :
  OsiNodeInfo(rhs.parent_)
{  
  numberChanged_ = rhs.numberChanged_;
  which_ = new unsigned int [numberChanged_];
  basis_ = new unsigned int [numberChanged_];
  int i;
  for (i=0;i<numberChanged_;i++) {
    which_[i]=rhs.which_[i];
    basis_[i]=rhs.basis_[i];
  }

  numberChangedBounds_ = rhs.numberChangedBounds_;
  variables_ = new int [numberChangedBounds_];
  newBounds_ = new int [numberChangedBounds_];

  for (i=0;i<numberChangedBounds_;i++) {
    variables_[i]=rhs.variables_[i];
    newBounds_[i]=rhs.newBounds_[i];
  }
}

OsiNodeInfo * 
OsiPartialNodeInfo::clone() const
{ 
  return (new OsiPartialNodeInfo(*this));
}


OsiPartialNodeInfo::~OsiPartialNodeInfo ()
{
  delete [] which_;
  delete [] basis_;
  delete [] variables_;
  delete [] newBounds_;
}

OsiNodeInfo * 
OsiPartialNodeInfo::applyToModel(OsiSolverInterface * model,int mode,
			      const int * integers,
			      OsiWarmStartBasis & basis) const
{
  if (mode) {
    // branch - do bounds
    int i;
    for (i=0;i<numberChanged_;i++) {
      unsigned int which = which_[i];
      unsigned int basisWord = basis_[i];
      unsigned int * update;
      if ((which&0x80000000)==0) {
	// structural
	update = (unsigned int *) basis.getStructuralStatus();
      } else {
	// row
	update = (unsigned int *) basis.getArtificialStatus();
	which &= 0x7fffffff;
      }
      update[which]=basisWord;
    }
    for (i=0;i<numberChangedBounds_;i++) {
      int variable = variables_[i];
      if ((variable&0x80000000)==0) {
	// lower bound changing
	model->setColLower(variable,newBounds_[i]);
      } else {
	// upper bound changing
	model->setColUpper(variable&0x7fffffff,newBounds_[i]);
      }
    }
    // apply cuts
    applyCuts(model);
    return NULL;
  } else {
    return parent_;
  }
}

// Trivial class for Branch and Bound

class OsiNode  {
 
public:
    
  // Default Constructor 
  OsiNode ();

  // Constructor from current state (and list of integers)
  // Also chooses branching variable (if none set to -1)
  // Strategy - 1 branch on largest min(up,down)
  //          - 0 branch opposite largest max(up,down)
  OsiNode (OsiSolverInterface * model,
	   OsiNode * lastNode,
	   int numberIntegers, int * integer,
	   OsiWarmStartBasis & lastws,
	   const int * lastLower, const int * lastUpper,
	   int strategy);
  
  // Copy constructor 
  OsiNode ( const OsiNode &);
   
  // Assignment operator 
  OsiNode & operator=( const OsiNode& rhs);

  // Destructor 
  ~OsiNode ();
  
  // Public data
  // Information to make basis and bounds
  OsiNodeInfo * nodeInfo_;
  OsiNodeInfo * nodeInfo() const
  {return nodeInfo_;};
  // Objective value
  double objectiveValue_;
  // Branching variable (0 is first integer)
  int variable_;
  // Way to branch - -1 down (first), 1 up, -2 down (second), 2 up (second)
  int way_;
  // Current value
  double value_;
  // Depth
  int depth_;
  // Number unsatisfied
  int numberUnsatisfied_;
};


OsiNode::OsiNode() :
  nodeInfo_(NULL),
  objectiveValue_(1.0e100),
  variable_(-100),
  way_(-1),
  value_(0.5),
  depth_(-1),
  numberUnsatisfied_(0)
{
}
OsiNode::OsiNode(OsiSolverInterface * model,
		 OsiNode * lastNode,
		 int numberIntegers, int * integer,
		 OsiWarmStartBasis & lastws,
		 const int * lastLower, const int * lastUpper,
		 int strategy) :
  nodeInfo_(NULL),
  objectiveValue_(1.0e100),
  variable_(-1),
  way_(-1),
  value_(0.5),
  depth_(-1),
  numberUnsatisfied_(0)
{
  
  const OsiWarmStartBasis* ws =
    dynamic_cast<const OsiWarmStartBasis*>(model->getWarmStart());
  
  assert (ws!=NULL); // make sure not volume
  objectiveValue_ = model->getObjSense()*model->getObjValue();
  int * lowerNow = new int [numberIntegers];
  int * upperNow = new int [numberIntegers];
  assert (upperNow!=NULL);
  const double * lower = model->getColLower();
  const double * upper = model->getColUpper();
  const double * solution = model->getColSolution();
  const double * reducedCost = model->getReducedCost();
  // Here could go cuts etc etc
  // For now just fix on reduced costs
  double bestObjective;
  model->getDblParam(OsiDualObjectiveLimit,bestObjective);
  double gap = bestObjective - model->getObjValue();
  double tolerance = 1.0e-6;
  int i;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    if (upper[iColumn]-lower[iColumn]>tolerance) {
      if (solution[iColumn]<lower[iColumn]+tolerance&&
	  reducedCost[iColumn]>gap)
	model->setColUpper(iColumn,lower[iColumn]);
      else if (solution[iColumn]>upper[iColumn]-tolerance&&
	       -reducedCost[iColumn]>gap)
	model->setColLower(iColumn,upper[iColumn]);
    }
    lowerNow[i]=(int)lower[iColumn];
    upperNow[i]=(int)upper[iColumn];
  }
  // Probing cuts
#ifdef CUTS
  CglProbing try2;
#endif
  OsiCuts cuts;
#ifdef CUTS
  try2.generateCuts(*model,cuts);
#endif
  int violated = 0;
  int numberColumnCuts = cuts.sizeColCuts();
  if (numberColumnCuts) {
    int i;
    for (i=0;i<numberColumnCuts;i++) {
      const OsiColCut * thisCut = cuts.colCutPtr(i);
      const OsiPackedVector & lbs = thisCut->lbs();
      const OsiPackedVector & ubs = thisCut->ubs();
      int j;
      int n;
      const int * which;
      const double * values;
      n = lbs.getNumElements();
      which = lbs.getIndices();
      values = lbs.getElements();
      for (j=0;j<n;j++) {
	int iColumn = which[j];
	double value=solution[iColumn];
	double current = lower[iColumn];
	assert(values[j]>current);
	model->setColLower(iColumn,values[j]);
	if (value<values[j]-1.0e-5)
	  violated = 1;
	if (values[j]>upper[iColumn]+1.0e-5) {
	  // infeasible
	  violated = 2;
	  objectiveValue_ = 1.0e100;
	  break;
	}
      }
      n = ubs.getNumElements();
      which = ubs.getIndices();
      values = ubs.getElements();
      for (j=0;j<n;j++) {
	int iColumn = which[j];
	double value=solution[iColumn];
	double current = upper[iColumn];
	assert(values[j]<current);
	model->setColUpper(iColumn,values[j]);
	if (value>values[j]+1.0e-5)
	  violated = 1;
	if (values[j]<lower[iColumn]-1.0e-5) {
	  // infeasible
	  violated = 2;
	  objectiveValue_ = 1.0e100;
	  break;
	}
      }
    }
    if (violated) {
      // return and re-solve
      delete [] lowerNow;
      delete [] upperNow;
      delete ws;
      if (violated==1)
	variable_ = numberIntegers;
      return;
    }
  }
  int numberArtificials = lastws.getNumArtificial();
  int numberStructurals = lastws.getNumStructural();
  int numberArtificialsNow = ws->getNumArtificial();
  if (!numberArtificials) {
    // full basis
    nodeInfo_=new OsiFullNodeInfo(*model,numberIntegers,integer);
  } else {
    // arrays for changes
    int * boundChanges = new int [ 2*numberIntegers];
    int * variables = new int [2*numberIntegers];
    // compute how big full and partial ones will be (approximately)
    int fullSize = ((numberArtificials+numberStructurals)>>4)
      + 2*numberIntegers;
    int partialSize=0;
    int maxBasisLength = ((max(numberArtificials,numberArtificialsNow)+15)>>4)+
      ((numberStructurals+15)>>4);
    unsigned int * whichBasis = new unsigned int [maxBasisLength]; 
    unsigned int * changeBasis = new unsigned int [maxBasisLength]; 
    // check number of changes
    int numberChanged = 0;
    int number;
    const unsigned int * old = 
      (const unsigned int *) lastws.getArtificialStatus();
    const unsigned int * now = 
      (const unsigned int *) ws->getArtificialStatus();
    //printf("rows now %d, was %d\n",numberArtificialsNow,numberArtificials);
    if (numberArtificials>=numberArtificialsNow) {
      number=(numberArtificialsNow+15)>>4;;
      for (i=0;i<number;i++) { 
	if (old[i]!=now[i]) {
	  whichBasis[numberChanged]=i|0x80000000;
	  changeBasis[numberChanged++]=now[i];
	}
      }
    } else {
      number=(numberArtificials+15)>>4;;
      for (i=0;i<number;i++) { 
	if (old[i]!=now[i]) {
	  whichBasis[numberChanged]=i|0x80000000;
	  changeBasis[numberChanged++]=now[i];
	}
      }
    }
    number=(numberStructurals+15)>>4;;
    old = (const unsigned int *) lastws.getStructuralStatus();
    now = (const unsigned int *) ws->getStructuralStatus();
    for (i=0;i<number;i++) { 
      if (old[i]!=now[i]) {
	whichBasis[numberChanged]=i;
	changeBasis[numberChanged++]=now[i];
      }
    }
    partialSize = 2*numberChanged;
    int numberChangedBounds=0;
    for (i=0;i<numberIntegers;i++) {
      if (lowerNow[i]!=lastLower[i]) {
	variables[numberChangedBounds]=integer[i];
	boundChanges[numberChangedBounds++]=lowerNow[i];
      }
      if (upperNow[i]!=lastUpper[i]) {
	variables[numberChangedBounds]=integer[i]|0x80000000;
	boundChanges[numberChangedBounds++]=upperNow[i];
      }
#ifdef OSI_DEBUG
      int iColumn = integer[i];
      if (lowerNow[i]!=lastLower[i])
	printf("lower on %d %d changed from %d to %d\n",
	       i,iColumn,lastLower[i],lowerNow[i]);
      if (upperNow[i]!=lastUpper[i])
	printf("upper on %d %d changed from %d to %d\n",
	       i,iColumn,lastUpper[i],upperNow[i]);
#endif
    }
    partialSize += 2*numberChangedBounds;
#ifdef OSI_DEBUG
    printf("%d changed bounds, %d changed basis\n",numberChangedBounds,
	   numberChanged);
#endif
    assert (numberChangedBounds&&numberChanged);
    if (numberChangedBounds==1&&numberChanged<=4) {
      // partial basis
      switch (numberChanged) {
      case 1:
	nodeInfo_=new Osi1NodeInfo(lastNode->nodeInfo_,
				   variables[0],boundChanges[0],
				   whichBasis,changeBasis);
	break;
      case 2:
	nodeInfo_=new Osi2NodeInfo(lastNode->nodeInfo_,
				   variables[0],boundChanges[0],
				   whichBasis,changeBasis);
	break;
      case 3:
	nodeInfo_=new Osi3NodeInfo(lastNode->nodeInfo_,
				   variables[0],boundChanges[0],
				   whichBasis,changeBasis);
	break;
      case 4:
	nodeInfo_=new Osi4NodeInfo(lastNode->nodeInfo_,
				   variables[0],boundChanges[0],
				   whichBasis,changeBasis);
	break;
      }
    } else if (partialSize*3<fullSize*2) {
      // partial
      nodeInfo_ = new OsiPartialNodeInfo (lastNode->nodeInfo_,
					  numberChangedBounds, variables,
					  boundChanges, numberChanged,
					  whichBasis,
					  changeBasis);
    } else {
      // full basis
      nodeInfo_=new OsiFullNodeInfo(*model,numberIntegers,integer,
				    lastNode->nodeInfo_);
    }
    delete [] boundChanges;
    delete [] variables;
    delete [] whichBasis;
    delete [] changeBasis;
  }
  // Gomory cuts as test
#ifdef CUTS
  CglGomory try1;
  try1.generateCuts(*model,cuts);
#endif
  nodeInfo_->addCuts(cuts);
  // Hard coded integer tolerance
#define INTEGER_TOLERANCE 1.0e-6
  variable_=-1;
  // This has hard coded integer tolerance
  double mostAway=INTEGER_TOLERANCE;
  numberUnsatisfied_ = 0;
  if (lastNode)
    depth_ = lastNode->depth_+1;
  else
    depth_ = 0;
  // Number of strong branching candidates
#define STRONG_BRANCHING 5
#ifdef STRONG_BRANCHING
  double upMovement[STRONG_BRANCHING];
  double downMovement[STRONG_BRANCHING];
  double solutionValue[STRONG_BRANCHING];
  int chosen[STRONG_BRANCHING];
  int iSmallest=0;
  // initialize distance from integer
  for (i=0;i<STRONG_BRANCHING;i++) {
    upMovement[i]=0.0;
    chosen[i]=-1;
  }
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    double value = solution[iColumn];
    value = max(value, lower[i]);
    value = min(value, upper[i]);
    double nearest = floor(value+0.5);
    if (fabs(value-nearest)>INTEGER_TOLERANCE) {
      numberUnsatisfied_++;
      if (fabs(value-nearest)>mostAway) {
	double away = fabs(value-nearest);
	if (away>upMovement[iSmallest]) {
	  //add to list
	  upMovement[iSmallest]=away;
	  solutionValue[iSmallest]=value;
	  chosen[iSmallest]=i;
	  int j;
	  iSmallest=-1;
	  double smallest = 1.0;
	  for (j=0;j<STRONG_BRANCHING;j++) {
	    if (upMovement[j]<smallest) {
	      smallest=upMovement[j];
	      iSmallest=j;
	    }
	  }
	}
      }
    }
  }
  int numberStrong=0;
  for (i=0;i<STRONG_BRANCHING;i++) {
    if (chosen[i]>=0) { 
      solutionValue[numberStrong]=solutionValue[i];
      chosen[numberStrong ++] = chosen[i];
    }
  }
  int trueSequence[STRONG_BRANCHING];
  for (i=0;i<numberStrong;i++) {
    int iInt = chosen[i];
    int iColumn = integer[iInt];
    trueSequence[i]=iColumn;
    double value = solutionValue[i]; // value of variable in original
    downMovement[i] = floor(value);
    upMovement[i]=ceil(value);
  }
  OsiClpSolverInterface * clpptr = 
    dynamic_cast< OsiClpSolverInterface*>(model);
  assert(clpptr);
  //clpptr=NULL;
  int returnCode=0;
  if (clpptr) {
    // Clp so use native method
    ClpSimplex simplex(*clpptr->getModelPtr());
    simplex.setBasis(*ws);
    returnCode = 
      simplex.strongBranching(numberStrong,trueSequence,
			      upMovement,downMovement,
			      false,false);
  } else {
    // Not Clp
    model->markHotStart();
    for (i=0;i<numberStrong;i++) {
      int iColumn = trueSequence[i];
      double objectiveChange;

      // try down
      double saveValue = upper[iColumn];
      model->setColUpper(iColumn,downMovement[i]);
      model->solveFromHotStart();
      model->setColUpper(iColumn,saveValue);
      if (model->isProvenOptimal()&&!model->isDualObjectiveLimitReached()) {
	objectiveChange = model->getObjSense()*model->getObjValue()
	  - objectiveValue_;
      } else {
	objectiveChange = 1.0e100;
      }
      downMovement[i]=objectiveChange;

      // try up
      saveValue = lower[iColumn];
      model->setColLower(iColumn,upMovement[i]);
      model->solveFromHotStart();
      model->setColLower(iColumn,saveValue);
      if (model->isProvenOptimal()&&!model->isDualObjectiveLimitReached()) {
	objectiveChange = model->getObjSense()*model->getObjValue()
	  - objectiveValue_;
      } else {
	objectiveChange = 1.0e100;
      }
      upMovement[i]=objectiveChange;

      /* Possibilities are:
	 Both sides feasible 
	 Neither side feasible - exit
	 One side feasible - exit
      */
      if (upMovement[i]<1.0e100) {
	if(downMovement[i]<1.0e100) {
	  // feasible - no action
	} else {
	  // up feasible, down infeasible
	  returnCode = 1;
	  break;
	}
      } else {
	if(downMovement[i]<1.0e100) {
	  // down feasible, up infeasible
	  returnCode = 1;
	  break;
	} else {
	  // neither side feasible
	  returnCode = 1;
	  break;
	}
      }
    }
    // Delete the snapshot
    model->unmarkHotStart();
  }
  if (returnCode<0) {
    // infeasible both ways
    objectiveValue_=1.0e100;
  } else if (returnCode>0) {
    // infeasible one way - find first
    for (i=0;i<numberStrong;i++) {
      if (downMovement[i]>1.0e50) {
	// set new lower
	int iColumn = trueSequence[i];
	model->setColLower(iColumn,ceil(solutionValue[i]));
	break;
      } else if (upMovement[i]>1.0e50) {
	// set new upper
	int iColumn = trueSequence[i];
	model->setColUpper(iColumn,floor(solutionValue[i]));
	break;
      }
    }
    variable_ = numberIntegers;
  } else {
    // all feasible - choose best bet
    // choose the one that makes most difference both ways
    double best = -1.0;
    double best2 = -1.0;
    for (i=0;i<numberStrong;i++) {
      int iInt = chosen[i];
      model->messageHandler()->message(OSI_BAB_STRONG,model->messages())
	<<i<<iInt<<downMovement[i]<<upMovement[i]
	<<solutionValue[i]
	<<OsiMessageEol;
      bool better = false;
      double up = upMovement[i];
      double down = downMovement[i];
      double value = solutionValue[i];
      // Modify a bit
      down += 0.01*(value-floor(value));
      up += 0.01*(ceil(value)-value);
      double test;
      if (!strategy)
	test = max(up,down); // before solution
      else
	test=min(up,down);
      if (test>best) {
	// smaller is better
	better=true;
      } else if (test>best-1.0e-5) {
	if (up+down>best2+1.0e-5) {
	  // smaller is about same, but larger is better
	  better=true;
	}
      }
      if (better) {
	best = test;
	best2 = up+down;
	variable_ = iInt;
	double value = solutionValue[i];
	value = max(value,(double) lowerNow[variable_]);
	value = min(value,(double) upperNow[variable_]);
	value_=value;
	if (up<=down)
	  way_=1; // up
	else
	  way_=-1; // down
      }
    }
  }
#else
  // not strong branching
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integer[i];
    double value = solution[iColumn];
    value = max(value, lower[i]);
    value = min(value, upper[i]);
    double nearest = floor(value+0.5);
    double distance = fabs(value-nearest);
    if (distance>INTEGER_TOLERANCE) {
      numberUnsatisfied_++;
      // If before solution - opposite worst
      distance = 1.0-distance;
      if (distance>mostAway) {
	mostAway=distance;
	variable_=i;
	value_=value;
	if (value<=nearest)
	  way_=1; // up
	else
	  way_=-1; // down
      }
    }
  }
#endif
  delete [] lowerNow;
  delete [] upperNow;
  delete ws;
}

OsiNode::OsiNode(const OsiNode & rhs) 
{  
  if (rhs.nodeInfo_)
    nodeInfo_ = rhs.nodeInfo_->clone();
  else
    nodeInfo_=NULL;
  objectiveValue_=rhs.objectiveValue_;
  variable_=rhs.variable_;
  way_=rhs.way_;
  value_=rhs.value_;
  depth_ = rhs.depth_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
}

OsiNode &
OsiNode::operator=(const OsiNode & rhs)
{
  if (this != &rhs) {
    delete nodeInfo_;
    if (nodeInfo_)
      nodeInfo_ = rhs.nodeInfo_->clone();
    else
      nodeInfo_ = NULL;
    objectiveValue_=rhs.objectiveValue_;
    variable_=rhs.variable_;
    way_=rhs.way_;
    value_=rhs.value_;
    depth_ = rhs.depth_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
  }
  return *this;
}


OsiNode::~OsiNode ()
{
  if (nodeInfo_&&!nodeInfo_->numberPointingToThis_)
    delete nodeInfo_;
}

#include <vector>

/*  These are alternative strategies for node traversal.  
    They can take data etc for fine tuning 
*/

class OsiCompareAny {
public:
  OsiCompareAny * test_;
  // Default Constructor 
  OsiCompareAny () {test_=NULL;};

  virtual ~OsiCompareAny() {};

  virtual bool test (OsiNode * x, OsiNode * y) {return true;};
  bool operator() (OsiNode * x, OsiNode * y) {
    return test(x,y);
  }
};
class OsiCompare {
public:
  OsiCompareAny * test_;
  // Default Constructor 
  OsiCompare () {test_=NULL;};

  virtual ~OsiCompare() {};

  bool operator() (OsiNode * x, OsiNode * y) {
    return test_->test(x,y);
  }
};
class OsiCompareDepth : public OsiCompareAny{
public:
  // Default Constructor 
  OsiCompareDepth () {test_=this;};

  ~OsiCompareDepth() {};
  virtual bool test (OsiNode * x, OsiNode * y) {
    return x->depth_ < y->depth_;
  }
};
class OsiCompareObjective  : public OsiCompareAny {
public:
  // Default Constructor 
  OsiCompareObjective () {test_=this;};

  ~OsiCompareObjective() {};
  virtual bool test (OsiNode * x, OsiNode * y) {
    return x->objectiveValue_ > y->objectiveValue_;
  }
};
// This is an example of a more complex rule with data
class OsiCompareExample  : public OsiCompareAny {
public:
  // Weight for each infeasibility
  double weight_;
  // Default Constructor 
  OsiCompareExample () : weight_(1.0) {test_=this;};
  // Constructor with weight
  OsiCompareExample (double weight) : weight_(weight) {test_=this;};

  ~OsiCompareExample() {};
  virtual bool test (OsiNode * x, OsiNode * y) {
    return x->objectiveValue_+ weight_*x->numberUnsatisfied_ > 
      y->objectiveValue_ + weight_*y->numberUnsatisfied_;
  }
};
// Heap/priority queue
template <class T, class OsiCompare>
class tree : public std::vector <T> {
  OsiCompare comparison_;
public:
  tree() {};
#if 1
  void setComparison(OsiCompareAny  &compare) {
    comparison_.test_ = &compare;
    make_heap(begin(), end(), comparison_);
  };
#endif
  const T& top() {return front();};
  void push(const T& x) {
    push_back(x);
    push_heap(begin(), end(), comparison_);
  };
  void pop() {
    pop_heap(begin(), end(), comparison_);
    pop_back();
  };
};
// Invoke solver's built-in enumeration algorithm
void 
OsiClpSolverInterface::branchAndBound() {
  // solve LP
  initialSolve();
  if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
    // This is a really simple Branch and Bound code - mainly
    // to test strong branching
    // I should look at STL more to allow other than depth first
    double bestObjective=1.0e50;
    int numberIntegers=0;
    int numberColumns = getNumCols();
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn))
	numberIntegers++;
    }
    if (!numberIntegers) {
      handler_->message(OSI_BAB_NOINT,messages_)
	<<OsiMessageEol;
      return;
    }
    // save number of rows
    int numberRowsAtContinuous = getNumRows();
    int * which = new int[numberIntegers]; // which variables are integer
    // Stuff so can do incremental nodes
    int * lowerBefore = new int[numberIntegers];
    int * upperBefore = new int[numberIntegers];
    OsiWarmStartBasis lastws;

    numberIntegers=0;
    const double * lower = getColLower();
    const double * upper = getColUpper();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( isInteger(iColumn)) {
	lowerBefore[numberIntegers]=(int) lower[iColumn];
	upperBefore[numberIntegers]=(int) upper[iColumn];
	which[numberIntegers++]=iColumn;
      }
    }
    activateRowCutDebugger("mod011");

    // empty tree
    tree <OsiNode *, OsiCompare> branchingTree;
    // various node selection criteria
    OsiCompareObjective compareObjective;
    OsiCompareDepth compareDepth;
    OsiCompareExample compareExample(50.0);
    // Start with depth first
    branchingTree.setComparison(compareDepth);
    int solutionFound=0;

    // Add continuous to it - forcing full node as lastws is empty;
    OsiNode * newNode = new OsiNode(this,NULL,
				    numberIntegers,which,lastws,
				    lowerBefore,upperBefore,
				    solutionFound);
    while (newNode->variable_==numberIntegers) {
      // we don't want to delete parents so do delayed delete
      OsiNode * save =newNode;
      resolve();
      if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
	newNode = new OsiNode(this,NULL,numberIntegers,which,lastws,
			      lowerBefore,upperBefore,solutionFound);
      } else {
	newNode=NULL;
	break;
      }
      delete save;
    }
    if (newNode&&newNode->objectiveValue_<bestObjective) 
      branchingTree.push(newNode); 

    // maximum depth for tree walkback
    int maximumDepth=10;
    OsiNodeInfo ** walkback = new OsiNodeInfo * [maximumDepth];
    // number of rows before each node application
    int * numberRowsNow = new int[maximumDepth];

    // For printing totals
    int numberIterations=0;
    int numberNodes =0;

    OsiFullNodeInfo * bestNode=NULL;
    // while until nothing on stack
    while (!branchingTree.empty()) {
      // last node
      OsiNode * node = branchingTree.top();
      branchingTree.pop();
      if (node->objectiveValue_ < bestObjective) {
	numberNodes++;
	int i;
	if (node->variable_>=0) {
	  int nNode=0;
	  OsiNodeInfo * nodeInfo = node->nodeInfo_;
	  while (nodeInfo) {
	    /* when working then just unwind until where new node
	       joins old node (for cuts?) */
	    //printf("nNode = %d, nodeInfo = %x\n",nNode,nodeInfo);
	    walkback[nNode++]=nodeInfo;
	    nodeInfo = nodeInfo->applyToModel(this,0,which,basis_);
	    if (nNode==maximumDepth) {
	      maximumDepth *= 2;
	      OsiNodeInfo ** temp = new OsiNodeInfo * [maximumDepth];
	      for (i=0;i<nNode;i++) 
		temp[i] = walkback[i];
	      delete [] walkback;
	      walkback = temp;
	      delete [] numberRowsNow;
	      numberRowsNow = new int[maximumDepth];
	    }
	  }
	  //int saveNumberNodes =nNode;
	  while (nNode) {
	    --nNode;
	    walkback[nNode]->applyToModel(this,1,which,basis_);
	    numberRowsNow[nNode] = getNumRows();
	  }
	  // now save situation for node differences
	  lastws = basis_;
	  const double * lower = getColLower();
	  const double * upper = getColUpper();
	  for (i=0;i<numberIntegers;i++) {
	    int iColumn = which[i];
	    lowerBefore[i]=(int)lower[iColumn];
	    upperBefore[i]=(int)upper[iColumn];
	  }
	  // do branching variable
	  bool deleteNode=false;
	  if (node->way_<0) {
	    setColUpper(which[node->variable_],floor(node->value_));
#ifdef OSI_DEBUG
	    printf("branching down on %d\n",node->variable_);
#endif
	    // now push back node if more to come
	    if (node->way_==-1) { 
	      node->way_=+2;	  // Swap direction
	      node->nodeInfo_->increment();
	      branchingTree.push(node);
	    } else {
	      deleteNode=true;
	    }
	  } else {
	    setColLower(which[node->variable_],ceil(node->value_));
#ifdef OSI_DEBUG
	    printf("branching up on %d\n",node->variable_);
#endif
	    // now push back node if more to come
	    if (node->way_==1) { 
	      node->way_=-2;	  // Swap direction
	      node->nodeInfo_->increment();
	      branchingTree.push(node);
	    } else {
	      deleteNode=true;
	    }
	  }
	  // solve
	  resolve();
	  numberIterations += getIterationCount();
	  if (!isIterationLimitReached()) {
	    if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
	      // We still have stack which applied cuts
	      // so go through finding which cuts can be dropped
	      /* NO NO - not that simple, anyway we would then
		 need to re-do incremental bases.  Do later as
		 garbage collection? - Ask Laci
	      */
#if 0
	      int numberRows = numberRowsAtContinuous;
	      int * drop = new int[getNumRows()-numberRows];
	      nNode = saveNumberNodes;
	      int nTotalDrop = 0;
	      while (nNode) {
		--nNode;
		int nDrop=0;
		for (i=numberRows;i<numberRowsNow[nNode];i++) {
		  if (basis_.getArtifStatus(i)==OsiWarmStartBasis::basic) {
		    drop[nTotalDrop+nDrop++]=i-numberRows;
		  }
		}
		if (nDrop) {
		  walkback[nNode]->deleteCuts(nDrop,drop+nTotalDrop);
		  // and modify indices for deleteRows
		  for (i=nTotalDrop;i<nTotalDrop+nDrop;i++) 
		    drop[i] += numberRows;
		  nTotalDrop += nDrop;
		}
		numberRows = numberRowsNow[nNode];
	      }
	      if (nTotalDrop) {
		// delete rows so basis will be correct
		deleteRows(nTotalDrop,drop);
		// should be optimal
		//printf("should be optimal\n");
		//resolve();
	      }
	      delete [] drop;
#endif
	      newNode = new OsiNode(this,node,numberIntegers,which,lastws,
				    lowerBefore,upperBefore,solutionFound);
	    } else {
	      newNode=NULL;
	    }
	    // something extra may have been fixed by strong branching
	    // if so go round again
	    while (newNode&&newNode->variable_==numberIntegers) {
	      // we don't want to delete parents so do delayed delete
	      OsiNode * save =newNode;
	      resolve();
	      if (isProvenOptimal()&&!isDualObjectiveLimitReached()) {
		newNode = new OsiNode(this,node,numberIntegers,which,lastws,
				      lowerBefore,upperBefore,solutionFound);
	      } else {
		newNode=NULL;
		break;
	      }
	      //save->nodeInfo_->decrement();
	      delete save;
	    }
	    if (newNode&&newNode->objectiveValue_<bestObjective) {
	      // push on stack
	      branchingTree.push(newNode);
	    } else {
	      delete newNode;
	    }
	    if (deleteNode)
	      delete node;
	  } else {
	    // maximum iterations - exit
	    delete node;
	    while (!branchingTree.empty()) {
	      // last node
	      OsiNode * node = branchingTree.top();
	      branchingTree.pop();
	      delete node;
	    }
	    handler_->message(OSI_BAB_MAXITS,messages_)
	      <<OsiMessageEol;
	    break;
	  }
	} else {
	  // integer solution - save
	  delete bestNode;
	  bestNode = new OsiFullNodeInfo(*this,numberIntegers,which,
					 node->nodeInfo());
	  // set cutoff (hard coded tolerance)
	  bestObjective = node->objectiveValue_-1.0e-4;
	  setDblParam(OsiDualObjectiveLimit,bestObjective);
	  handler_->message(OSI_BAB_SOLUTION,messages_)
	    <<node->objectiveValue_<<numberIterations
	    <<numberNodes
	    <<OsiMessageEol;
	  branchingTree.setComparison(compareExample);
	  branchingTree.setComparison(compareDepth);
	  solutionFound=1;
	}
	// take off rows 
	int numberNow = getNumRows();
	int numberToDelete = numberNow - numberRowsAtContinuous;
	if (numberToDelete) {
	  int * delRows = new int[numberToDelete];
	  for (i=0;i<numberToDelete;i++) {
	    delRows[i]=i+numberRowsAtContinuous;
	  }
	  deleteRows(numberToDelete,delRows);
	  delete [] delRows;
	}
      }
    }
    handler_->message(OSI_BAB_END,messages_)
	  <<numberIterations
	  <<numberNodes
	  <<OsiMessageEol;
    if (bestNode) {
      // we have a solution - restore
      bestNode->applyToModel(this,1,which,basis_);
      delete bestNode;
      resolve();
    }
    delete [] which;
    delete [] lowerBefore;
    delete [] upperBefore;
    delete [] walkback;
    delete [] numberRowsNow;
  } else {
    handler_->message(OSI_BAB_INFEAS,messages_)
    <<OsiMessageEol;
    throw CoinError("The LP relaxation is infeasible or too expensive",
		    "branchAndBound", "OsiClpSolverInterface");
  }
}
