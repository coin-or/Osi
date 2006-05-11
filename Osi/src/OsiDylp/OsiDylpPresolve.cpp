/*! \legal
  Copyright (C) 2002, 2003, 2004.
  Lou Hafer, International Business Machines Corporation and others.
  All Rights Reserved.

  This file is a portion of the COIN/OSI interface for dylp.
*/

#ifdef _MSC_VER

/* Turn off compiler warning about long names */

#  pragma warning(disable:4786)

#endif // _MSC_VER

/* Cut name lengths for readability. */

#define ODSI OsiDylpSolverInterface

/*
  \file OsiDylpPresolve.cpp

  \brief Implementation of `native' presolve for dylp using CoinPresolve.

  `Native' presolve for dylp. The basic outline, implemented within
  OsiDylpSolverInterface::initialSolve, is as follows:
  <ol>
    <li>
    Create an empty CoinPresolve object and load it from the ODSI object.
    Implemented by OsiDylpSolverInterface::initialisePresolve.

    <li>
    Apply the presolve transforms.
    Implemented by OsiDylpSolverInterface::doPresolve.

    <li>
    Consider whether it's worth proceeding --- if presolve had little or no
    effect, better to bail now and save the effort of loading the presolved
    system, postsolving, and reloading the original system.
    Implemented by OsiDylpSolverInterface::evalPresolve.

    <li>
    Save the current constraint system.
    Implemented by OsiDylpSolverInterface::saveOriginalSys.

    <li>
    Load the presolved problem.
    Implemented by OsiDylpSolverInterface::installPresolve

    <li>
    Run dylp (cold start).

    <li>
    Create a postsolve object from the presolve object. Destroy the presolve
    object in the process, after transferring its CoinPrePostsolveMatrix
    component to the CoinPostsolveMatrix object.
    Implemented by CoinPostsolveMatrix::assignPresolveToPostsolve

    <li>
    Load the postsolve object with the basis for the presolved problem.
    This is done by creating a warm start object and using it to set the status
    in the postsolve object.
    Implemented by CoinPrePostsolveMatrix::setStatus

    <li>
    Apply the postsolve transforms to convert the basis back to the
    original system. It's worth mentioning here that the only reason that
    postsolve transforms the constraint system is to keep the bookkeeping
    straight --- you'll notice there's no provision to extract the constraint
    system from the postsolve object. Postsolve will also handle things like
    primal and dual solution, but OsiDylp makes a call to dylp that will
    recreate all of this.
    Implemented by OsiDylpSolverInterface::doPostsolve.

    <li>
    Restore the saved constraint system. Make the postsolved basis the
    active basis. Destroy the postsolve object.
    Implemented by OsiDylpSolverInterface::installPostsolve.

    <li>
    Run dylp (warm start). This is (well, should be) a trivial (no-pivot)
    execution, simply to set up dylp's data structures. Given various bits
    of information loss and numerical twitchiness, a few pivots is reasonable.
  </ol>

  Memory consumption is a nontrivial issue here. Midway through step 5), we
  have the original consys, plus column- and row-major copies in preObj_,
  and we're building a new consys. Now I see why it's worthwhile to destroy
  chunks of the data in preObj_. Give this some attention.

  A thought concerning presolve for resolve(). Here, there's some utility in
  returning a basis that's suitably modified for the presolved problem, and
  we have an optimal basis to start from. But presumably presolve should not
  change the status of variables that remain.  Why not wait until presolve is
  completed, and modify the basis based on the information in
  originalColumns_ and originalRows_.
*/


//#define PRESOLVE_DEBUG 1
//#define PRESOLVE_CONSISTENCY 1

#include <OsiDylpSolverInterface.hpp>
#include <OsiDylpWarmStartBasis.hpp>
#include "OsiDylpMessages.hpp"
#include "CoinPresolveMatrix.hpp"

namespace {
  char sccsid[] UNUSED = "@(#)OsiDylpPresolve.cpp	1.6	11/06/04" ;
  char cvsid[] UNUSED = "$Id$" ;
}

#include "CoinPresolveEmpty.hpp"
#include "CoinPresolveFixed.hpp"
#include "CoinPresolveSingleton.hpp"
#include "CoinPresolveDoubleton.hpp"
#include "CoinPresolveTripleton.hpp"
#include "CoinPresolveZeros.hpp"
#include "CoinPresolveSubst.hpp"
#include "CoinPresolveForcing.hpp"
#include "CoinPresolveDual.hpp"
#include "CoinPresolveTighten.hpp"
#include "CoinPresolveUseless.hpp"
#include "CoinPresolveDupcol.hpp"
#include "CoinPresolveImpliedFree.hpp"
#include "CoinPresolveIsolated.hpp"

#if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
/*
  Debug stuff. This block covers the include and the namespace.
*/
#include "CoinPresolvePsdebug.hpp"

namespace {	// Anonymous namespace for debug routines

/*
  A control routine for debug checks --- keeps down the clutter in doPresolve.
  Each time it's called, it prints a list of transforms applied since the last
  call, then does the requested checks. Checks are controlled by two magic
  numbers:

  chkMtx: 1<<0	consistency of column- and row-major representations
	  1<<1	check column-major representation for duplicate entries
	  1<<2	check row-major representation for duplicate entries
	  1<<3	check column-major representation for explicit zero
	  1<<4	check row-major representation for explicit zero
	  1<<5	check column-major bulk storage management links
	  1<<6	check row-major bulk storage management links

  chkSol: 1<<0	check primal variables for NaN/Inf
	  1<<1	check primal variables for feasibility (variable bounds)
	  1<<2	check primal variable status vs. solution and bounds

	  1<<8	check constraint lhs for NaN/Inf
	  1<<9	check constraint lhs for accuracy, feasibility (row bounds)

  These will likely grow in the future. Mstrix checks work only if
  PRESOLVE_CONSISTENCY is defined. Solution checks work only if PRESOLVE_DEBUG
  is defined.

*/
void check_and_tell (ioid chn, const CoinPresolveMatrix *preObj_,
		     int chkMtx, int chkSol,
		     const CoinPresolveAction *first,
		     const CoinPresolveAction *&mark)

{ const CoinPresolveAction *current ;

  outfmt(chn,true,"PRESOLVE: ") ;
  if (first == 0 && mark == 0)
  { outfmt(chn,true,"checking prior to start.") ; }
  else
  if (first == mark)
  { outfmt(chn,true,"no transforms since last check.") ; }
  else
  { outfmt(chn,true,"applied") ;
    for (current = first ;
	 current != mark && current != 0 ;
	 current = current->next)
    { outfmt(chn,true," %s",current->name()) ; }
    mark = first ; }
  outfmt(chn,true,"\n") ;

# if !PRESOLVE_CONSISTENCY
  chkMtx = 0 ;
# endif
# if !PRESOLVE_DEBUG
  chkSol = 0 ;
# endif

  if (chkMtx&(1<<0))
  { presolve_consistent(preObj_) ; }

  if (chkMtx&((1<<1)|(1<<2)))
  { bool chkCol = chkMtx&((1<<1)) ;
    bool chkRow = chkMtx&((1<<2)) ;
    presolve_no_dups(preObj_,chkCol,chkRow) ; }

  if (chkMtx&((1<<3)|(1<<4)))
  { bool chkCol = chkMtx&((1<<3)) ;
    bool chkRow = chkMtx&((1<<4)) ;
    presolve_no_zeros(preObj_,chkCol,chkRow) ; }

  if (chkMtx&((1<<5)|(1<<6)))
  { bool chkCol = chkMtx&((1<<5)) ;
    bool chkRow = chkMtx&((1<<6)) ;
    presolve_links_ok(preObj_,chkCol,chkRow) ; }

/*
  As the routine is currently written, we need sol_ and colstat_ to do primal
  variable checks. In addition, we need acts_ to do constraint lhs checks.
  This could be more cleanly separated in presolve_check_sol with a bit of
  thought.
*/
  if (chkSol)
  { bool colSolPresent = (preObj_->sol_ && preObj_->colstat_) ;
    bool rowActPresent = (preObj_->acts_) ;
    int chkColSol = 0 ;
    int chkRowAct = 0 ;
    int chkStatus = 0 ;
    if (colSolPresent)
    { chkColSol = chkSol&((1<<1)|(1<<0)) ;
      chkStatus = (chkSol&(1<<2))>>2 ;
      if (rowActPresent)
      { chkRowAct = (chkSol&((1<<8)|(1<<9)))>>8 ; } }
    if (chkColSol+chkStatus+chkRowAct > 0)
    { presolve_check_sol(preObj_,chkColSol,chkRowAct,chkStatus) ; } }

  return ; }

} // end anonymous namespace for debug routines

#endif		// PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY

/*
  This routine creates and initializes the CoinPresolve object.

  keepIntegers: true to keep integer variable type information; false to
		treat all variables as continuous.
*/

CoinPresolveMatrix *ODSI::initialisePresolve (bool keepIntegers)

{ int m = getNumRows() ;
  int n = getNumCols() ;
  int numElems = getNumElements() ;

  CoinPresolveMatrix *preObj = new CoinPresolveMatrix(n,m,numElems) ;
  preObj->messageHandler()->setLogLevel(messageHandler()->logLevel()) ;

  const CoinPackedMatrix *pkMtx = getMatrixByCol() ;
  preObj->setMatrix(pkMtx) ;


  const double *dvec ;
  double val,exp ;

  dvec = getColLower() ;
  preObj->setColLower(dvec,-1) ;
  dvec = getColUpper() ;
  preObj->setColUpper(dvec,-1) ;

  dvec = getRowLower() ;
  preObj->setRowLower(dvec,-1) ;
  dvec = getRowUpper() ;
  preObj->setRowUpper(dvec,-1) ;

  dvec = getObjCoefficients() ;
  preObj->setCost(dvec,-1) ;
  val = getObjSense() ;
  preObj->setObjSense(val) ;
  getDblParam(OsiObjOffset,val) ;
  preObj->setObjOffset(val) ;

  getDblParam(OsiPrimalTolerance,val) ;
  exp = ((int) (.5 + log10((double) n))) - 2 ;
  if (exp > 0)
  { val *= pow(10,(double) exp) ; }
  preObj->setPrimalTolerance(val) ;
  preObj->setFeasibilityTolerance(1000*val) ;

  getDblParam(OsiDualTolerance,val) ;
  exp = ((int) (.5 + log10((double) m))) - 2 ;
  if (exp > 0)
  { val *= pow(10,(double) exp) ; }
  preObj->setDualTolerance(val) ;

  if (keepIntegers_)
  { unsigned char *variableType = new unsigned char [m] ;
    bool anyInteger = false ;
    for (int j = 0 ; j < m ; j++)
    { if (isInteger(j))
      { variableType[j] = 1 ;
	anyInteger = true ; }
      else
      { variableType[j] = 0 ; } }
    preObj->setVariableType(variableType,n) ;
    preObj->setAnyInteger(anyInteger) ;
    delete[] variableType ; }
  else
  { preObj->setVariableType(false,n) ;
    preObj->setAnyInteger(false) ; }

  preObj->setAnyProhibited(false) ;

  return (preObj) ; }



/*
  This routine orchestrates the application of presolve transforms to reduce
  the constraint system held in the presolve object (preObj_) to the
  presolved problem. Shamelessly adapted from OsiPresolve. Various comments
  questioning the wisdom of bits of code are scattered below.

  I've stripped the consistency checks so I can easily see structure while
  I'm composing the code. I will want to put them back in later. If I want
  to do that, I should consider moving them to CoinPresolveHelperFunctions.cpp
  or CoinPresolveDebug.cpp.
*/

void ODSI::doPresolve ()

{ 
/*
  Status codes used by CoinPresolveMatrix. These really should be symbolic.
  Fake it here.
*/
  const int feasibleStatus = 0 ;
  const int infeasibleStatus = 1 ;
  const int unboundedStatus = 2 ;
/*
  Some definitions to make it easy to control the presolve transforms that
  are applied. Unless you're debugging, all of these should be true.
*/
  const bool slackdoubleton = true ;
  const bool doubleton = true ;
  const bool tripleton = true ;
  const bool forcing = true ;
  const bool ifree = true ;
  const bool zerocost = true ;
  const bool dupcol = true ;
  const bool duprow = true ;
  const bool dual = true ;

  handler_->message(ODSI_PRESOL_STATS,messages_)
      << "Before presolve"
      << preObj_->getNumRows() << preObj_->getNumCols()
      << preObj_->getNumElems() << CoinMessageEol ;

/*
  We start with no postsolve transforms. Assume we're feasible.
*/
  postActions_ = 0 ;
  preObj_->setStatus(feasibleStatus) ;

# if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
/*
  chkMtx and chkSol are magic numbers to control checks. See the
  documentation with check_and_tell. chkMtx is nonfunctional unless
  PRESOLVE_CONSISTENCY is defined.  chkSol is nonfunctional unless
  PRESOLVE_DEBUG is defined.
*/
  const CoinPresolveAction *dbgActionMark = 0 ;
  int chkMtx = 0x7f ;
  int chkSol = 0x00 ;

  check_and_tell(local_logchn,
		 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
# endif

/*
  `Fix' variables before we get into the main transform loop. This transform
  physically removes from the problem variables that already have equal upper
  and lower bounds.
*/
  postActions_ = make_fixed(preObj_,postActions_) ;
# if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
  check_and_tell(local_logchn,
		 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
# endif
/*
  If we have integer variables, skip the presolve checks based on dual
  values.

  I'm unconvinced that integrality forces us to skip these, but that's
  probably because my understanding of the dual presolve checks is
  incomplete.
*/
  bool doDualStuff = (preObj_->anyInteger() == false && dual == true) ;
/*
  If we're feasible, set up for the main presolve transform loop by putting
  all (non-prohibited) rows and columns on the ToDo lists.

  Track the reduction in size of the constraint system (lastDropped). Also,
  CoinPresolveMatrix maintains a pass number which is used to control
  transform activity. Initialise it to 0.
*/
  if (preObj_->status() == feasibleStatus)
  { preObj_->initColsToDo() ;
    preObj_->initRowsToDo() ;
    int colCnt,rowCnt ;
    int colsDropped,rowsDropped ;
    int lastRowsDropped = 0 ;
    int lastColsDropped = 0 ;
    double rowShrink,colShrink ;
    preObj_->setPass(0) ;
/*
  Open the main presolve transform loop. There are two phases to each
  iteration:  a minor loop that does inexpensive presolve transforms until
  convergence, then a single application of the expensive transforms.
*/
    int currentPass ;
    for (currentPass = 0 ; currentPass < passLimit_ ; currentPass++)
    { const CoinPresolveAction *const lastAction = postActions_ ;

#     ifdef PRESOLVE_DEBUG
      outfmt(local_logchn,true,"Starting major pass %d\n",currentPass) ;
#     endif

/*
  implied_free_action looks for a constraint i that implicitly bounds a
  variable x<j> (think bounds tightening with arc consistency).  If column j
  is singleton (sole nonzero coefficient is a<ij>) then row i and column j
  can be removed from the problem without complications.  But if column j has
  two or more coefficients, complications ensue: we have to add multiples of
  row i to the other rows to eliminate x<j>.  fill_level puts an initial
  limit on the number of coefficients in a column. If there's exactly one
  other coefficient a<kj>, we can eliminate x<j> from row k without fillin.
  After that, fillin will, in general, occur.
*/
      int fill_level = 2 ;
/*
  Apply inexpensive transforms until convergence.

  In the original code, the criteria for `convergence' are
      (lastAction == lastSimpleAction && fill_level > 0)

  The first checks to see if we actually did anything (in which case, new
  postsolve transforms have been added).

  The second catches a hack that (in a nutshell) allows us to know when we've
  run low on doubleton columns for implied_free_action and its subordinate
  subst_constraint_action. A limit of 2 (and *only* 2) will be converted to
  -3 in subst_constraint_action. When -3 is passed in, it triggers i_f_a and
  s_c_a to consider all columns, not just those in ColsToDo. The -3 will come
  back as 3, and will not change again. (Except for this special case, a call
  to i_f_a will not change fill_level.)

  I'm asking why I can't simply look at numColsToDo and numRowsToDo. All the
  transforms queue to both NextColsToDo and NextRowsToDo (some indirectly due
  to triggered transforms). Except for the special case of implied_free,
  they all scan from ColsToDo or RowsToDo. Seems to me that if we break on
      (numColsToDo == 0 && numRowsToDo == 0 && fill_level > 0)
  we'll actually save one round of function calls.
*/
      while (1)
      { 
/*
  For some reason there's a hardwired limit on the number of
  slack_doubleton_actions performed per call.  notFinished comes back true if
  the number of opportunities exceeds the limit.
*/
	if (slackdoubleton == true)
	{ bool notFinished = true ;
	  while (notFinished) 
	  { postActions_ =
		slack_doubleton_action::presolve(preObj_,
						 postActions_,notFinished) ;
#	    if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	    check_and_tell(local_logchn,
			   preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	    endif
	  }
	  if (preObj_->status() != feasibleStatus)
	    break ; }
	if (doubleton == true)
	{ postActions_ = doubleton_action::presolve(preObj_,postActions_) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ; }
	if (tripleton == true)
	{ postActions_ = tripleton_action::presolve(preObj_,postActions_) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ; }
	if (zerocost == true)
	{ postActions_ = do_tighten_action::presolve(preObj_,postActions_) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ; }
	if (forcing == true)
	{ postActions_ =
	      forcing_constraint_action::presolve(preObj_,postActions_) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ; }
	if (ifree)
	{ postActions_ =
	      implied_free_action::presolve(preObj_,postActions_,fill_level) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ; }
/*
  Step preObj's ToDo arrays, moving queued rows and columns from NextToDo to
  ToDo and resetting changed flags.
*/
	preObj_->stepColsToDo() ;
	preObj_->stepRowsToDo() ;
/*
  Break out of the loop if no rows or columns are queued, unless we've
  triggered the special case for i_f_a.
*/
	if (preObj_->numberColsToDo() == 0 &&
		preObj_->numberRowsToDo() == 0 && fill_level > 0)
	  break ; }
/*
  That's it for the computationally cheap transforms. Reset the ToDo lists and
  try the more expensive transforms.
*/
      preObj_->initColsToDo() ;
      preObj_->initRowsToDo() ;
/*
  Attempt to fix variables at bound by propagating bounds on reduced costs.
  The routine does not consult the ToDo lists, so looking for the addition of
  postsolve actions is the only way to determine if another iteration might be
  useful. If we manage to fix variables, try to remove implied free variables
  but disable substitution of equalities into other constraints
  (fill_level = 0).
*/
      if (doDualStuff)
      { for (int iter = 0 ; iter < 5 ; iter++)
	{ const CoinPresolveAction *const marker = postActions_ ;
	  postActions_ = remove_dual_action::presolve(preObj_,postActions_) ;
#	  if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	  check_and_tell(local_logchn,
			 preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	  endif
	  if (preObj_->status() != feasibleStatus)
	    break ;
	  if (ifree)
	  { int fill_level = 0 ;
	    postActions_ =
		implied_free_action::presolve(preObj_,
					      postActions_,fill_level) ;
#	    if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	    check_and_tell(local_logchn,
			   preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	    endif
	    if (preObj_->status() != feasibleStatus)
	      break ; }
	  if (postActions_ == marker)
	    break ; } }
/*
  Check for duplicate columns and rows.
*/
      if (dupcol)
      { postActions_ = dupcol_action::presolve(preObj_,postActions_) ;
#	if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	check_and_tell(local_logchn,
		       preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	endif
	if (preObj_->status() != feasibleStatus)
	  break ; }
      if (duprow)
      { postActions_ = duprow_action::presolve(preObj_,postActions_) ;
#	if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
	check_and_tell(local_logchn,
		       preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#	endif
	if (preObj_->status() != feasibleStatus)
	  break ; }
/*
  See if we should continue. The requirement is that we've substantially
  reduced the problem. Checking postActions_ == lastAction is a quick test
  for any action at all. If there's action, see if it's enough to continue.
  For now, make the criteria a 5% reduction in rows or columns.
*/
      rowsDropped = preObj_->countEmptyRows() ;
      colsDropped = preObj_->countEmptyCols() ;
      colCnt = preObj_->getNumCols() ;
      rowCnt = preObj_->getNumRows() ;
      rowShrink = ((double) (rowsDropped-lastRowsDropped))/rowCnt ;
      colShrink = ((double) (colsDropped-lastColsDropped))/colCnt ;

      handler_->message(ODSI_PRESOL_PASS,messages_)
	  << currentPass << rowsDropped-lastRowsDropped << rowShrink
	  << colsDropped-lastColsDropped << colShrink << CoinMessageEol ;

      if (postActions_ == lastAction) break ;
      if (rowShrink < .05 && colShrink < .05) break ;

      lastRowsDropped = rowsDropped ;
      lastColsDropped = colsDropped ; } }
/*
  That's it, we're done with the presolve major loop. If we're still
  feasible, it's time for final cleanup: drop zero coefficients from the
  matrix, then drop empty rows and columns.
*/
  if (preObj_->status() == feasibleStatus)
  { postActions_ = drop_zero_coefficients(preObj_,postActions_) ;
#   if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
    check_and_tell(local_logchn,
	           preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#   endif
    postActions_ = drop_empty_cols_action::presolve(preObj_,postActions_) ;
#   if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
    chkMtx &= 0x2a ;
    check_and_tell(local_logchn,
		   preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#   endif
    postActions_ = drop_empty_rows_action::presolve(preObj_,postActions_) ;
#   if PRESOLVE_DEBUG || PRESOLVE_CONSISTENCY
    check_and_tell(local_logchn,
		   preObj_,chkMtx,chkSol,postActions_,dbgActionMark) ;
#   endif
    handler_->message(ODSI_PRESOL_STATS,messages_)
      << "After presolve"
      << preObj_->getNumRows() << preObj_->getNumCols()
      << preObj_->getNumElems() << CoinMessageEol ;
    }

/*
  Whoops --- we've come up infeasible or unbounded. Let the user know.
*/
  else
  { CoinMessageHandler *handler = preObj_->messageHandler() ;
    if (preObj_->status() == infeasibleStatus)
    { handler->message(COIN_PRESOLVE_INFEAS,preObj_->messages())
	  << preObj_->feasibilityTolerance_ << CoinMessageEol ; }
    else
    if (preObj_->status() == unboundedStatus)
    { handler->message(COIN_PRESOLVE_UNBOUND,preObj_->messages())
	  << CoinMessageEol ; }
    else
    { handler->message(COIN_PRESOLVE_INFEASUNBOUND,preObj_->messages())
	  << CoinMessageEol ; } }

  return ; }



/*
  This routine has the responsibility of deciding whether presolve has done
  enough to make it worth installing. Or too much.

  For now, take the attitude that if we don't manage to eliminate any
  constraints, then presolve likely wasn't worth the effort. If we come
  back infeasible or unbounded, well, better let dylp check it on its own.

  At the other end of the spectrum, a B&C code will often pass in a problem in
  which a majority (or all) of the variables are fixed. Presolve may well
  reduce the problem to 0x0. In this case, there's literally nothing to do.
*/

bool ODSI::evalPresolve ()

{ int origCons = getNumRows() ;
  int presolCons = preObj_->getNumRows() ;
  int presolVars = preObj_->getNumCols() ;
  int presolStatus = preObj_->status() ;

/*
  Presolve returns infeasible or unbounded.
*/
  if (presolStatus != 0)
    return (false) ;
/*
  Presolve didn't manage to reduce the number of constraints.
*/
  if (origCons <= presolCons)
    return (false) ;
/*
  Presolve worked all too well; the constraint system is 0x0. I don't think
  it's possible for presolve to return 0xn (it should fix otherwise
  unconstrained bounded variables at some bound or declare unboundedness).
  But it never hurts to check.
*/
  if (presolCons == 0 && presolVars == 0)
    return (false) ;

  return (true) ; }



/*
  This routine saves the original constraint system. Trivial, but broken out
  in anticipation more functionality might be useful in the future. From an
  efficiency perspective, this is a wonderful idea. From a space perspective,
  it's not so wonderful.
*/

void ODSI::saveOriginalSys ()

{
/*
  If we don't have a constraint system, we're deeply confused. Also, the value
  of lpprob->consys (if lpprob exists) must equal the value of consys at
  all times.
*/
  assert(consys) ;

  savedConsys_ = consys ;
  consys = 0 ;
  if (lpprob) lpprob->consys = 0 ;
/*
  Presolve should be invisible to clients, and one reasonable expectation is
  that pointers to structural vectors obtained from ODSI will remain constant
  across calls to the solver. However, the presolved problem is structurally
  different, so we need to hide these while dylp works on the presolved
  problem.
*/
  saved_col_obj = _col_obj ; _col_obj = 0 ;
  saved_row_rhs = _row_rhs ; _row_rhs = 0 ;
  saved_row_lower = _row_lower ; _row_lower = 0 ;
  saved_row_upper = _row_upper ; _row_upper = 0 ;
  saved_row_sense = _row_sense ; _row_sense = 0 ;
  saved_row_range = _row_range ; _row_range = 0 ;
  saved_matrix_by_col = _matrix_by_col ; _matrix_by_col = 0 ;
  saved_matrix_by_row = _matrix_by_row ; _matrix_by_row = 0 ;

  return ; }



/*
  This routine has the responsibility of installing a presolved problem in
  the ODSI object. Simply a matter of retrieving parameters and calling the
  right version of load_problem. We also need to replace the lpprob_struct.
*/


void ODSI::installPresolve ()

{ int n = preObj_->getNumCols() ;
  int m = preObj_->getNumRows() ;

  const CoinBigIndex *colStarts = preObj_->getColStarts() ;
  const int *colLens = preObj_->getColLengths() ;
  const int *rowIndices = preObj_->getRowIndicesByCol() ;
  const double *coeffs = preObj_->getElementsByCol() ;

  const double *vlbs = preObj_->getColLower() ;
  const double *vubs = preObj_->getColUpper() ;
  const double *obj = preObj_->getCost() ;

  const double *clbs = preObj_->getRowLower() ;
  const double *cubs = preObj_->getRowUpper() ;

  double *rhs = new double[m] ;
  double *rhslow = new double[m] ;
  contyp_enum *ctyp = new contyp_enum[m] ;

  gen_rowparms(m,rhs,rhslow,ctyp,clbs,cubs) ;

  load_problem(n,m,colStarts,colLens,rowIndices,coeffs,
	       vlbs,vubs,obj,ctyp,rhs,rhslow) ;
  construct_lpprob() ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ;

  return ; }






/*
  This routine is responsible for setting up for postsolve. It creates 
  the postsolve object and does the major part of the initialisation by
  cannibalising the presolve object. It then loads solution information for the
  presolved problem: status, primal and dual solutions, and primal and dual
  activity.

  The call to assignPresolvetoPostsolve will destroy preObj and set the pointer
  to null.
*/

CoinPostsolveMatrix *ODSI::initialisePostsolve (CoinPresolveMatrix *&preObj)

{ CoinPostsolveMatrix *postObj = new CoinPostsolveMatrix(0,0,0) ;
  postObj->assignPresolveToPostsolve(preObj) ;

  CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis *>(getWarmStart()) ;
  postObj->setStatus(ws) ;
  delete ws ;

  const double *temp ;
  temp = getColSolution() ;
  postObj->setColSolution(temp,-1) ;
  temp = getRowActivity() ;
  postObj->setRowActivity(temp,-1) ;
  temp = getRowPrice() ;
  postObj->setRowPrice(temp,-1) ;
  temp = getReducedCost() ;
  postObj->setReducedCost(temp,-1) ;

  return (postObj) ; }


/*
  This routine is responsible for applying the list of postsolve actions
  in order to restore the original system. After each transform is applied,
  it's destroyed.
*/

void ODSI::doPostsolve ()

{ handler_->message(ODSI_POSTSOL,messages_) << "start" << CoinMessageEol ;
# if PRESOLVE_CONSISTENCY
  presolve_check_threads(postObj_) ;
  presolve_check_free_list(postObj_) ;
# endif
# if PRESOLVE_DEBUG
  presolve_check_sol(postObj_) ;
  presolve_check_reduced_costs(postObj_) ;
  presolve_check_duals(postObj_) ;
# endif
  while (postActions_ != 0)
  { const CoinPresolveAction *postAction = postActions_ ;
    postActions_ = postAction->next ;
    handler_->message(ODSI_POSTSOL_ACT,messages_)
      << postAction->name() << CoinMessageEol ;
    postAction->postsolve(postObj_) ;
#   if PRESOLVE_CONSISTENCY
    presolve_check_threads(postObj_) ;
    presolve_check_free_list(postObj_) ;
#   endif
#   if PRESOLVE_DEBUG
    presolve_check_sol(postObj_) ;
    presolve_check_reduced_costs(postObj_) ;
    presolve_check_duals(postObj_) ;
#   endif
    delete postAction ; }
  handler_->message(ODSI_POSTSOL,messages_) << "complete" << CoinMessageEol ;
  
  return ; }



/*
  This routine has the responsibility of restoring the original constraint
  system, setting the corresponding basis as the new active basis for
  this ODSI object.

  There's a really annoying problem here --- CoinWarmStartBasis does not
  support superbasic status, and postsolve will return superbasic status.
*/

void ODSI::installPostsolve ()

{ 
/*
  Restore the constraint system --- trivially easy, just destroy the presolved
  system and replace it with the saved copy of the original. We need to claim
  ownership of the solver.
*/
  if (flgon(lpprob->ctlopts,lpctlDYVALID))
    dylp_owner = this ;
  destruct_problem(true) ;
  consys = savedConsys_ ;
  savedConsys_ = 0 ;
/*
  Restore cached structural vectors.
*/
  _col_obj = saved_col_obj ; saved_col_obj = 0 ;
  _row_rhs = saved_row_rhs ; saved_row_rhs = 0 ;
  _row_lower = saved_row_lower ; saved_row_lower = 0 ;
  _row_upper = saved_row_upper ; saved_row_upper = 0 ;
  _row_sense = saved_row_sense ; saved_row_sense = 0 ;
  _row_range = saved_row_range ; saved_row_range = 0 ;
  _matrix_by_col = saved_matrix_by_col ; saved_matrix_by_col = 0 ;
  _matrix_by_row = saved_matrix_by_row ; saved_matrix_by_row = 0 ;
/*
  Set the basis.
*/
  CoinWarmStart *postBasis =
      dynamic_cast<CoinWarmStart *>(postObj_->getStatus()) ;
  if (setWarmStart(postBasis) == false)
  { throw CoinError("Could not install postsolve basis.",
		    "installPostsolve","OsiDylpSolverInterface") ; }
  delete postBasis ;
/*
  And destroy the postsolve object.
*/
  delete postObj_ ;
  postObj_ = 0 ;

  return ; }



/*
  This routine deletes data structures associated with presolve. It's a utility
  cleanup routine.
*/

void ODSI::destruct_presolve ()

{ if (preObj_)
  { delete preObj_ ;
    preObj_ = 0 ; }
  if (postObj_)
  { delete postObj_ ;
    postObj_ = 0 ; }

  while (postActions_ != 0)
  { const CoinPresolveAction *postAction = postActions_ ;
    postActions_ = postAction->next ;
    delete postAction ; }
/*
  Presence of savedConsys_ means that we called saveOriginalSys.
*/
  if (savedConsys_)
  { consys_free(savedConsys_) ; savedConsys_ = 0 ;
    _col_obj = saved_col_obj ; saved_col_obj = 0 ;
    _row_rhs = saved_row_rhs ; saved_row_rhs = 0 ;
    _row_lower = saved_row_lower ; saved_row_lower = 0 ;
    _row_upper = saved_row_upper ; saved_row_upper = 0 ;
    _row_sense = saved_row_sense ; saved_row_sense = 0 ;
    _row_range = saved_row_range ; saved_row_range = 0 ;
    _matrix_by_col = saved_matrix_by_col ; saved_matrix_by_col = 0 ;
    _matrix_by_row = saved_matrix_by_row ; saved_matrix_by_row = 0 ; }

  return ; }
