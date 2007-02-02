/*! \legal
  Copyright (C) 2002, 2003, 2004.
  Lou Hafer, Stephen Tse, International Business Machines Corporation and
  others. All Rights Reserved.

  This file is a portion of the COIN/OSI interface for dylp.
*/

#ifdef _MSC_VER

/* Turn off compiler warning about long names */
# pragma warning(disable:4786)

/*
  MS C++ doesn't want to cope with this set of templates. Since they're just
  paranoid checks, it's easiest to simply disable them.
*/

# define assert_same

/*
  In general, templates do not seem to be a strong point for MS C++. To get
  around this, there are a number of macros defined below.  For Sun Workshop
  and Gnu programming environments, these will expand to template functions.
  For MS, they expand to some code plus calls to STL functions. Apparently MS
  C++ doesn't have trouble with the definitions of the templates, just
  instantiations.
*/

#endif	// _MSC_VER

/*
  Right now, there's no way to ask (or set) the value that COIN uses for
  infinity. Portions of COIN code depend on having a finite infinity ---
  specifically, DBL_MAX --- and this is hardwired into the code. This
  definition tries to anticipate the future, when one can inquire.
  CoinFinite.hpp defines COIN_DBL_MAX, CoinFinite, and CoinIsnan, hiding
  platform idiosyncracies.

  Dylp can work with both finite infinity (DBL_MAX) and infinite infinity (IEEE
  infinity).
  
  The value of ODSI::odsiInfinity is initialised to CoinInfinity, and dylp
  is told to use odsiInfinity as its internal value for infinity.
*/

#include "CoinFinite.hpp"
const double CoinInfinity = COIN_DBL_MAX ;

#include "CoinTypes.hpp"

/* Cut name lengths for readability. */

#define ODSI OsiDylpSolverInterface
#define OSI OsiSolverInterface

#include <string>
#include <cassert>
#include <sstream>
#include <sstream>
#include "CoinTime.hpp"
#include <OsiColCut.hpp>
#include <OsiRowCut.hpp>
#include <OsiRowCutDebugger.hpp>
#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpMessages.hpp"
#include "OsiDylpWarmStartBasis.hpp"
#include "OsiPresolve.hpp"

namespace {
  char sccsid[] UNUSED = "@(#)OsiDylpSolverInterface.cpp	1.20	11/13/04" ;
  char cvsid[] UNUSED = "$Id$" ;
}

/*! \brief Define to enable implicit read of options file

  If ODSI_IMPLICIT_SPC is defined `1' at compile time, when
  readMps(const char*,const char*) is asked to read <problem>.mps, it will
  first look for <problem>.spc and, if found, process it as a dylp options
  file.

  Disabled by default. Guaranteed not to work in Windows environments.
*/

#define ODSI_IMPLICIT_SPC 0

/*! \brief Define to control behaviour when a request is made to access a
	   stale solution.

  A solution is `fresh' if the problem has not been modified since the most
  recent call to dylp. The OSI specification says that problem modifications
  should invalidate the current solution. In practice this nicety is
  sometimes stretched a bit. Defining ODSI_STRICTLY_FRESH will cause ODSI to
  throw an exception if asked for stale solution data. Otherise, it will
  ignore the problem and return stale data. At any log level greater than 0,
  you'll get a warning.
*/

#define ODSI_STRICTLY_FRESH 1

/*
  The following symbols are useful only for detailed debugging.

  ODSI_TRACK_FRESH	track how a solution is made stale/fresh
  ODSI_TRACK_SOLVERS	track creation, use, and deletion of ODSI objects
  ODSI_TRACK_ACTIVE	track creation and deletion of activeBasis (the
			active warm start object)

 #define ODSI_TRACK_FRESH 1
 #define ODSI_TRACK_SOLVERS 1
 #define ODSI_TRACK_ACTIVE 1
*/

/*! \brief Define to enable paranoid checks.

  When non-zero, enables various paranoid checks. 
    1: bounds checks on indices
    2: test for equality after structure copy operations and test for
       basis internal consistency.

  An error will cause a throw. Configuration should set this symbol to 0 for
  an optimised build, 2 if --enable-osidylp-paranoia is requested.

  In particular, this symbol must be defined as >= 1 in order for OsiCbc(dylp)
  to pass the OsiCbc unit test.
*/

#ifndef ODSI_PARANOIA
# define ODSI_PARANOIA 1
#endif


/*! \brief Define to enable statistics collection in dylp

  When defined, ODSI will pass an lpstats_struct to dylp to be filled with
  statistics. Information is collected on a per-call basis only.  In order to
  do so, OsiDylpSolverInterface must be built with the compile-time symbol \c
  ODSI_STATISTICS defined and the dylp library must be built with \c
  DYLP_STATISTICS defined. Normally this is handled through options to the
  configure script.

  #define ODSI_STATISTICS
*/

/*! \brief Define to produce informational messages

  When defined, informational printing is enabled. Nothing will be printed,
  however, unless the log level is sufficiently high. Normally handled through
  options to the configure script.

  #define ODSI_INFOMSGS
*/

/*!
  \file OsiDylpSolverInterface.cpp

  \brief Implementation of COIN OSI layer for dylp.

  This file contains the implementation of a COIN OSI layer for dylp, an lp
  solver originally written for the bonsaiG MILP code.

  More information on the COIN/OR project and the OSI layer specification
  can be found at http://www.coin-or.org.

  <h3>Interface Principles for Implementors</h3>

  In addition to the principles described in the documentation for the
  OsiDylpSolverInterface object in OsiDylpSolverInterface.hpp, implementors
  should be aware of the following:

  <strong>Problem Structures</strong>:
  Within dylp, a constraint system is held in a row- and column-linked
  structure called a \c consys_struct. The lp problem (constraint system plus
  additional problem information --- primal and dual variables, basis, status,
  <i>etc.</i>) is held in a structure called an \c lpprob_struct.
  
  Down in the private section of the class, OsiDylpSolverInterface::consys is
  the constraint system, and OsiDylpSolverInterface::lpprob is the
  \c lpprob_struct that's used to pass the problem to \c dylp and
  return the results.
  The structures OsiDylpSolverInterface::initialSolveOptions,
  OsiDylpSolverInterface::resolveOptions, OsiDylpSolverInterface::tolerances,
  and OsiDylpSolverInterface::statistics respectively hold control parameters
  tolerances, and statistics collected as the problem is solved.

  <strong>Options and Tolerances</strong>:
  The \c initialSolveOptions, \c resolveOptions, and \c tolerances structures
  are created and initialised by the ODSI constructor, so they're valid from
  the moment the ODSI object is created. Option and tolerance settings are
  preserved when a new problem is loaded or assigned.  The
  OsiDylpSolverInterface::reset() method will reset them to their defaults.

  When a call is made to the solver, a local copy of the options and
  tolerances is made, then tweaked with a call to \c dy_checkdefaults. This
  partially subverts attempts to set tolerances using
  OsiDylpSolverInterface::dylp_controlfile() or any of the various ODSI
  parameter set methods.

  <strong>Statistics Collection</strong>:
  Dylp can collect detailed statistics on its performance. In order to do so,
  OsiDylpSolverInterface must be built with the compile-time symbol
  \c ODSI_STATISTICS defined and the dylp library must be built with
  \c DYLP_STATISTICS defined. Normally this is handled through options to
  the configure script.

  Statistics are collected on a per-call basis only; there is no facility to
  accumulate statistics over multiple calls to dylp. The statistics structure
  must also be sized to the problem. For these reasons, its creation is
  delayed until the actual call to the solver. Any existing statistics
  structure is deleted before the new one is created. Statistics structures are
  not copied at assignment or cloning.

  <strong>Constraint System Existence</strong>:
  An invariant in the interface is that if
  OsiDylpSolverInterface::getNumCols() or
  OsiDylpSolverInterface::getNumRows() returns a value greater than 0, then a
  valid constraint system structure (consys) exists, and attached to it are
  valid vectors for the objective coefficients (obj), variable type (vtyp)
  and upper and lower bounds (vub, vlb), and constraint type (ctyp) and
  right-hand side (rhs, rhslow).

  <strong>LP Problem Existence</strong>:
  The LP problem structure is not created until it's needed. A valid LP
  problem structure will exist after the first call to the solver
  (OsiDylpSolverInterface::initialSolve) or after a warm start object
  is loaded (OsiDylpSolverInterface::setWarmStart).

  <strong>Caching</strong>:
  Since vectors returned by OsiSolverInterface "get" functions are constant
  (not modifiable by users) and it is convenient for users to make repeated
  calls for the same information, a cache mechanism is used when necessary
  to avoid repeated expensive conversions in the interface. All cache
  variables are prefixed with underscore (_col*, _row*, _matrix*).  In
  general, modifying the problem causes the cache to be invalidated.

  Most cached values are <i>not</i> replicated when an ODSI is cloned or
  assigned. The exceptions are #_col_x, #_row_price, and #_objval, which must
  be copied because the functions #setColSolution() and #setRowPrice() use
  them directly to hold the values specified by the client.

  <strong>`Fresh' Solutions</strong>

  The OSI specification says that changes to the problem structure should
  invalidate the current solution. However, the de facto standards (the
  OSI unitTest and the OsiClp implementation) are more relevant, particularly
  for solvers that expect to work in Cbc. Cbc pushes the rules pretty hard.
  OsiDylp offers a compile-time option, ODSI_STRICTLY_FRESH. When this symbol
  is defined, ODSI will throw an exception if the user tries to use a stale
  solution. If it's not defined, ODSI will return whatever information
  it has. In any event, you'll get a warning at any log level greater than 0.
  ODSI_STRICTLY_FRESH is enabled by default.

  The actual behaviour is a bit more complicated, due to the extreme cases.
  The OSI unitTest requires that an empty solver (no constraint system) return
  null when asked for solution vectors. Also, the client can load in a primal
  or dual solution. ODSI keeps these in the cached solution vectors (and they
  are totally irrelevant to dylp). Similarly, when a problem is loaded, a
  pessimal primal solution is generated and cached. The actual tests are:
  <ul>
    <li>
    If a cached solution vector exists, return it. This implies that ODSI
    must be scrupulous about invalidating the cache when changes are made
    to the problem. The converse is also true; the cache should never be
    invalidated unless there's been a change in the problem.
    </li>
    <li>
    If the solver is empty (no consys), return null.
    </li>
    <li>
    If #solnIsFresh is true, rebuild cached solution data as needed from
    lpprob and return it. If #solnIsFresh is false, act according to
    ODSI_STRICTLY_FRESH, as described above.
    </li>
  </ul>

  <strong>Index Base</strong>:
  OsiSolverInterface indexes variables and constraints from 0 (the standard
  approach for C/C++) while dylp indexes them from 1 (for robustness; see
  consys.h). Some caution is needed in the interface to avoid off-by-one
  errors. Similarly, arrays (vectors) for k elements in dylp occupy k+1
  space with vector[0] unused.  ODSI::idx, ODSI::inv, and ODSI::inv_vec are
  used to make these trivial conversions explicit.

  <strong>Copy and Assert</strong>:
  To clone, each part of dylp's data structures (lpprob, consys, basis) is
  copied (ODSI::copy*) and then (optionally) verified (ODSI::assert_same*
  and assert). Helpers for copying and verifying primitive types (integer
  and double) and array types are implemented with C++ templates. (If you
  skipped it, this is a good moment to check the source file and read the
  note at the beginning about MS C++ and templates.)

  <strong>Constraint Grooming</strong>
  The OSI layer (or at least the OSI test suite) expects that the constraint
  system read back from the solver will be the same as the constraint system
  given to the solver.  Dylp expects that any grooming of the constraint
  system will be done before it's called, and furthermore expects that this
  grooming will include conversion of >= constraints to <= constraints.  The
  solution here is to do this grooming in a wrapper, do_lp.
  do_lp only implements the conversion from >= to <=, which is
  reversible. dylp can tolerate empty constraints.  do_lp also embodies a
  rudimentary strategy that attempts to recover from numerical inaccuracy by
  refactoring the basis more often (on the premise that this will reduce
  numerical inaccuracy in calculations involving the basis inverse).

*/

namespace {
/*
  A little print helper routine for the ODSI_start_enum type.
*/

const char *startString (ODSI_start_enum start)

{ switch (start)
  { case startCold:
    { return ("cold") ; }
    case startWarm:
    { return ("warm") ; }
    case startHot:
    { return ("hot") ; }
    default:
    { return ("!invalid!") ; } } }

}	// end unnamed file-local namespace

using std::string ;
using std::vector ;


extern "C"
{

#include "dy_cmdint.h"

#ifndef DYLP_ERRMSGDIR
/*
  This is the correct path to find dy_errmsgs.txt, assuming the default
  COIN directory structure for Osi and Dylp, to wit:
    COIN
      Osi
	src
	  OsiDylp
      DyLP
	src
	  Dylp
*/
# define DYLP_ERRMSGDIR "../../../DyLP/src/Dylp"
#endif

extern void dy_initbasis(int concnt, int factor_freq, double zero_tol) ;
extern void dy_freebasis() ;

/*
  These routines work with the literal tree in DylpLib (a legacy, what can I
  say). DO NOT MIX strings handled with stralloc/strfree with strings handled
  using malloc/free.  Needless to say, std::string is another beast
  entirely.
*/

extern char *stralloc(const char *str) ;
extern void strfree(const char *str) ;

}

/*!
  \defgroup DylpIO dylp i/o control variables
  \brief Variables controlling dylp file and terminal i/o.

  Dylp is capable of generating a great deal of output, but the control
  mechanisms are not a good fit for C++ objects and the OSI framework.  Four
  global variables, dy_cmdchn, dy_cmdecho, dy_logchn, and dy_gtxecho, control
  command input/echoing and log message output/echoing.

  dy_gtxecho can be controlled using the OsiDoReducePrint hint; it
  is set to true whenever the print level is specified as an absolute integer
  such that (print level) & 0x10 != 0

  By default, solver command (option) files are handled in the following
  manner: When an MPS file example.mps is read, the solver looks for an
  options file example.spc in the same directory. If it exists, it is
  processed. This is a local file open/close, well-defined and clean, but
  with little user control. The commands will not be echoed.

  The user can process and arbitrary command file by calling dylp_controlfile
  (a non-OSI function).

  Log output is disabled by default. The user can enable logging by calling
  dylp_logfile (a non-OSI function).

  For historical reasons, dylp uses a private i/o subsystem for file and
  terminal i/o, built on top of the C stdio library.
*/

//@{
/*! \var ioid dy_cmdchn
    \brief ioid used to read option (.spc) files
*/
/*! \var ioid dy_logchn
    \brief ioid used for logging to a file
*/
/*! \var bool dy_cmdecho
    \brief controls echoing of commands from option files to the terminal
*/
/*! \var bool dy_gtxecho
    \brief controls echoing of generated text to the terminal
*/

ioid dy_cmdchn = IOID_NOSTRM,
     dy_logchn = IOID_NOSTRM ;

bool dy_cmdecho = false,
     dy_gtxecho = false ;

//@}

/*!
  \defgroup DylpResidual dylp residual control variables
  \brief Variables controlling persistent dylp subsystems

  These variables are used to control initialisation and shutdown of the i/o
  and basis maintenance subsystems, and to track the state of the dylp solver.
  
  The i/o subsystem is initialised when the first ODSI object is created
  (constructors, reference_count == 1), and shut down when the last ODSI
  object is destroyed (destructor, reference_count == 0).

  The basis maintenance package is initialised at the first attempt to solve
  an lp (initialSolve, basis_ready == false, or resolve, basis_ready ==
  false) and shut down when the last ODSI object is destroyed (destructor,
  reference_count == 0, basis_ready == true).
*/

//@{
/*! \var int OsiDylpSolverInterface::reference_count
    \brief Tracks the number of ODSI objects in existence.
*/
/*! \var bool OsiDylpSolverInterface::basis_ready
    \brief Tracks basis initialisation

    Initialisation of the basis data structures is delayed until the user
    actually attempts to solve an lp.
*/
/*! \var OsiDylpSolverInterface *OsiDylpSolverInterface::dylp_owner
    \brief ODSI instance controlling dylp

    Records the ODSI instance which owns retained state in dylp. Used to
    relinquish ownership at destruction or when another ODSI instance wants
    to use the solver.
*/

int ODSI::reference_count = 0 ;
bool ODSI::basis_ready = false ;
ODSI *ODSI::dylp_owner = 0 ;

//@}



/*! \defgroup VectorHelpers Vector Helper Functions
    \brief Vector helper functions

  Dylp uses 1-based indexing rather than the usual 0-based indexing; idx,
  idx_vec, inv, and inv_vec make this conversion explicit.
  See the \link OsiDylpSolverInterface.cpp file comments \endlink for a brief
  description.
  This group also provides a function for conversion between the OSI notion of
  a packed vector and the dylp notion of a packed vector.
*/
//@{

/*! \brief Convert 0-based index to 1-based index */

inline int ODSI::idx (int i) { return i+1 ; }

/*! \brief Convert 0-based vector pointer to 1-based vector pointer

  For cases where it's inconvenient to adjust indices, the alternative is to
  adjust the pointer to the vector so it points to vector[-1]. Be careful!
*/

template<class T> inline T* ODSI::idx_vec (T* vec) { return vec-1 ; }

/*! \brief Convert 1-based index to 0-based index */

inline int ODSI::inv (int i) { return i-1 ; }

/*! \brief Convert 1-based vector pointer to 0-based vector pointer

  For cases where it's inconvenient to adjust indices, the alternative is to
  adjust the pointer to the vector so it points to vector[1].
*/

template<class T> inline T* ODSI::inv_vec (T* vec) { return vec+1 ; }

/*
  MS C++ doesn't like the template functions above. To keep it happy, wrap
  them in macros which will accomplish the same thing.
*/

#ifdef _MSC_VER
#  define IDX_VEC(zz_type_zz,zz_vec_zz) (zz_vec_zz-1)
#  define INV_VEC(zz_type_zz,zz_vec_zz) (zz_vec_zz+1)
#else
#  define IDX_VEC(zz_type_zz,zz_vec_zz) idx_vec<zz_type_zz>(zz_vec_zz)
#  define INV_VEC(zz_type_zz,zz_vec_zz) inv_vec<zz_type_zz>(zz_vec_zz)
#endif

/*! \brief Load a CoinShallowPackedVector into a new dylp pkvec.

  pkvec's created with this routine must eventually be freed with pkvec_free.
  Note that 0-based indexing is used for the entries of the coeff vector
  (i.e., coeffs[0] is valid).

  \param dimension The length of the vector if it were to be unpacked.
*/

pkvec_struct* ODSI::packed_vector (const CoinShallowPackedVector src,
				   int dimension)

{ int n = src.getNumElements() ;
  pkvec_struct* dst = pkvec_new(n) ;

  assert(dst) ;

  if (n == 0) return (dst) ;

  packed_vector(src,dimension,dst) ;

  return (dst) ; }


/*! \brief Load a CoinShallowPackedVector to an existing dylp pkvec.

  Note that 0-based indexing is used for the entries of the coeff vector
  (i.e., coeffs[0] is valid).

  \param dimension The length of the vector if it were to be unpacked.
*/

void ODSI::packed_vector (const CoinShallowPackedVector src, int dimension,
			  pkvec_struct *dst)

{ int n = src.getNumElements() ;

  dst->cnt = n ;
  dst->dim = dimension ;

  if (n == 0) return ;

  const int* indices = src.getIndices() ;
  const double* elements = src.getElements() ;
  pkcoeff_struct* coeffs = dst->coeffs ;

  for (int i = 0 ; i < n ; i++)
  { coeffs[i].ndx = idx(indices[i]) ;
    coeffs[i].val = elements[i] ; }

  return ; }


//@} // VectorHelpers



/*! \defgroup CopyHelpers Copy Helpers

  Copy template functions for simple vectors and fixed size objects, and
  specializations for various complex structures.
*/
//@{

/* Templates for simple copies */

/*! \brief Copy a structure byte-wise.

  Destination structure is allocated and returned to caller.
*/

template<class T> inline T* ODSI::copy (const T* src)

{ if (!src) return 0 ;
  T* dst = new T ;
  memcpy(dst, src, sizeof(T)) ;
  return dst ;
}

#ifdef _MSC_VER
#  define CLONE(zz_type_zz,zz_src_zz,zz_dst_zz) \
	    if (zz_src_zz) \
	    { zz_dst_zz = new zz_type_zz ; \
	      memcpy(zz_dst_zz,zz_src_zz,sizeof(zz_type_zz)) ; } \
	    else \
	    { zz_dst_zz = NULL ; }
#else
#  define CLONE(zz_type_zz,zz_src_zz,zz_dst_zz) \
	    zz_dst_zz = copy<zz_type_zz>(zz_src_zz) ;
#endif
  

/*! \brief Copy a primitive array of values element-wise.

  Caller supplies the destination vector.
*/

template<class T> inline void ODSI::copy (const T* src, T* dst, int n)

{ if (!dst || !src || n == 0) return ;
  int size = sizeof(T) * n ;
  memcpy(dst,src,size) ;
}

#ifdef _MSC_VER
# define COPY_VEC(zz_type_zz,zz_src_zz,zz_dst_zz,zz_n_zz) \
  	   if (zz_dst_zz && zz_src_zz && zz_n_zz > 0 ) \
	     std::copy(zz_src_zz,zz_src_zz+zz_n_zz,zz_dst_zz)
#else
# define COPY_VEC(zz_type_zz,zz_src_zz,zz_dst_zz,zz_n_zz) \
	   copy<zz_type_zz>(zz_src_zz,zz_dst_zz,zz_n_zz)
#endif


/*! \brief Copy a primitive array of values element-wise.

  Destination vector is allocated and returned to caller.
*/

template<class T> inline T* ODSI::copy (const T* src, int n)

{ if (!src || n == 0) return 0 ;
  T* dst = new T[n] ;
  copy(src, dst, n) ;
  return dst ;
}

#ifdef _MSC_VER
#  define CLONE_VEC(zz_type_zz,zz_src_zz,zz_dst_zz,zz_n_zz) \
	    if (zz_src_zz && zz_n_zz > 0) \
	    { zz_dst_zz = new zz_type_zz[zz_n_zz] ; \
	      std::copy(zz_src_zz,zz_src_zz+zz_n_zz,zz_dst_zz) ; } \
	    else \
	    { zz_dst_zz = NULL ; }
#else
#  define CLONE_VEC(zz_type_zz,zz_src_zz,zz_dst_zz,zz_n_zz) \
	    zz_dst_zz = copy<zz_type_zz>(zz_src_zz,zz_n_zz)
#endif


/* Specialised copy functions for complex dylp structures. */


/*! \brief Specialised copy function for a dylp basis_struct

  This routine expects that the basis will come already provisioned with a
  basis element vector of the proper capacity.
*/
  
inline void ODSI::copy_basis (const basis_struct* src, basis_struct* dst)

{ if (!src) return ;

  dst->len = src->len ;
  if (dst->el == 0)
  { throw CoinError("No basis element vector",
		    "copy_basis","OsiDylpSolverInterface") ; }
  memcpy(dst->el,src->el,idx(src->len)*sizeof(basisel_struct)) ;

# if ODSI_PARANOIA >= 2
  assert_same(*dst, *src, true) ;
# endif

  return ; }


/*! \brief Specialised copy function for a dylp basis_struct

  This routine will allocate a basis element vector of the proper capacity.
  NOTE that the capacity must equal or exceed lpprob_struct->rowsze for the
  lpprob_struct that this basis is attached to.  We use CALLOC here to
  minimize mixing of memory allocation with new & alloc.  This structure
  stands a good chance of being reallocated within dylp. Recall that this is
  1-based addressing, hence the need to allocate (and copy) one additional
  element.
*/

inline basis_struct* ODSI::copy_basis (const basis_struct* src, int dstsze)

{ if (!src) return 0 ;
  basis_struct* dst = new basis_struct ;
  dst->el = (basisel_struct *) CALLOC(idx(dstsze),sizeof(basisel_struct)) ;
  copy_basis(src,dst) ;
  return dst ; }


/*! \brief Specialised copy function for a dylp lpprob_struct */

lpprob_struct* ODSI::copy_lpprob (const lpprob_struct* src)

{ if (!src) return 0 ;

  int col_count = idx(src->colsze) ;
  int row_count = idx(src->rowsze) ;

  lpprob_struct* dst = NULL;
  CLONE(lpprob_struct,src,dst);

  dst->basis = copy_basis(src->basis,row_count) ;
  CLONE_VEC(flags,src->status,dst->status,col_count);
  CLONE_VEC(double,src->x,dst->x,row_count);
  CLONE_VEC(double,src->y,dst->y,row_count);
  CLONE_VEC(bool,src->actvars,dst->actvars,col_count);

# if ODSI_PARANOIA >= 2
  assert_same(*dst, *src, true) ;
# endif

  return dst ; }

//@} // CopyHelpers



/* Cache functions

   These properly belong in the DestructorHelpers group, but they need to be
   up here so that they are declared before the first use.
*/

/*! \brief Destroy cached values
    \ingroup DestructorHelpers

  Destroy cached copies of computed values for columns. If structure is true,
  structural values are destroyed along with solution values.
  See the \link OsiDylpSolverInterface.cpp comments \endlink at the head of
  the file for a brief description of the ODSI cache mechanism.
*/

inline void ODSI::destruct_col_cache (bool structure)

{ delete [] _col_x ; _col_x = 0 ;
  delete [] _col_cbar ; _col_cbar = 0 ;
  if (structure == true)
  { delete [] _col_obj ; _col_obj = 0 ; }

  return ; }


/*! \brief Destroy cached values
    \ingroup DestructorHelpers

  Destroy cached copies of computed values for rows. If structure is true,
  structural values are destroyed along with solution values.
  See the \link OsiDylpSolverInterface.cpp comments \endlink at the head of
  the file for a brief description of the ODSI cache mechanism.
*/

void ODSI::destruct_row_cache (bool structure)

{ delete [] _row_price ; _row_price = 0 ;
  delete [] _row_lhs ; _row_lhs = 0 ;

  if (structure == true)
  { delete [] _row_lower ; _row_lower = 0 ;
    delete [] _row_range ; _row_range = 0 ;
    delete [] _row_rhs ; _row_rhs = 0 ;
    delete [] _row_sense ; _row_sense = 0 ;
    delete [] _row_upper ; _row_upper = 0 ; }
  
  return ; }

/*! \brief Destroy cached values
    \ingroup DestructorHelpers

  Destroy cached copies of computed values. If either rowStructure or
  colStructure is true, take it as an indication we need to delete the cached
  copies of the constraint system. If you want to remove structural aspects
  of the row or column cache, without affecting the matrices, call them
  individually.
  See the \link OsiDylpSolverInterface.cpp file comments \endlink
  for a brief description of the ODSI cache mechanism.
*/
  
inline void ODSI::destruct_cache (bool rowStructure, bool colStructure)

{ destruct_row_cache(rowStructure) ;
  destruct_col_cache(colStructure) ;

  if (rowStructure == true || colStructure == true)
  { delete _matrix_by_row ; _matrix_by_row = 0 ;
    delete _matrix_by_col ; _matrix_by_col = 0 ; }

  return ; }




/*! \defgroup ConsysHelpers Helper Functions for Problem Modification
    \brief Private helper functions for problem modification

    This group of functions implement conversions between the OSI notion of
    a constraint system and the dylp notion of a constraint system, and
    calculations related to modification of the primal solution or objective.
*/
//@{

/*! \brief Determine dylp constraint type from constraint bounds. */

inline contyp_enum ODSI::bound_to_type (double lower, double upper)

{ 

  if (upper == lower) return contypEQ ;
  bool finite_low = (lower > -odsiInfinity) ;
  bool finite_up = (upper < odsiInfinity) ;
  if (finite_low && finite_up) return contypRNG ;
  if (finite_low) return contypGE ;
  if (finite_up) return contypLE ;
  return contypNB ; }


/*! \brief Convert OSI row sense to dylp constraint type */

inline contyp_enum ODSI::sense_to_type (char sense)

{ switch (sense)
  { case 'E': return contypEQ ;
    case 'G': return contypGE ;
    case 'L': return contypLE ;
    case 'N': return contypNB ;
    case 'R': return contypRNG ;
    default: { assert(0) ; return contypINV ; } } }


/*! \brief Convert dylp constraint type to OSI row sense */

inline char ODSI::type_to_sense (contyp_enum ctypi)

{ switch (ctypi)
  { case contypNB:
    { return 'N' ; }
    case contypGE:
    { return 'G' ; }
    case contypEQ:
    { return 'E' ; }
    case contypLE:
    { return 'L' ; }
    case contypRNG:
    { return 'R' ; }
    case contypINV:
    default: 
    { assert(0) ;
      return '?' ; } } }


/*! \brief Convert OSI sense/rhs/range description for a single constraint to
	   dylp type/rhs/rhslow form
  
  This routine generates the dylp rhs, rhslow, and ctyp values corresponding
  to the supplied OSI sense, rhs, and range values.
*/

inline void
ODSI::gen_rowiparms (contyp_enum* ctypi, double* rhsi, double* rhslowi,
		     char sensei, double rhsini, double rangei)

{ *ctypi = sense_to_type(sensei) ;
  switch (*ctypi)
  { case contypEQ:
    { *rhslowi = 0 ;
      *rhsi = rhsini ;
      break ; }
    case contypLE:
    { *rhslowi = 0 ;
      *rhsi = rhsini ;
      break ; }
    case contypGE:
    { *rhslowi = 0 ;
      *rhsi = rhsini ;
      break ; }
    case contypRNG:
    { *rhslowi = rhsini-rangei ;
      *rhsi = rhsini ;
      break ; }
    case contypNB:
    { *rhslowi = 0 ;
      *rhsi = 0 ;
      break ; }
    default:
    { assert(false) ; } }

  return ; }


/*! \brief Convert OSI upper/lower bound description for a single constraint to
	   dylp type/rhs/rhslow form
  
  This routine generates the dylp rhs, rhslow, and ctyp values corresponding
  to the supplied constraint upper and lower bound values.
*/

inline void
ODSI::gen_rowiparms (contyp_enum* ctypi, double* rhsi, double* rhslowi,
		     double rowlbi, double rowubi)

{ *ctypi = bound_to_type(rowlbi,rowubi) ;
  switch (*ctypi)
  { case contypEQ:
    case contypLE:
    { *rhslowi = 0 ;
      *rhsi = rowubi ;
      break ; }
    case contypGE:
    { *rhslowi = 0 ;
      *rhsi = rowlbi ;
      break ; }
    case contypRNG:
    { *rhslowi = rowlbi ;
      *rhsi = rowubi ;
      break ; }
    case contypNB:
    { *rhslowi = 0 ;
      *rhsi = 0 ;
      break ; }
    default:
    { assert(false) ; } }

  return ; }


/*! \brief Add a column to the constraint matrix.

  Convert the OSI packed vector to a dylp packed vector and install the
  column. The remaining parameters are expected to be in the correct dylp
  form.

  Note that the ODSI client may be building a constraint system from scratch
  using calls to addCol/addRow. There's no guarantee that consys exists yet.
*/

void ODSI::add_col (const CoinPackedVectorBase& coin_colj,
		    vartyp_enum vtypj, double vlbj, double vubj, double objj)

{ pkvec_struct *pk_colj = packed_vector(coin_colj,getNumRows()) ;

/*
  Do we need a consys?
*/
  if (!consys) construct_consys(0,0) ;
/*
  Generate a name. Not a trivial issue; unique names are critical to the
  ability to write out a system in mps format.
*/
  std::ostringstream nme ;
  nme << "var" ;
  nme << addedColCnt++ ;
  std::string colname = nme.str() ;
  pk_colj->nme = colname.c_str() ;
/*
  Add the column.
*/
  bool r = consys_addcol_pk(consys,vtypj,pk_colj,objj,vlbj,vubj) ;
  pkvec_free(pk_colj) ;
  if (!r)
  { lp_retval = lpFATAL ; }
/*
  After adding a column, the best we can do is a warm start. Any current
  solution is no longer valid, and we need to delete the structural side
  of the column cache.
*/
  resolveOptions->forcewarm = true ;
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::add_col: new column." << std::endl ;
# endif

  destruct_cache(false,true) ; }


/*! \brief Add a row to the constraint matrix.

  Convert the OSI packed vector to a dylp packed vector and install the row.
  The remaining parameters are expected to be in the correct dylp form.

  Note that the ODSI client may be building a constraint system from scratch
  using calls to addCol/addRow. There's no guarantee that consys exists yet.
*/

void ODSI::add_row (const CoinPackedVectorBase &coin_rowi, char clazzi,
		    contyp_enum ctypi, double rhsi, double rhslowi)

{ pkvec_struct *pk_rowi = packed_vector(coin_rowi,getNumCols()) ;

/*
  Do we need a consys?
*/
  if (!consys) construct_consys(0,0) ;
/*
  Generate a name. Not a trivial issue; unique names are critical to the
  ability to write out a system in mps format.
*/
  std::ostringstream nme ;
  if (clazzi == 'a')
  { nme << "con" ; }
  else
  { nme << "cut" ; }
  nme << addedRowCnt++ ;
  std::string rowname = nme.str() ;
  pk_rowi->nme = rowname.c_str() ;
/*
  Add the row.
*/
  bool r = consys_addrow_pk(consys,clazzi,ctypi,
				   pk_rowi,rhsi,rhslowi,0,0) ;
  pkvec_free(pk_rowi) ;
  if (!r)
  { lp_retval = lpFATAL ; }
/*
  After adding a constraint, the best we can do is a warm start. Any current
  solution is no longer valid, and we need to delete the structural side of the
  row cache.
*/
  resolveOptions->forcewarm = true ;
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::add_row: new row." << std::endl ;
# endif

  destruct_cache(true,false) ; }


/*! \brief Establish a pessimistic primal solution.

  This routine sets up a pessimistic basic primal solution. The logicals are
  assigned to the basis, and the architecturals are set to the (finite) bound
  that produces the worst objective function value. We can't claim this is
  a worst-case bound on the objective because we're restricted to choosing a
  finite bound for the nonbasic variables, if one exists.

  Feasibility (with respect to the architectural constraints) is not even on
  the radar.

  The interface design guarantees that if a consys_struct exists, then
  the objective and bounds vectors are already attached to it.
*/

void ODSI::pessimal_primal ()

{ 

/*
  No columns, no solution. A non-zero value guarantees a consys structure.
*/
  int n = getNumCols() ;
  if (n == 0) return ;
  assert(consys && consys->obj && consys->vub && consys->vlb) ;
  int m = getNumRows() ;

  double *obj = consys->obj ;
  double *vlb = consys->vlb ;
  double *vub = consys->vub ;

  double lbj,ubj,xj ;
  CoinWarmStartBasis::Status statj ;
/*
  We have columns. Replace or allocate the cached primal solution vector.
*/
  if (_col_x) delete[] _col_x ;
  _col_x = new double[n] ;
/*
  Create a basis of the appropriate size.
*/
  OsiDylpWarmStartBasis *pessbasis =
    dynamic_cast<OsiDylpWarmStartBasis *>(getEmptyWarmStart()) ;
  pessbasis->setSize(n,m) ;
/*
  Walk the objective vector, taking values for variables from the appropriate
  bounds vector. The objective attached to consys is always a minimisation
  objective.
*/
  for (int j = 1 ; j <= n ; j++)
  { lbj = vlb[j] ;
    ubj = vub[j] ;
    if (lbj > -odsiInfinity && ubj < odsiInfinity)
    { if (obj[j] > 0)
      { xj = ubj ;
	statj = CoinWarmStartBasis::atUpperBound ; }
      else
      { xj = lbj ;
	statj = CoinWarmStartBasis::atLowerBound ; } }
    else
    if (lbj > -odsiInfinity)
    { xj = lbj ;
      statj = CoinWarmStartBasis::atLowerBound ; }
    else
    if (ubj < odsiInfinity)
    { xj = ubj ;
      statj = CoinWarmStartBasis::atLowerBound ; }
    else
    { xj = 0 ;
      statj = CoinWarmStartBasis::isFree ; }
    _col_x[inv(j)] = xj ;
    pessbasis->setStructStatus(inv(j),statj) ; }
/*
  Set the logicals (aka artificials) to be basic and mark all constraints as
  active.
*/
  for (int i = 1 ; i <= m ; i++)
  { pessbasis->setArtifStatus(inv(i),CoinWarmStartBasis::basic) ;
    pessbasis->setConStatus(inv(i),CoinWarmStartBasis::atLowerBound) ; }
/*
  Set pessbasis to be the active basis and return.
*/
# if ODSI_TRACK_ACTIVE > 0
  std::cout
    << "ODSI(" << std::hex << this
    << ")::pessimal_primal: replacing active basis "
    << activeBasis << " with " << pessbasis << std::dec
    << "." << std::endl ;
# endif
  delete activeBasis ;
  activeIsModified = true ;
  activeBasis = pessbasis ;
  
  return ; }

/*! \brief Calculate objective value using cached objective function and
	   primal solution.

  This routine refreshes the cached objective value.  It's intended for use
  after client changes to the the objective coefficients or primal solution
  values.  Note that the cached objective value is updated directly from the
  solver return value whenever the solver is called.
*/

void ODSI::calc_objval ()

{ int n = getNumCols() ;

/*
  The easy case --- if we have no variables, the objective must be 0.
*/
  if (n == 0)
  { _objval = 0 ;
    return ; }
/*
  We have variables. Grab the primal solution and objective function and
  calculate their dot product.
*/
  const double *obj = getObjCoefficients() ;
  const double *x = getColSolution() ;
  _objval = 0.0 ;
  for (int j = 0 ; j < n ; j++) _objval += x[j]*obj[j] ;
  setcleanzero(_objval,tolerances->cost) ;

  return ; }



//@} // ConsysHelpers



/* OsiDylpSolverInterface constructors and related subroutines. */


/*! \defgroup ConstructorHelpers Helper functions for problem construction */
//@{

/*! \brief Construct a dylp lpprob_struct (LP problem).  */

void ODSI::construct_lpprob ()

{ lpprob = new lpprob_struct ;
  memset(lpprob, 0, sizeof(lpprob_struct)) ;
  setflg(lpprob->ctlopts,lpctlNOFREE) ;
  lpprob->phase = dyINV ;
  lpprob->consys = consys ;
  lpprob->rowsze = consys->rowsze ;
  lpprob->colsze = consys->colsze ;
  
  return ; }


/*! \brief Construct and initialise default options and tolerances.
*/

inline void ODSI::construct_options ()

{ 
/*
  Acquire the default options and tolerances from dylp. Set the OSI options
  and tolerances to match. Turn dylp's printing options down to `catatonic'.
*/
  delete initialSolveOptions ;
  initialSolveOptions = new lpopts_struct ;
  delete tolerances ;
  tolerances = new lptols_struct ;
  dy_defaults(&initialSolveOptions,&tolerances) ;
  tolerances->inf = odsiInfinity ;
  delete resolveOptions ;
  CLONE(lpopts_struct,initialSolveOptions,resolveOptions);
  dy_setprintopts(0,initialSolveOptions) ;
  dy_setprintopts(0,resolveOptions) ;

  setIntParam(OsiMaxNumIteration,3*initialSolveOptions->iterlim) ;
  setIntParam(OsiMaxNumIterationHotStart,3*initialSolveOptions->iterlim) ;
  setDblParam(OsiDualTolerance,tolerances->dfeas_scale*tolerances->cost) ;
  setDblParam(OsiPrimalTolerance,tolerances->pfeas_scale*tolerances->zero) ;
/*
  Differentiate options for initial solution vs. reoptimising.
*/
  initialSolveOptions->forcecold = true ;
  initialSolveOptions->fullsys = true ;
  resolveOptions->forcecold = false ;
  resolveOptions->fullsys = false ;

  return ; }

/*! \brief Construct an empty constraint system of the specified size.

  Construct an empty constraint system of the specified capacity, with all
  the appropriate options and attached vectors. If the specified capacity
  is 0, you'll get dylp's (moderately large) default values (capacity for a
  thousand or so constraints and variables).
*/

void ODSI::construct_consys (int cols, int rows)

{
/*
  Set parameters to specify the appropriate options and attached vectors.
  
  parts specifies (in order) an objective, variable upper bounds, variable
  lower bounds, constraint rhs, range constraint rhs, variable type, and
  constraint type.

  options specifies that consys should issue a warning if requested to attach
  a vector that's already attached.
*/
  flags parts = CONSYS_OBJ | 
		CONSYS_VUB | CONSYS_VLB | CONSYS_VTYP |
		CONSYS_RHS | CONSYS_RHSLOW | CONSYS_CTYP ;
  flags opts = CONSYS_WRNATT ;
/*
  The first parameter to consys is the constraint system name, which is
  never conveniently available in the OSI frame.
*/
  consys = consys_create(0,parts,opts,rows,cols,odsiInfinity) ;
  if (consys == 0)
  { lp_retval = lpFATAL ; }

  return ; }


/*! \brief Initialise dylp i/o subsystems

  Initialise dylp i/o subsystems for normal and error output. This should be
  done only once.
*/

void ODSI::dylp_ioinit ()

{ if (reference_count > 1) return ;

  string errfile = string(DYLP_ERRMSGDIR)+string("/dy_errmsgs.txt") ;
# ifdef ODSI_INFOMSGS
  errinit(const_cast<char *>(errfile.c_str()),0,true) ;
# else
  errinit(const_cast<char *>(errfile.c_str()),0,false) ;
# endif
  bool r1 UNUSED = dyio_ioinit() ;
  assert(r1) ;
}

/*! \brief Load a problem description

  This routine expects a CoinMpsIO object, and proceeds to extract the
  information from the object. This permits problem, constraint, and
  variable names to be passed in to dylp, which makes for much friendlier
  output.

  The routine is built on load_problem(matrix,sense/rhs/range), so the
  first thing it does is unpack the expected operands from the CoinMpsIO
  object.

  The next action is to destroy the existing problem. Existing options and
  tolerances settings are not affected.

  The routine then builds an empty consys_struct of the proper size with a
  call to construct_consys. The main body of the routine fills the constraint
  matrix in two loops. In the first, calls to consys_addrow_pk insert empty
  constraints and establish the constraint type, rhs, and (optional) range.
  Then calls to consys_addcol_pk insert the variables, along with their
  bounds and objective coefficient.

  Finally, a primal solution is invented, and the objective is set accordingly.
*/

void ODSI::load_problem (const CoinMpsIO &mps)

{ int n = mps.getNumCols() ;
  int m = mps.getNumRows() ;
/*
  Unpack the operands from the CoinMpsIO object.
*/
  const CoinPackedMatrix matrix = *mps.getMatrixByCol() ;
  const double* col_lower = mps.getColLower() ;
  const double* col_upper = mps.getColUpper() ;
  const double* obj = mps.getObjCoefficients() ;
  const char *sense = mps.getRowSense() ;
  const double* rhsin = mps.getRightHandSide() ;
  const double* range = mps.getRowRange() ;

  double *rhs = new double[m] ;
  double *rhslow = new double[m] ; 
  contyp_enum *ctyp = new contyp_enum[m] ;

  gen_rowparms(m,rhs,rhslow,ctyp,sense,rhsin,range) ;
/*
  Free the existing problem structures, preserving options and tolerances.
*/
  destruct_problem(true) ;
/*
  Create an empty consys_struct. Load the constraint system name, objective
  name, and objective offset.
*/
  construct_consys(n,m) ;
  setStrParam(OsiProbName,mps.getProblemName()) ;
  if (consys->objnme) strfree(consys->objnme) ;
  consys->objnme = stralloc(mps.getObjectiveName()) ;
  setDblParam(OsiObjOffset,mps.objectiveOffset()) ;
/*
  First loop: Insert empty constraints into the new constraint system.
*/
  pkvec_struct* rowi = pkvec_new(0) ;
  assert(rowi) ;

  bool r = true ;
  for (int i = 0 ; i < m ; i++)
  { rowi->nme = const_cast<char *>(mps.rowName(i)) ;
    r = consys_addrow_pk(consys,'a',ctyp[i],rowi,rhs[i],rhslow[i],0,0) ;
    if (!r)
    { lp_retval = lpFATAL ;
      break ; } }
  
  if (rowi) pkvec_free(rowi) ;
  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
/*
  Take a moment and set up a vector of variable types. It's a bit of a pain
  that CoinMpsIO doesn't distinguish binary from general integer, but we can
  do it here by checking bounds. The test is drawn from OSI::isBinary and
  will include variables fixed at 0 or 1.
*/
  const char *const intvars = mps.integerColumns() ;
  vartyp_enum *const vtyp = new vartyp_enum[n] ;
  if (intvars)
  { for (int j = 0 ; j < n ; j++)
    { if (intvars[j])
      { if ((col_lower[j] == 0.0 || col_lower[j] == 1.0) &&
	    (col_upper[j] == 0.0 || col_upper[j] == 1.0))
	{ vtyp[j] = vartypBIN ; }
	else
	{ vtyp[j] = vartypINT ; } }
      else
      { vtyp[j] = vartypCON ; } } }
  else
  { for (int j = 0 ; j < n ; j++) vtyp[j] = vartypCON ; }
/*
  Second loop. Insert the coefficients by column. If we need to create a
  column-ordered copy, take advantage and cache it.
*/
  const CoinPackedMatrix *matrix2 ;
  if (matrix.isColOrdered())
  { matrix2 = &matrix ; }
  else
  { _matrix_by_col = new CoinPackedMatrix ;
    _matrix_by_col->reverseOrderedCopyOf(matrix) ;
    matrix2 = _matrix_by_col ; }
  
  pkvec_struct* colj = pkvec_new(m) ;

  for (int j = 0 ; j < n ; j++)
  { const CoinShallowPackedVector coin_col = matrix2->getVector(j) ;
    packed_vector(coin_col,n,colj) ;
    colj->nme = const_cast<char *>(mps.columnName(j)) ;
    r = consys_addcol_pk(consys,vtyp[j],
			 colj,obj[j],col_lower[j],col_upper[j]) ;
    if (!r)
    { break ; } }

  pkvec_free(colj) ;
  delete[] vtyp ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
  assert(matrix2->isEquivalent(*getMatrixByCol())) ;

/*
  Construct a pessimal solution and we're done.
*/
  pessimal_primal() ;
  calc_objval() ;

  return ; }

/*! \brief Load a problem description

  This routine expects a constraint system described in terms of an OSI
  packed matrix and dylp constraint type (ctyp) and right-hand-side (rhs,
  rhslow) vectors.

  The first action is to destroy the existing problem. Existing options and
  tolerances settings are not affected.

  Next the routine builds an empty consys_struct of the proper size with a
  call to construct_consys. The main body of the routine fills the constraint
  matrix in two loops. In the first, calls to consys_addrow_pk insert empty
  constraints and establish the constraint type, rhs, and (optional) range.
  Then calls to consys_addcol_pk insert the variables, along with their
  bounds and objective coefficient.

  Finally, a primal solution is invented, and the objective is set accordingly.
*/

void ODSI::load_problem (const CoinPackedMatrix& matrix,
    const double* col_lower, const double* col_upper, const double* obj,
    const contyp_enum *ctyp, const double* rhs, const double* rhslow)

{ 
/*
  Free the existing problem structures, preserving options and tolerances.
*/
  destruct_problem(true) ;
/*
  We're expecting to work with a column-ordered matrix, so make sure that's
  true. If we need to create a column-ordered copy, take advantage and cache
  it.
*/
  const CoinPackedMatrix *matrix2 ;
  if (matrix.isColOrdered())
  { matrix2 = &matrix ; }
  else
  { _matrix_by_col = new CoinPackedMatrix ;
    _matrix_by_col->reverseOrderedCopyOf(matrix) ;
    matrix2 = _matrix_by_col ; }
/*
  Create an empty consys_struct.
*/
  int m = matrix2->getNumRows() ;
  int n = matrix2->getNumCols() ;

  construct_consys(n,m) ;
/*
  First loop: Insert empty constraints into the new constraint system.
*/
  pkvec_struct* rowi = pkvec_new(0) ;
  assert(rowi) ;

  bool r = true ;
  for (int i = 0 ; i < m ; i++)
  { rowi->nme = 0 ;
    r = consys_addrow_pk(consys,'a',ctyp[i],rowi,rhs[i],rhslow[i],0,0) ;
    if (!r)
    { break ; } }

  if (rowi) pkvec_free(rowi) ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
/*
  Second loop. Insert the coefficients by column.
*/
  pkvec_struct* colj = pkvec_new(m) ;

  for (int j = 0 ; j < n ; j++)
  { const CoinShallowPackedVector coin_col = matrix2->getVector(j) ;
    packed_vector(coin_col,n,colj) ;
    double objj = obj?obj[j]:0 ;
    double vlbj = col_lower?col_lower[j]:0 ;
    double vubj = col_upper?col_upper[j]:odsiInfinity ;
    colj->nme = 0 ;
    r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    if (!r)
    { break ; } }

  pkvec_free(colj) ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
  assert(matrix2->isEquivalent(*getMatrixByCol())) ;
/*
  Construct a pessimal solution and we're done.
*/
  pessimal_primal() ;
  calc_objval() ;

  return ; }


/*! \brief Load a problem description

  This routine expects a constraint system described using a standard column-
  major packed matrix structure and dylp constraint sense and right-hand-side
  (rhs, rhslow) vectors.

  The matrix description consists of row and column size and four arrays.
  Coefficients and corresponding row indices are given in value and index,
  respectively (the bulk storage). The vectors start and len contain the
  starting position and length, respectively, of each column in the bulk
  storage arrays. By passing in an explicit len array, we can handle both
  fully packed matrices (the usual case) and loosely packed matrices
  (as produced by presolve).

  The first action is to free the existing problem.  Existing options and
  tolerances settings are saved, if they exist, and will be restored after
  the new problem is loaded.

  Next the routine builds an empty consys_struct of the proper size with a
  call to construct_consys. The main body of the routine fills the constraint
  matrix in two loops. In the first, calls to consys_addrow_pk insert empty
  constraints and establish the constraint type and rhs. Then calls to
  consys_addcol_pk insert the variables, along with their bounds and
  objective coefficient.

  Finally, an lpprob_struct is created to hold the new problem, and options
  and tolerances are reattached or initialised from defaults.
*/

void ODSI::load_problem (const int colcnt, const int rowcnt,
	   const int *start, const int *len,
	   const int *index, const double *value,
	   const double* col_lower, const double* col_upper, const double* obj,
	   const contyp_enum *ctyp, const double* rhs, const double* rhslow)

{ 
/*
  Free the existing problem structures, preserving options and tolerances.
*/
  destruct_problem(true) ;
/*
  Create an empty consys_struct. Request the standard attached vectors:
  objective, variable upper & lower bounds, variable and constraint types,
  and constraint right-hand-side and range (rhslow).
*/
  construct_consys(rowcnt,colcnt) ;
/*
  First loop: Insert empty constraints into the new constraint system.
*/
  pkvec_struct* rowi = pkvec_new(0) ;
  assert(rowi) ;

  bool r = true ;
  for (int i = 0 ; i < rowcnt ; i++)
  { rowi->nme = 0 ;
    r = consys_addrow_pk(consys,'a',ctyp[i],rowi,rhs[i],rhslow[i],0,0) ;
    if (!r)
    { break ; } }

  if (rowi) pkvec_free(rowi) ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
/*
  Second loop. Insert the coefficients by column. The size of colj is a gross
  overestimate, but it saves the trouble of constantly reallocating it.
*/
  pkvec_struct *colj = pkvec_new(rowcnt) ;
  assert(colj) ;
  colj->dim = rowcnt ;
  pkcoeff_struct *coeffs = colj->coeffs ;

  for (int j = 0 ; j < colcnt ; j++)
  { int startj = start[j] ;
    int lenj = (len)?(len[j]):(start[j+1]-startj) ;

    for (int ndx = 0 ; ndx < lenj ; ndx++)
    { coeffs[ndx].ndx = idx(index[startj+ndx]) ;
      coeffs[ndx].val = value[startj+ndx] ; }
    colj->cnt = lenj ;
  
    double objj = obj?obj[j]:0 ;
    double vlbj = col_lower?col_lower[j]:0 ;
    double vubj = col_upper?col_upper[j]:odsiInfinity ;
    colj->nme = 0 ;
    r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    if (!r)
    { break ; } }

  if (colj) pkvec_free(colj) ;
  if (!r)
  { lp_retval = lpFATAL ;
    return ; }
/*
  Construct a pessimal solution and we're done.
*/
  pessimal_primal() ;
  calc_objval() ;

  return ; }



/*! \brief Generate dylp constraint parameters from upper and lower bounds.

  This routine generates the dylp rhs, rhslow, and ctyp vectors given
  constraint upper and lower bound vectors. rowlb[i] and rowub[i] default to
  -inf and inf, respectively. 
*/

void ODSI::gen_rowparms (int rowcnt,
			 double *rhs, double *rhslow, contyp_enum *ctyp,
			 const double *rowlb, const double *rowub)

{ double rowlbi,rowubi ;

  for (int i = 0 ; i < rowcnt ; i++)
  { rowlbi = rowlb?rowlb[i]:-odsiInfinity ;
    rowubi = rowub?rowub[i]:odsiInfinity ;
    ctyp[i] = bound_to_type(rowlbi,rowubi) ;

    switch (ctyp[i])
    { case contypEQ:
      case contypLE:
      { rhs[i] = rowubi ; 
	rhslow[i] = 0 ;
	break ; }
      case contypGE:
      { rhs[i] = rowlbi ;
	rhslow[i] = 0 ;
	break ; }
      case contypRNG:
      { rhs[i] = rowubi ; 
	rhslow[i] = rowlbi ;
	break ; }
      case contypNB:
      { rhs[i] = odsiInfinity ;
	rhslow[i] = -odsiInfinity ;
	break ; }
      default:
      { assert(0) ; } } }

  return ; }


/*! \brief Generate dylp constraint parameters from rhs, range, and row sense

  This routine generates the dylp rhs, rhslow, and ctyp vectors given
  constraint right-hand-side (rhs), range, and sense vectors. Rhs and range
  default to 0, and sense defaults to 'G' (>=).
*/

void ODSI::gen_rowparms (int rowcnt,
			 double *rhs, double *rhslow, contyp_enum *ctyp,
		 const char *sense, const double *rhsin, const double *range)

{ double rhsi,rangei ;
  char sensei ;

  for (int i = 0 ; i < rowcnt ; i++)
  { rhsi = rhsin?rhsin[i]:0 ;
    rangei = range?range[i]:0 ;
    sensei = sense?sense[i]:'G' ;
    
    switch (sensei)
    { case 'E':
      { rhs[i] = rhsi ;
	rhslow[i] = 0 ;
	ctyp[i] = contypEQ ;
	break ; }
      case 'L':
      { rhs[i] = rhsi ;
	rhslow[i] = 0 ;
	ctyp[i] = contypLE ;
	break ; }
      case 'G':
      { rhs[i] = rhsi ;
	rhslow[i] = 0 ;
	ctyp[i] = contypGE ;
	break ; }
      case 'R':
      { rhs[i] = rhsi ;
	rhslow[i] = rhsi-rangei ;
	ctyp[i] = contypRNG ;
	break ; }
      case 'N':
      { rhs[i] = 0 ;
	rhslow[i] = 0 ;
	ctyp[i] = contypNB ;
	break ; }
      default:
      { assert(0) ; } } }

  return ; }

//@} // ConstructorHelpers



/*! \defgroup ODSIConstructorsDestructors ODSI Constructors and Destructors
    \brief ODSI Constructors and destructors

    ODSI provides a default constructor and a copy constructor, and a default
    destructor.
*/
//@{

/*!
  Creates the shell of an ODSI object. The options and tolerances structures
  will be instantiated and initialised, but nothing else specific to dylp.
  Message and hint structures inherited from OsiSolverInterface are
  initialised.

  Creation of the first ODSI object triggers initialisation of the i/o
  subsystem.
*/

ODSI::OsiDylpSolverInterface ()

  : initialSolveOptions(0),
    resolveOptions(0),
    tolerances(0),
    consys(0),
    lpprob(0),
    statistics(0),

    local_outchn(IOID_NOSTRM),
    local_logchn(IOID_NOSTRM),
    initial_gtxecho(false),
    resolve_gtxecho(false),
    lp_retval(lpINV),
    obj_sense(1.0),
    odsiInfinity(CoinInfinity),

    solvername("dylp"),
    mps_debug(0),
    hotstart_fallback(0),
    activeBasis(0),
    activeIsModified(false),
    solnIsFresh(false),
    addedColCnt(0),
    addedRowCnt(0),

    _objval(0),
    _col_obj(0),
    _col_x(0),
    _col_cbar(0),
    _row_rhs(0),
    _row_lower(0),
    _row_upper(0),
    _row_sense(0),
    _row_range(0),
    _row_lhs(0),
    _row_price(0),
    _matrix_by_col(0),
    _matrix_by_row(0),

    preObj_(0),
    postActions_(0),
    postObj_(0),
    passLimit_(5),
    keepIntegers_(false),
    savedConsys_(0),
    saved_col_obj(0),
    saved_row_rhs(0),
    saved_row_lower(0),
    saved_row_upper(0),
    saved_row_sense(0),
    saved_row_range(0),
    saved_matrix_by_col(0),
    saved_matrix_by_row(0)

{
/*
  Replace the OSI default messages with ODSI messages.
*/
  setOsiDylpMessages(CoinMessages::us_en) ;
/*
  Clear the hint info_ array, then turn on the hints that should be on by
  default.
*/
  for (int i = 0 ; i < OsiLastHintParam ; i++) info_[i] = 0 ;
  setHintParam(OsiDoPresolveInInitial,true,OsiForceDo,0) ;
/*
  Acquire the default options and tolerances.
*/
  construct_options() ;
/*
  Increment the ODSI existence count. If this is the first interface in
  existence, initialise the i/o package.
*/
  reference_count++ ;

  if (reference_count == 1)
  { dylp_ioinit() ;
    CoinRelFltEq eq ;
    assert(eq(odsiInfinity, odsiInfinity)) ; }

# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << "): default constructor." << std::endl ;
# endif

  return ; }


/*!
  A true copy --- no data structures are shared with the original. The
  statistics structure, hot start structure, and presolve structures are not
  copied (but presolve control information is copied). The copy does not
  inherit open files from the original. Cached information is not replicated.
*/

ODSI::OsiDylpSolverInterface (const OsiDylpSolverInterface& src)

  : OsiSolverInterface(src),
    statistics(0),		// do not copy statistics
    local_outchn(IOID_NOSTRM),  // do not copy open file descriptors
    local_logchn(IOID_NOSTRM),
    initial_gtxecho(src.initial_gtxecho),
    resolve_gtxecho(src.resolve_gtxecho),
    lp_retval(src.lp_retval),
    obj_sense(src.obj_sense),
    odsiInfinity(src.odsiInfinity),

    solvername(src.solvername),
    mps_debug(src.mps_debug),
    hotstart_fallback(0),	// do not copy hot start information
    activeBasis(0),
    activeIsModified(false),
    addedColCnt(src.addedColCnt),
    addedRowCnt(src.addedRowCnt),
  
    _objval(src._objval),
    _col_obj(0),
    _col_x(0),
    _col_cbar(0),
    _row_rhs(0),
    _row_lower(0),
    _row_upper(0),
    _row_sense(0),
    _row_range(0),
    _row_lhs(0),
    _row_price(0),
    _matrix_by_col(0),
    _matrix_by_row(0),

    preObj_(0),			// do not copy presolve structures
    postActions_(0),
    postObj_(0),
    passLimit_(src.passLimit_),
    keepIntegers_(src.keepIntegers_),
    savedConsys_(0),
    saved_col_obj(0),
    saved_row_rhs(0),
    saved_row_lower(0),
    saved_row_upper(0),
    saved_row_sense(0),
    saved_row_range(0),
    saved_matrix_by_col(0),
    saved_matrix_by_row(0)

{ if (src.consys)
  { bool r UNUSED = consys_dupsys(src.consys,&consys,src.consys->parts) ;
    assert(r) ; }
  else
  { consys = 0 ; }
  if (src.lpprob)
  { lpprob = copy_lpprob(src.lpprob) ;
    lpprob->consys = consys ; }
  else
  { lpprob = 0 ; }
  solnIsFresh = src.solnIsFresh ;

  CLONE(lpopts_struct,src.initialSolveOptions,initialSolveOptions) ;
  CLONE(lpopts_struct,src.resolveOptions,resolveOptions) ;
  CLONE(lptols_struct,src.tolerances,tolerances) ;

  if (src.activeBasis)
  { activeBasis = src.activeBasis->clone() ;
    activeIsModified = src.activeIsModified ; }

  int n = getNumCols() ;
  int m = getNumRows() ;
  CLONE_VEC(double,src._col_x,_col_x,n) ;
  CLONE_VEC(double,src._row_price,_row_price,m) ;

  COPY_VEC(void *,&src.info_[0],&info_[0],OsiLastHintParam) ;

  reference_count++ ;

# if ODSI_PARANOIA >= 2
  assert_same(*this, src, true) ;
# endif
# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << "): copy from "
    << &src << std::dec << "." << std::endl ;
# endif

}


/*!
  An alternate copy constructor; the parameter is ignored (and assumed true).
*/

inline OsiSolverInterface* ODSI::clone (bool copyData) const

{ 
# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec << "): cloning ("
    << ((this->solnIsFresh == true)?"fresh":"stale") << ")." << std::endl ;
# endif

  if (copyData)
  { return new OsiDylpSolverInterface(*this) ; }
  else
  { return new OsiDylpSolverInterface() ; }
}


/*! Assignment operator.

  Much as for the copy constructor, except that we trash the contents of the
  lhs object before we start. As with copy, statistics and hot start
  information are not copied over to the assignment target, nor does it inherit
  open files or terminal echo settings. Cached information is not replicated.
  Presolve information simply should not exist at any point where one can
  assign one ODSI object to another.
*/

OsiDylpSolverInterface &ODSI::operator= (const OsiDylpSolverInterface &rhs)

{ if (this != &rhs)
  { destruct_problem(false) ;
  
    OSI::operator=(rhs) ;

    if (rhs.consys)
    { bool r UNUSED = consys_dupsys(rhs.consys,&consys,rhs.consys->parts) ;
      assert(r) ; }
    else
    { consys = 0 ; }
    if (rhs.lpprob)
    { lpprob = copy_lpprob(rhs.lpprob) ;
      lpprob->consys = consys ; }
    else
    { lpprob = 0 ; }
    solnIsFresh = rhs.solnIsFresh ;

    CLONE(lpopts_struct,rhs.initialSolveOptions,initialSolveOptions) ;
    CLONE(lpopts_struct,rhs.resolveOptions,resolveOptions) ;
    CLONE(lptols_struct,rhs.tolerances,tolerances) ;

    lp_retval = rhs.lp_retval ;
    obj_sense = rhs.obj_sense ;
    odsiInfinity = rhs.odsiInfinity ;
    mps_debug = rhs.mps_debug ;

    if (rhs.activeBasis)
    { activeBasis = rhs.activeBasis->clone() ;
      activeIsModified = rhs.activeIsModified ; }

    addedColCnt = rhs.addedColCnt ;
    addedRowCnt = rhs.addedRowCnt ;

    _objval = rhs._objval ;
    _col_obj = 0 ;
    _col_x = 0 ;
    _col_cbar = 0 ;
    _row_rhs = 0 ;
    _row_lower = 0 ;
    _row_upper = 0 ;
    _row_sense = 0 ;
    _row_range = 0 ;
    _row_price = 0 ;
    _row_lhs = 0 ;
    _matrix_by_col = 0 ;
    _matrix_by_row = 0 ;

    preObj_ = 0 ;
    postActions_ = 0 ;
    postObj_ = 0 ;
    passLimit_ = rhs.passLimit_ ;
    keepIntegers_ = rhs.keepIntegers_ ;
    savedConsys_ = 0 ;
    saved_col_obj = 0 ;
    saved_row_rhs = 0 ;
    saved_row_lower = 0 ;
    saved_row_upper = 0 ;
    saved_row_sense = 0 ;
    saved_row_range = 0 ;
    saved_matrix_by_col = 0 ;
    saved_matrix_by_row = 0 ;

    int n = getNumCols() ;
    int m = getNumRows() ;
    COPY_VEC(double,rhs._col_x,_col_x,n) ;
    COPY_VEC(double,rhs._row_price,_row_price,m) ;

    COPY_VEC(void *,&rhs.info_[0],&info_[0],OsiLastHintParam) ;

    reference_count++ ;

#   if ODSI_PARANOIA >= 2
    assert_same(*this, rhs, true) ;
#   endif
  }
# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << "): assign from "
    << &rhs << std::dec << "." << std::endl ;
# endif

  return (*this) ; }


//@}



/* OsiDylpSolverInterface destructor and related subroutines. */

/*! \defgroup DestructorHelpers Destructor Helpers */
//@{

/*! \brief Destroy LP problem. 

  This routine will free the problem-specific structures in an ODSI object.
  It relinquishes ownership of the solver, if necessary, and then the
  lpprob_struct, the main structure passed to dylp, is destroyed.
  Consys_free will free the constraint system, and dy_freesoln will free the
  remaining data structures associated with lpprob. There remains only to
  destroy the lpprob itself. To finish the job, free the options, tolerances,
  and statistics structures and call destruct_cache to free the cache
  structures.

  Because the lpprob structure is created only when dylp is actually
  called, it's possible that the interface contains only a constraint system.

  With the addition of presolve capability, it's also possible to have an
  lpprob but (apparently) no consys. This comes about if the client, for some
  reason, solves the lp and then decides to try again using presolve.
  Presolve will save the original constraint system and then restore it. Note
  that there's no point in saving any existing basis, even when applying
  presolve in resolve(), as we'll need a modified basis to get started.

  If the base ODSI object will be retained, set preserve_interface to true.
  This will retain the current options and tolerance settings.
*/

void ODSI::destruct_problem (bool preserve_interface)

{ 
/*
  If this object claims ownership of the solver, it should possess an lpprob
  structure created when the solver was called.
*/
  assert((dylp_owner != this) || (dylp_owner == this && lpprob)) ;

  if (dylp_owner == this)
  { detach_dylp() ; }

  if (lpprob)
  { assert(lpprob->consys == consys) ;
    if (consys)
    { consys_free(consys) ;
      consys = 0 ; }
    dy_freesoln(lpprob) ;
    delete lpprob ;
    lpprob = 0 ; }
  else
  if (consys)
  { consys_free(consys) ;
    consys = 0 ; }
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::destruct_problem." << std::endl ;
# endif
  addedColCnt = 0 ;
  addedRowCnt = 0 ;
 
  if (hotstart_fallback)
  { delete hotstart_fallback ;
    hotstart_fallback = 0 ; }
  if (activeBasis)
  { 
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::destruct_problem: deleting active basis "
      << activeBasis << std::dec << "." << std::endl ;
#   endif
    delete activeBasis ;
    activeBasis = 0 ;
    activeIsModified = false ; }

  destruct_cache(true,true) ;

  if (preserve_interface == false)
  { if (initialSolveOptions)
    { delete initialSolveOptions ;
      initialSolveOptions = 0 ; }
    if (resolveOptions)
    { delete resolveOptions ;
      resolveOptions = 0 ; }
    if (tolerances)
    { delete tolerances ;
      tolerances = 0 ; } }
# ifdef ODSI_STATISTICS
  if (statistics) dy_freestats(&statistics) ;
# endif
  
  return ; }


/*! \brief Detach an ODSI instance from the dylp solver

  Dylp retains static data structures from one call to the next, in
  anticipation of reoptimizing following problem modification. When an ODSI
  instance is destroyed or another instance needs to use dylp, the current
  instance must be detached. Setting lpctlONLYFREE and phase == dyDONE tells
  dylp the call is solely to free data structures.

  Note that we need to take care not to change the lpprob control flags
  (other than lpctlDYVALID, which will be cleared by dylp), as this could be
  a temporary detach and the lpprob will be used again by its owner.
*/

void ODSI::detach_dylp ()

{ assert(dylp_owner == this && lpprob && lpprob->consys) ;
  flags save_flags = getflg(lpprob->ctlopts,lpctlNOFREE|lpctlONLYFREE) ;
  clrflg(lpprob->ctlopts,lpctlNOFREE) ;
  setflg(lpprob->ctlopts,lpctlONLYFREE) ;
  lpprob->phase = dyDONE ;
/*
  Either of the initialSolve or resolve option blocks are ok here; dylp does
  not look.
*/
# ifdef ODSI_INFOMSGS
  CoinMessageHandler *hdl = messageHandler() ; 
  hdl->message(ODSI_DETACH,messages_)
    << (int) reinterpret_cast<CoinIntPtr>(this)
    << CoinMessageEol ;
# endif
  dylp(lpprob,initialSolveOptions,tolerances,statistics) ;
  clrflg(lpprob->ctlopts,lpctlONLYFREE) ;
  setflg(lpprob->ctlopts,save_flags) ;
  dylp_owner = 0 ; }

//@}

/*! \ingroup ODSIConstructorsDestructors
    \brief Destructor for OsiDylpSolverInterface

    Destruction of the last ODSI object triggers shutdown of the basis
    maintenance and i/o subsystems and deallocation of all dylp internal
    data structures.
*/

ODSI::~OsiDylpSolverInterface ()

{
/*
  Destroy the problem-specific structures used by dylp, and close the log
  and output files, if open.
*/
  destruct_presolve() ;
  destruct_problem(false) ;
  if (dyio_isactive(local_logchn)) (void) dyio_closefile(local_logchn) ;
  if (dyio_isactive(local_outchn)) (void) dyio_closefile(local_outchn) ;

  reference_count-- ;
  if (reference_count == 0)
  { if (basis_ready == true)
    { dy_freebasis() ;
      basis_ready = false ; }
    dyio_ioterm() ;
    errterm() ; }

# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << "): destructor." << std::endl ;
# endif

  return ; }

/*! \ingroup ODSIConstructorsDestructors
    \brief Reset the SI to the state produced by the default constructor.

    When it's inconvenient to destroy one SI and construct a new one, this
    function can be used to reset an existing SI to the same state produced
    by the default constructor.
*/

void ODSI::reset ()

{ 
/*
  Destroy the problem-specific structures used by dylp, and close the log
  file, if open.
*/
  destruct_presolve() ;
  destruct_problem(false) ;
  if (dyio_isactive(local_logchn))
  { (void) dyio_closefile(local_logchn) ;
    local_logchn = IOID_NOSTRM ; }
  if (dyio_isactive(local_outchn))
  { (void) dyio_closefile(local_outchn) ;
    local_outchn = IOID_NOSTRM ; }
/*
  That takes care of cleaning out the ODSI object. Call setInitialData to
  reset the underlying OSI object.
*/
  setInitialData() ;
/*
  Now the equivalent for ODSI: ensure default values for various fields in the
  object. Regrettably, messages_ is not a pointer, so we can't preserve the
  CoinMessages structure over the reset.
*/
  initial_gtxecho = false ;
  resolve_gtxecho = false ;
  lp_retval = lpINV ;
  setObjSense(1.0) ;
  mps_debug = false ;
  construct_options() ;
  setOsiDylpMessages(CoinMessages::us_en) ;
  for (int i = 0 ; i < OsiLastHintParam ; i++) info_[i] = 0 ;
  setHintParam(OsiDoPresolveInInitial,true,OsiForceDo,0) ;

# if ODSI_TRACK_SOLVERS > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << "): reset." << std::endl ;
# endif

  return ; }



/*! \defgroup ProbAdjust Functions for Problem Modification
    \brief Public functions for problem modification

  This group of functions comprise the public interface for adjusting the
  optimization problem seen by the dylp solver.
*/
//@{

inline void ODSI::setContinuous (int j)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(j,true,"setContinuous") ;
# endif

  if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }
/*
  Keep up with the bookkeeping.
*/
  vartyp_enum vtyp = consys->vtyp[idx(j)] ;
  switch (vtyp)
  { case vartypBIN:
    { consys->binvcnt-- ;
      break ; }
    case vartypINT:
    { consys->intvcnt-- ;
      break ; }
    default:
    { break ; } }
/*
  Set the new type. Changing the variable type should not affect the value of
  the current lp solution (which sees everything as continuous) hence there's
  no need to mark the solution as stale.
*/
  consys->vtyp[idx(j)] = vartypCON ;
  
  return ; }


inline void ODSI::setContinuous (const int* indices, int len)

{ for (int i = 0 ; i < len ; i++) setContinuous(indices[i]) ;

  return ; }
  

inline void ODSI::setInteger (int j)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(j,true,"setInteger") ;
# endif

  if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }
/*
  Keep up with the bookkeeping.
*/
  vartyp_enum vtyp = consys->vtyp[idx(j)] ;
  switch (vtyp)
  { case vartypBIN:
    { consys->binvcnt-- ;
      break ; }
    case vartypINT:
    { consys->intvcnt-- ;
      break ; }
    default:
    { break ; } }
/*
  Set the new type. Changing the variable type should not affect the value of
  the current lp solution (which sees everything as continuous) hence there's
  no need to mark the solution as stale.
*/
  if (getColLower()[j] == 0.0 && getColUpper()[j] == 1.0)
  { consys->vtyp[idx(j)] = vartypBIN ;
    consys->binvcnt++ ; }
  else
  { consys->vtyp[idx(j)] = vartypINT ;
    consys->intvcnt++ ; }

  return ; }


inline void ODSI::setInteger (const int* indices, int len)

{ for (int i = 0 ; i < len ; i++) setInteger(indices[i]) ;

  return ; }


/*!
  If the variable is binary, and the new bound is not 0 or 1, the type of the
  variable is automatically converted to general integer. If the variable is
  integer, and the new bound is not, the type of the variable is converted to
  continuous.

  Other conversions are too ambiguous --- suppose a general integer variable
  were temporarily bounded between 0 and 1. Should the type really be
  converted from general integer to binary? Similarly, bounding a continuous
  variable to an integer range should not cause it to become integer.
*/

inline void ODSI::setColLower (int i, double val)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(i,true,"setColLower") ;
# endif

  if (!consys->vlb)
  { bool r = consys_attach(consys,CONSYS_VLB,sizeof(double),
				  reinterpret_cast<void **>(&consys->vlb)) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }

  double primalTol ;
  (void) getDblParam(OsiPrimalTolerance,primalTol) ;

/*
  Make a clean integer value if the variable type is integer.
*/
  double cleanval ;
  if (isInteger(i))
  { cleanval = ceil(val-primalTol) ; }
  else
  { cleanval = val ; }
/*
  Change the bound. In general, this can result in a change in the optimal
  solution, but we'll be punctilious and only mark the solution as stale if the
  new bound conflicts with the primal solution value. But ... if the solution's
  already stale, don't check further.
*/
  consys->vlb[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlLBNDCHG) ;

  if (solnIsFresh == true)
  { const double *xvals = getColSolution() ;
    if (xvals[i] < val-primalTol)
    { solnIsFresh = false ;
      destruct_col_cache(false) ;
#     if ODSI_TRACK_FRESH > 0
      std::cout
	<< "ODSI(" << std::hex << this << std::dec
	<< ")::setColLower: new bound "
	<< val << " exceeds current value " << xvals[i]
	<< " by " << (val - xvals[i]) << "." << std::endl ;
#     endif
    } }

# if 0
  if (isInteger(i))
  { if (floor(val) != val)
      setContinuous(i) ;
    else
    if (isBinary(i) && !(val == 0.0 || val == 1.0))
      setInteger(i) ; }
# endif
}

/*!
  See the comments with setColLower re. automatic variable type conversions.
  In short, binary to general integer and integer to continuous are the only
  ones that will occur.
*/
inline void ODSI::setColUpper (int i, double val)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(i,true,"setColUpper") ;
# endif

  if (!consys->vub)
  { bool r = consys_attach(consys,CONSYS_VUB,sizeof(double),
			   reinterpret_cast<void **>(&consys->vub)) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }

  double primalTol ;
  (void) getDblParam(OsiPrimalTolerance,primalTol) ;

/*
  Make a clean integer value if the variable type is integer.
*/
  double cleanval ;
  if (isInteger(i))
  { cleanval = floor(val+primalTol) ; }
  else
  { cleanval = val ; }
/*
  Change the bound. In general, this can result in a change in the optimal
  solution, but we'll be punctilious and only mark the solution as stale if the
  new bound conflicts with the primal solution value.
*/
  consys->vub[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlUBNDCHG) ;

  if (solnIsFresh == true)
  { const double *xvals = getColSolution() ;
    if (xvals[i] > val+primalTol)
    { solnIsFresh = false ;
      destruct_col_cache(false) ;
#     if ODSI_TRACK_FRESH > 0
      std::cout
	<< "ODSI(" << std::hex << this << std::dec
	<< ")::setColUpper: new bound "
	<< val << " exceeds current value " << xvals[i]
	<< " by " << (xvals[i]-val) << "." << std::endl ;
#     endif
    } }

# if 0
  if (isInteger(i))
  { if (floor(val) != val)
      setContinuous(i) ;
    else
    if (isBinary(i) && !(val == 0.0 || val == 1.0))
      setInteger(i) ; }
# endif
}


/*!
  A call to this routine destroys all cached row values and any cached
  solution.
*/

void ODSI::setRowType (int i, char sense, double rhs, double range)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(i,false,"setRowType") ;
# endif
/*
  Install the change. In general, this can change the optimal solution.
*/
  int k = idx(i) ;
  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		sense,rhs,range) ;
  if (resolveOptions) resolveOptions->forcewarm = true ;
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::setRowType: new row type."
    << std::endl ;
# endif
/*
  Destroy cached values. We need to clear the structural side of the row
  cache.
*/
  destruct_row_cache(true) ;
  destruct_col_cache(false) ; }


/*!
  A call to this routine destroys all cached row values and any cached solution
  values.
*/

void ODSI::setRowUpper (int i, double val)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(i,false,"setRowUpper") ;
# endif

  int k = idx(i) ;
  double clbi = -odsiInfinity ;

  switch (consys->ctyp[k])
  { case contypEQ:
    case contypGE:
    { clbi = consys->rhs[k] ;
      break ; }
    case contypRNG:
    { clbi = consys->rhslow[k] ;
      break ; }
    case contypLE:
    case contypNB:
    { clbi = -odsiInfinity ;
      break ; }
    default:
    { assert(false) ; } }
/*
  Install the change. In general, this can change the optimal solution.
*/
  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		clbi,val) ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlRHSCHG) ;
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::setRowUpper: new bound."
    << std::endl ;
# endif

  destruct_row_cache(true) ;
  destruct_col_cache(false) ; }


/*!
  A call to this routine destroys all cached row values and any cached solution
  values.
*/

void ODSI::setRowLower (int i, double val)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(i,false,"setRowLower") ;
# endif

  int k = idx(i) ;
  double cubi ;
  contyp_enum ctypi = consys->ctyp[k] ;

  if (ctypi == contypGE || ctypi == contypNB)
    cubi = odsiInfinity ;
  else
    cubi = consys->rhs[k] ;
/*
  Install the change. In general, this can change the optimal solution.
*/
  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		val,cubi) ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlRHSCHG) ;
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::setRowLower: new bound."
    << std::endl ;
# endif

  destruct_row_cache(true) ;
  destruct_col_cache(false) ; }


/*!
  Add a row to the constraint system given the coefficients and upper and lower
  bounds on the left-hand-side.

  A call to this routine destroys all cached values. Unlike deleteRows,
  however, we don't need to fiddle with the basis, as dylp will automatically
  pick up the added constraints.
*/

inline void ODSI::addRow (const CoinPackedVectorBase &row,
			  double rlb, double rub)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,rlb,rub) ;
  add_row(row,'a',ctypi,rhsi,rhslowi) ;

  return ; }


/*!
  Add a row to the constraint system given the coefficients, constraint sense,
  right-hand-side value, and range.

  A call to this routine destroys all cached values. Unlike deleteRows,
  however, we don't need to fiddle with the basis, as dylp will automatically
  pick up the added constraints.
*/

inline void ODSI::addRow (const CoinPackedVectorBase &coin_row, 
			  char sense, double rhs, double range)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,sense,rhs,range) ;
  add_row(coin_row,'a',ctypi,rhsi,rhslowi) ;

  return ; }


/*!
  One cut at a time, please, when it comes to rows.

  A call to this routine destroys all cached values. Unlike deleteRows,
  however, we don't need to fiddle with the basis, as dylp will automatically
  pick up the added constraints.
*/

void ODSI::applyRowCut (const OsiRowCut &cut)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,cut.lb(),cut.ub()) ;
  add_row(cut.row(),'c',ctypi,rhsi,rhslowi) ;
  
  return ; }


/*!
  A call to this routine destroys all cached values.

  To make this work properly (in the absence of a more efficient routine
  implemented in the consys package), the trick is to sort the indices
  from smallest to largest, and delete in reverse order. It has the
  added advantage of being marginally more efficient.

  Note that we must also keep the active basis up to date, or remove it if we
  cannot guarantee validity.
*/

void ODSI::deleteRows (int count, const int* rows)

{ if (count <= 0) return ;
/*
  Sort the indices and do the deletions, one by one. This invalidates any
  existing solution in lpprob.
*/
  vector<int> lclrows = vector<int>(&rows[0],&rows[count]) ;

  if (count > 1) std::sort(lclrows.begin(),lclrows.end()) ;

  for (int k = count-1 ; k >= 0 ; k--)
  { int i = idx(lclrows[k]) ;
    bool r = consys_delrow_stable(consys,i) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::deleteRows: deleted " << count
    << " rows." << std::endl ;
# endif
/*
  Now, see if there's an active basis. If so, check that all the constraints
  to be deleted are slack. If they are, we can delete them from activeBasis
  and still guarantee a valid basis. If not, throw away activeBasis.
*/
  if (activeBasis)
  { bool allslack = true ;
    OsiDylpWarmStartBasis *odwsb =
      dynamic_cast<OsiDylpWarmStartBasis *>(activeBasis) ;
    for (int k = count-1 ; k >= 0 ; k--)
    { int i = lclrows[k] ;
      if (odwsb->getArtifStatus(i) != CoinWarmStartBasis::basic)
      { allslack = false ;
	break ; } }
    if (allslack == true)
    { odwsb->compressRows(count,rows) ;
      activeIsModified = true ;
      resolveOptions->forcewarm = true ; }
    else
    { 
#     if ODSI_TRACK_ACTIVE > 0
      std::cout
	<< "ODSI(" << std::hex << this
	<< ")::deleteRows: deleted tight constraints, deleting basis "
	<< activeBasis << std::dec
	<< "." << std::endl ;
#     endif
      delete activeBasis ;
      activeBasis = 0 ;
      activeIsModified = false ; } }

  destruct_cache(true,false) ; }



/*!
  Change a coefficient in the objective function. In general, this can change
  the optimal solution.
*/

void ODSI::setObjCoeff (int j, double objj)

{ 
# if ODSI_PARANOIA >= 1
  indexCheck(j,true,"setObjCoeff") ;
# endif
  
  consys->obj[idx(j)] = getObjSense()*objj ;
  if (_col_obj) _col_obj[j] = objj ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ;

  return ; }
  

/*!
  Change the sense of the objective; use 1.0 to request the objective be
  minimised, -1.0 to request it be maximised.
*/

void ODSI::setObjSense (double val)
/*
  The `natural' action of OSI (and dylp) is minimisation. Maximisation is
  accomplished as min -cx. So using -1 for maximisation is more natural
  than you'd think at first glance.

  Changing the objective sense will in general change the optimal solution.
*/

{ int n = getNumCols() ;

  if (n > 0 && val != obj_sense)
  { double *tmpobj = INV_VEC(double,consys->obj) ;
    std::transform(tmpobj,tmpobj+n,tmpobj,std::negate<double>()) ;
    if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ;
    solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::setObjSense: changing to "
    << ((val < 0)?"minimisation":"maximisation") << "." << std::endl ;
# endif
  }
  
  obj_sense = val ;
  
  return ; }


/*!
  Add a column to the constraint system given the coefficients, upper and lower
  bounds, and objective function coefficient. The variable type defaults to
  continuous.

  A call to this routine destroys all cached values.
*/

inline void ODSI::addCol (const CoinPackedVectorBase& coin_col,
			  double vlb, double vub, double obj)

{ add_col(coin_col,vartypCON,vlb,vub,obj) ; }

/*!
  Perhaps better named `applyColCuts', as an OsiColCut object can contain
  any number of adjustments to variable bounds. The routine simply steps
  through the new bounds, tightening existing bounds as necessary. Note that
  a bound is never loosened by this routine.
*/

void ODSI::applyColCut (const OsiColCut &cut)

{ const double* old_lower = getColLower() ;
  const double* old_upper = getColUpper() ;
  const CoinPackedVector& new_lower = cut.lbs() ;
  const CoinPackedVector& new_upper = cut.ubs() ;

  int i ;
  int n1 = new_lower.getNumElements() ;
  int n2 = new_upper.getNumElements() ;

  for (i = 0 ; i < n1 ; i++)
  { int j = new_lower.getIndices()[i] ;
    double x = new_lower.getElements()[i] ;
    if (x > old_lower[j]) setColLower(j,x) ; }

  for (i = 0 ; i < n2 ; i++)
  { int j = new_upper.getIndices()[i] ;
    double x = new_upper.getElements()[i] ;
    if (x < old_upper[j]) setColUpper(j,x) ; }
}


/*!
  A call to this routine destroys all cached values.

  To make this work properly (in the absence of a more efficient routine
  implemented in the consys package), the trick is to sort the indices
  from largest to smallest, and delete in the same order. It has the
  added advantage of being marginally more efficient.

  Note that we must also keep the active basis up to date, or remove it if we
  cannot guarantee validity.
*/

void ODSI::deleteCols (int count, const int* cols)

{ if (count <= 0) return ;
/*
  Sort the indices and do the deletions, one by one. This invalidates the
  current optimal solution.
*/
  vector<int> lclcols = vector<int>(&cols[0],&cols[count]) ;

  if (count > 1) std::sort(lclcols.begin(),lclcols.end()) ;

  for (int k = 0 ; k < count ; k++)
  { int j = idx(lclcols[k]) ;
    bool r = consys_delcol(consys, j) ;
    if (!r)
    { lp_retval = lpFATAL ;
      return ; } }
  solnIsFresh = false ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::deleteCols: deleted " << count
    << "columns." << std::endl ;
# endif
/*
  Now, see if there's an active basis. If so, check that all the variables to
  be deleted are nonbasic. If they are, we can delete them from activeBasis
  and still guarantee a valid basis. If not, throw away activeBasis.
*/
  if (activeBasis)
  { bool allnonbasic = true ;
    OsiDylpWarmStartBasis *odwsb =
      dynamic_cast<OsiDylpWarmStartBasis *>(activeBasis) ;
    for (int k = count-1 ; k >= 0 ; k--)
    { int j = lclcols[k] ;
      if (odwsb->getStructStatus(j) == CoinWarmStartBasis::basic)
      { allnonbasic = false ;
	break ; } }
    if (allnonbasic == true)
    { odwsb->deleteColumns(count,cols) ;
      activeIsModified = true ;
      resolveOptions->forcewarm = true ; }
    else
    { 
#     if ODSI_TRACK_ACTIVE > 0
      std::cout
	<< "ODSI(" << std::hex << this
	<< ")::deleteCols: deleted basic variables, deleting basis "
	<< activeBasis << std::dec
	<< "." << std::endl ;
#     endif
      delete activeBasis ;
      activeBasis = 0 ;
      activeIsModified = false ; } }

  destruct_cache(false,true) ; }


/*!
  This function will change the primal solution recorded in the interface,
  but has no effect on the actions of the solver. Its primary utility is to
  establish a set of values that can be used to calculate an objective value
  prior to invoking dylp. dylp is not prepared to accept a set of values for
  the primal variables. It firmly believes it should calculate these values
  given a basis, rhs, and bounds.
*/

void ODSI::setColSolution (const double* solution)

{ int n = getNumCols() ;
/*
  You can't set what you don't have. Otherwise, replace the current solution.
*/
  if (n == 0) return ;

  if (_col_x) delete[] _col_x ;
  _col_x = new double[n] ;
  
  COPY_VEC(double,solution,_col_x,n) ;
  calc_objval() ;

  return ; }


/*!
  This function will change the dual solution recorded in the interface, but
  has no effect on the actions of the solver. Its primary utility is to
  establish a set of values that can be used in calculation prior to invoking
  dylp. dylp is not prepared to accept a set of values for the dual
  variables. It firmly believes it should calculate these values given a
  basis and objective coefficients.
*/

void ODSI::setRowPrice (const double* price)

{ int m = getNumRows() ;

/*
  You can't set what you don't have. Otherwise, replace the current solution.
*/
  if (m == 0) return ;

  if (_row_price) delete[] _row_price ;
  _row_price = new double[m] ;
  
  COPY_VEC(double,price,_row_price,m) ;

  return ; }

//@} // ProbAdjust



#ifndef _MSC_VER

/*! \defgroup CopyVerifiers Copy Verification Routines
    \brief Private helper functions to verify correctness of copies

  This group of helper functions verify the correctness of copies of the
  data structures used in ODSI objects. For uniformity, all functions take
  a parameter, `exact', which influences the rigour of the test for some types.
*/
//@{


/*! \brief Compare two double values for equality.

  Exact bit-for-bit equality (i.e., d1 == d2) is usually not what is wanted
  when comparing floating point values. And even when it is, NaN doesn't
  equal NaN. Why is NaN not simply an error? It happens on occasion that we're
  cloning some structure and, for one reason or another, some values are
  essentially uninitialised. For example, we've added cuts but not yet
  reoptimised. The dual vector has been extended, but the values are junk.
  On occasion, dylp will deliberately initialise parts of a vector with NaN
  as a guard, and check later to see if it shows up in honest values.

  The inexact test looks for equality within a tolerance of 1.0e-10, scaled
  by the maximum of the two values.  To do a toleranced test, both values
  must be finite.
*/

void ODSI::assert_same (double d1, double d2, bool exact)

{ if (d1 == d2) return ;
  if (CoinIsnan(d1) && CoinIsnan(d2)) return ;

  assert(!exact && CoinFinite(d1) && CoinFinite(d2)) ;

  static const double epsilon UNUSED = 1.e-10 ;
  double tol UNUSED = std::max(fabs(d1),fabs(d2))+1 ;
  double diff UNUSED = fabs(d1 - d2) ;

  assert(diff <= tol*epsilon) ; }


/*! \brief Byte-wise equality comparison of a structure

  Byte-for-byte comparison of a structure.
*/

template<class T>
  void ODSI::assert_same (const T& t1, const T& t2, bool exact)


{ if (&t1 == &t2) return ;
  assert(memcmp(&t1,&t2,sizeof(T)) == 0) ; }


/*! \brief Element-wise comparison of two primitive arrays.  */

template<class T>
  void ODSI::assert_same (const T* t1, const T* t2, int n, bool exact)

{ if (t1 == t2) return ;
  for (int i=0 ; i<n ; i++) assert_same(t1[i],t2[i],exact) ; }


/*! \brief Verify copy of a dylp basis structure

  Check that two dylp basis structures are the same. They should have the
  same length and identical basis arrays. Basis_struct's use 1-based indexing,
  hence location 0 is unused.
*/

void ODSI::assert_same (const basis_struct& b1, const basis_struct& b2,
			bool exact)
{ assert(exact) ;                // be explicit
  if (&b1 == &b2) return ;
  assert(b1.len == b2.len) ;

  int size UNUSED = b1.len*sizeof(basisel_struct) ;
  assert(memcmp(inv_vec<basisel_struct>(b1.el),
		inv_vec<basisel_struct>(b2.el),size) == 0) ; }


/*! \brief Verify copy of a dylp constraint bound

  This structure contains upper and lower bounds on the value of the
  left-hand-side of a constraint, calculated from upper and lower bounds on the
  variables. Unused when dylp is operating solely as an lp solver.
*/

void ODSI::assert_same (const conbnd_struct& c1, const conbnd_struct& c2,
			bool exact)

{ assert(exact) ;                // be explicit
  if (&c1 == &c2) return ;

  //-- consys.h::conbnd_struct
  assert(c1.revs == c2.revs) ;
  assert(c1.inf == c2.inf) ;
  assert(c1.bnd == c2.bnd) ; }


/*! \brief Verify copy of a dylp constraint matrix (consys_struct)

  Compare two consys_struct's for equality. There are two levels of
  comparison. Exact looks for exact equivalence of content and size. Inexact
  checks that each consys_struct describes the same mathematical system, but
  ignores the constraint system name, various counts, and the allocated size
  of the structure.
*/

void ODSI::assert_same (const consys_struct& c1, const consys_struct& c2,
			bool exact)

{ if (&c1 == &c2) return ;
/*
  Simple fields: flags, counts, indices.
*/
  assert(c1.parts == c2.parts) ;
  assert(c1.opts == c2.opts) ;
  assert(c1.varcnt == c2.varcnt) ;
  assert(c1.archvcnt == c2.archvcnt) ;
  assert(c1.logvcnt == c2.logvcnt) ;
  assert(c1.intvcnt == c2.intvcnt) ;
  assert(c1.binvcnt == c2.binvcnt) ;
  assert(c1.maxcollen == c2.maxcollen) ;
  assert(!exact || (c1.maxcolndx == c2.maxcolndx)) ;
  assert(c1.concnt == c2.concnt) ;
  assert(c1.archccnt == c2.archccnt) ;
  assert(c1.cutccnt == c2.cutccnt) ;
  assert(c1.maxrowlen == c2.maxrowlen) ;
  assert(!exact || (c1.maxrowndx == c2.maxrowndx)) ;
  assert(c1.objndx == c2.objndx) ;
  assert(c1.xzndx == c2.xzndx) ;
/*
  Associated vectors.
*/
  int var_count = ODSI::idx(c1.varcnt) ;
  int con_count = ODSI::idx(c1.concnt) ;

  assert_same(c1.obj, c2.obj, var_count, exact) ;
  assert_same(c1.vtyp, c2.vtyp, var_count, true) ;
  assert_same(c1.vub, c2.vub, var_count, exact) ;
  assert_same(c1.vlb, c2.vlb, var_count, exact) ;
  assert_same(c1.colscale, c2.colscale, var_count, exact) ;

  assert_same(c1.rhs, c2.rhs, con_count, exact) ;
  assert_same(c1.rhslow, c2.rhslow, con_count, exact) ;
  assert_same(c1.ctyp, c2.ctyp, con_count, true) ;
  assert_same(c1.cub, c2.cub, con_count, exact) ;
  assert_same(c1.clb, c2.clb, con_count, exact) ;
  assert_same(c1.rowscale, c2.rowscale, con_count, exact) ;
/*
  These can differ and have no mathematical effect. Test them only if the
  client has requested exact equality.
*/
  if (exact)
  { assert(c1.nme == c2.nme) ;
    assert(c1.colsze == c2.colsze) ;
    assert(c1.rowsze == c2.rowsze) ;
    assert(c1.objnme == c2.objnme) ; } }


/*! \brief  Verify copy of a dylp lp problem (lpprob_struct)

  Essentially a matter of checking each individual component. Exact requires
  precise equivalence. Inexact allows for differences in the allocated capacity
  of the structures.
*/

void ODSI::assert_same (const lpprob_struct& l1, const lpprob_struct& l2,
			bool exact)

{ if (&l1 == &l2) return ;
/*
  Simple fields: flags, LP phase & return code, objective value, iteration
  count.
*/
  assert(l1.ctlopts == l2.ctlopts) ;
  assert(l1.phase == l2.phase) ;
  assert(l1.lpret == l2.lpret) ;
  assert_same(l1.obj,l2.obj,exact) ;
  assert(l1.iters == l2.iters) ;
/*
  The allocated capacity can differ unless the client requires exact
  equivalence.
*/
  if (exact)
  { assert(l1.colsze == l2.colsze) ;
    assert(l1.rowsze == l2.rowsze) ; }
/*
  Vectors: x, y, status
*/
  int min_col = idx(std::min(l1.colsze, l2.colsze)) ;
  int min_row = idx(std::min(l1.rowsze, l2.rowsze)) ;

  assert_same(l1.status,l2.status,min_col,exact) ;
  assert_same(l1.x,l2.x,min_row, exact) ;
  assert_same(l1.y,l2.y,min_row, exact) ;
/*
  Complex structures: the constraint system and basis.
*/
  assert_same(*l1.consys,*l2.consys,exact) ;
  assert_same(*l1.basis,*l2.basis,exact) ; }


/*! \brief Verify copy of an ODSI object.

  Verify equivalence by checking each component in turn.

  Note that statistics, hot start information, and cached data are never
  copied when an ODSI object is cloned or assigned, and open file descriptors
  are not inherited.
*/

void ODSI::assert_same (const OsiDylpSolverInterface& o1,
			const OsiDylpSolverInterface& o2, bool exact)

{ 
/*
  The same chunk of memory?
*/
  if (&o1 == &o2) return ;
/*
  These three fields are static. If they don't match, we're really
  confused.
*/
  assert(o1.reference_count == o2.reference_count) ;
  assert(o1.basis_ready == o2.basis_ready) ;
  assert(o1.dylp_owner == o2.dylp_owner) ;
/*
  Test the rest of the simple data fields.
*/
  assert(o1.initial_gtxecho == o2.initial_gtxecho) ;
  assert(o1.resolve_gtxecho == o2.resolve_gtxecho) ;
  assert(o1.lp_retval == o2.lp_retval) ;
  assert(o1.obj_sense == o2.obj_sense) ;
  assert(o1.odsiInfinity == o2.odsiInfinity) ;
  assert(o1.solvername == o2.solvername) ;
  assert(o1.mps_debug == o2.mps_debug) ;
/*
  Options and tolerances should be byte-for-byte identical.
*/
  assert_same(*o1.initialSolveOptions, *o2.initialSolveOptions, true) ;
  assert_same(*o1.resolveOptions, *o2.resolveOptions, true) ;
  assert_same(*o1.tolerances, *o2.tolerances, true) ;
/*
  The info_ array. The best we can do here is entry by entry equality.
*/
  assert_same(o1.info_,o2.info_,OsiLastHintParam,exact) ;
/*
  Now the more complicated structures. If the lpprob exists, it should
  reference a constraint system, and it should be the same constraint system
  as the consys field. The call to assert_same will check both the lpprob and
  the consys structures.
*/
  if (o1.lpprob)
  { assert(o2.lpprob) ;
    assert(o1.lpprob->consys && (o1.consys == o1.lpprob->consys)) ;
    assert(o2.lpprob->consys && (o2.consys == o2.lpprob->consys)) ;
    assert_same(*o1.lpprob,*o2.lpprob,exact) ; }
  else
  { assert(!o2.lpprob) ;
    assert_same(*o1.consys,*o2.consys,exact) ; }
/*
  A belt & suspenders check: retrieve the matrix as a COINPackedMatrix and
  use the COIN isEquivalent routine to check for equality.
*/
  const CoinPackedMatrix* m1 = o1.getMatrixByCol() ;
  const CoinPackedMatrix* m2 UNUSED = o2.getMatrixByCol() ;
  if (m1)
  { assert(m2 || m1->isEquivalent(*m2)) ; }
  else
  { assert(!m2) ; }
}

//@} // CopyVerifiers

#endif /* ! _MSC_VER */



/*
  Problem specification routines: MPS and control files, and loading a problem
  from OSI data structures.
*/

/*! \defgroup Main_ Options and Tolerances
    \brief Dylp Options and Tolerances Structures

  Dylp obtains options and tolerances through two global pointers. Changes to
  individual options and tolerances can be made by modifying values in the
  current options or tolerances structure. Wholesale changes can be made by
  redirecting main_lpopts or main_lptols to point to another structure.
*/
//@{

/*! \var lpopts_struct* main_lpopts
    \brief Points to the active dylp options structure
*/
/*! \var lptols_struct* main_lptols
    \brief Points to the active dylp tolerances structure
*/
  
#ifdef _MSC_VER
extern "C" lpopts_struct* main_lpopts ;
extern "C" lptols_struct* main_lptols ;
#endif
  
lpopts_struct* main_lpopts ;     // just for cmdint.c::process_cmds
lptols_struct* main_lptols ;     // just for cmdint.c::process_cmds

//@}

/*! \defgroup FileHelpers File I/O Helper Routines */
//@{

/*! \brief Construct a new name from filename, with ext1 removed and
	   ext2 added.

  A utility routine for building related file names. Starting with filename,
  ext1 is stripped, if present, and ext2 is added. Either of ext1 or ext2
  can be null or "" (the null string), in which case nothing is stripped or
  added, respectively.
*/

string ODSI::make_filename (const char *filename,
			    const char *ext1, const char *ext2)

{ string basename(filename) ;
  string ext1str(ext1), ext2str(ext2) ;

/*
  Extensions must start with ".", so fix if necessary.
*/
  if (ext1 && strlen(ext1) > 0 && ext1[0] != '.')
    ext1str.insert(ext1str.begin(),'.') ;
  if (ext2 && strlen(ext2) > 0 && ext2[0] != '.')
    ext2str.insert(ext2str.begin(),'.') ;
/*
  Strip ext1 from filename, if it exists, leaving the base name.
*/
  if (ext1 && strlen(ext1) > 0)
  { string tmpname(filename) ;
    string::size_type ext1pos = tmpname.rfind(ext1str) ;
    if (ext1pos < tmpname.npos)
    { basename = tmpname.substr(0,ext1pos) ; } }
/*
  Add ext2, if it exists.
*/
  if (ext2 && strlen(ext2) > 0) basename += ext2str ;

  return (basename) ; }


//@}



/*! \defgroup ODSILoadProb Methods to Load a Problem

   \brief Methods to load a problem into dylp

   This group of functions provides methods to load a problem from an MPS
   file, from an OSI packed matrix structure, or from a standard column-major
   packed matrix structure. A call to any of these functions results in the
   destruction of the existing problem and the creation of a new problem.

   All functions will preserve existing options and tolerances across the
   problem change. In general, if you want to reset options and tolerances
   to defaults between problems, you're best off to destroy the current ODSI
   object and create a new one. The additional overhead is small.
   
   loadProblem and assignProblem consider all variables to be continuous.
   assignProblem destroys its arguments; loadProblem leaves them unchanged.

   The readMps and writeMps routines extend the default OSI
   behaviour by allowing for multiple extensions on file names. They apply
   the following rule: any specified extension is added to the specified
   filename, with the exception that if the extension already exists, it is
   not duplicated. For example, if the extension is "mps", the file name
   myproblem.mps would be unchanged, while the filename myproblem.ver1 would
   become myproblem.ver1.mps.

   The extension will default to "mps" unless explicitly set to null (either
   the empty string, "", or 0).
*/
//@{



/*!
  Read a problem definition in MPS format.

  readMps will recognize continuous, binary, and general integer variables.

  If compiled with #ODSI_IMPLICIT_SPC defined as `1', before reading the MPS
  problem file, readMps looks for a dylp options file (extension ".spc"). The
  name is constructed by stripping the extension (if present) from the given
  file name and adding the extension ".spc".

  A `worst case' (and typically infeasible) primal solution is constructed
  which provides a valid bound on the objective.
*/

int ODSI::readMps (const char* basename, const char* ext)

/*
  This routine uses the COIN MPS input code, which allows me to sidestep
  the issues that dylp's MPS reader has when confronted with an MPS file
  that requires fixed-field processing or contains a constant objective
  function offset.

  Parameters:
    basename:	base name for mps file
    ext:	file extension (no leading `.')

  Returns: -1 if the MPS file could not be opened, otherwise the number of
	   errors encountered while reading the file.
*/

{ int errcnt ;
  CoinMpsIO mps ;
  CoinMessageHandler *mps_handler = mps.messageHandler() ;
  string filename ;

  if (mps_debug)
  { mps_handler->setLogLevel(handler_->logLevel()) ; }
  else
  { mps_handler->setLogLevel(0) ; }

# if ODSI_IMPLICIT_SPC == 1
/*
  See if there's an associated .spc file and read it first, if present.
*/
  filename = make_filename(basename,ext,"spc") ;
  dylp_controlfile(filename.c_str(),true,false) ;
# endif
/*
  Make sure the MPS reader has the same idea of infinity as dylp, then
  attempt to read the MPS file.
*/
  mps.setInfinity(odsiInfinity) ;
  filename = make_filename(basename,ext,ext) ;
  errcnt = mps.readMps(filename.c_str(),0) ;
  handler_->message(ODSI_MPSFILEIO,messages_)
    << filename << "read" << errcnt << CoinMessageEol ;
  if (errcnt != 0) return (errcnt) ;
/*
  Load the problem and build a worst-case solution. To take full advantage
  of information in the MPS file, it's easiest to pass the mps object to
  load_problem.
*/
  load_problem(mps) ;
/*
  Done!
*/
  return (0) ; }

/*!
  Read a problem definition in MPS format, including SOS information.

  Functionality is as for readMps(const char*,const char*), with the added
  feature that any special ordered sets (SOS) defined in the MPS file are
  returned through the \p numberSets and \p set parameters.
*/

int ODSI::readMps (const char* basename, const char* ext,
		   int &numberSets, CoinSet **&sets)

/*
  This routine uses the COIN MPS input code, which allows me to sidestep
  the issues that dylp's MPS reader has when confronted with an MPS file
  that requires fixed-field processing or contains a constant objective
  function offset. Note the extension in parameters to return SOS sets
  defined in the MPS file.

  Parameters:
    basename:	base name for mps file
    ext:	file extension (no leading `.')
    numberSets:	the number of SOS sets
    sets:	pointer to array of SOS sets

  Returns: -1 if the MPS file could not be opened, otherwise the number of
	   errors encountered while reading the file.
*/

{ int errcnt ;
  CoinMpsIO mps ;
  CoinMessageHandler *mps_handler = mps.messageHandler() ;
  string filename ;

  if (mps_debug)
  { mps_handler->setLogLevel(handler_->logLevel()) ; }
  else
  { mps_handler->setLogLevel(0) ; }

# if ODSI_IMPLICIT_SPC == 1
/*
  See if there's an associated .spc file and read it first, if present.
*/
  filename = make_filename(basename,ext,"spc") ;
  dylp_controlfile(filename.c_str(),true,false) ;
# endif
/*
  Make sure the MPS reader has the same idea of infinity as dylp, then
  attempt to read the MPS file.
*/
  mps.setInfinity(odsiInfinity) ;
  filename = make_filename(basename,ext,ext) ;
  errcnt = mps.readMps(filename.c_str(),0,numberSets,sets) ;
  handler_->message(ODSI_MPSFILEIO,messages_)
    << filename << "read" << errcnt << CoinMessageEol ;
  if (errcnt != 0) return (errcnt) ;
/*
  Load the problem and build a worst-case solution. To take full advantage
  of information in the MPS file, it's easiest to pass the mps object to
  load_problem.
*/
  load_problem(mps) ;
/*
  Done!
*/
  return (0) ; }

/*!
  Write a problem to the file basename.ext in MPS format.

  `sense' controls the objective sense in the MPS file. -1.0 forces
  maximisation, 1.0 forces minimisation. 0.0 leaves it up to the solver.
*/

void ODSI::writeMps (const char *basename,
		     const char *ext, double objsense) const

{ string filename = make_filename(basename,ext,ext) ;
  
  CoinMpsIO mps ;
  CoinMessageHandler *mps_handler = mps.messageHandler() ;

  if (mps_debug)
  { mps_handler->setLogLevel(handler_->logLevel()) ; }
  else
  { mps_handler->setLogLevel(0) ; }

  double *outputobj ;
  const double *const solverobj = getObjCoefficients() ;
  
  int n = getNumCols(),
      m = getNumRows() ;
/*
  If the client has no opinion about the sense of the objective, go with the
  solver's notion. If the client's desired sense conflicts with the solver's
  sense, negate the objective.
*/
  if (objsense == 0) objsense = getObjSense() ;
  if (objsense == getObjSense())
  { outputobj = const_cast<double *>(solverobj) ; }
  else
  { outputobj = new double[n] ;
    std::transform(solverobj,solverobj+n,outputobj,std::negate<double>()) ; }

  mps.setProblemName(consys->nme) ;

  char *vartyp = new char[n] ;
  typedef const char *charp ;
  const char **colnames = new charp[n],
	     **rownames = new charp[m] ;
  int i,j ;

  for (j = 0 ; j < n ; j++) vartyp[j] = isInteger(j) ;

  for (i = 0 ; i < m ; i++)
    rownames[i] = consys_nme(consys,'c',idx(i),false,0) ;
  
  for (j = 0 ; j < n ; j++)
    colnames[j] = consys_nme(consys,'v',idx(j),false,0) ;

  mps.setMpsData(*getMatrixByRow(),odsiInfinity,
		 getColLower(),getColUpper(),outputobj,vartyp,
		 getRowLower(),getRowUpper(),colnames,rownames) ;
/*
  We really need to work on symbolic names for these magic numbers.
*/
  int errcnt = mps.writeMps(filename.c_str(),0,0,2) ;
  handler_->message(ODSI_MPSFILEIO,messages_)
    << filename << "written" << errcnt << CoinMessageEol ;
/*
  Free up the arrays we allocated.
*/
  delete[] vartyp ;
  delete[] colnames ;
  delete[] rownames ;
  if (outputobj != solverobj) delete[] outputobj ;

  return ; }


/*!
  Constraints are described using an COIN packed matrix and upper and lower
  bounds for the left-hand-side.  The routine's parameters are destroyed once
  they have been copied to ODSI structures. Existing options and tolerances
  are preserved if they exist.
*/

inline void
ODSI::assignProblem (CoinPackedMatrix*& matrix,
		     double*& col_lower, double*& col_upper, double*& obj, 
		     double*& row_lower, double*& row_upper)

{ loadProblem(*matrix,col_lower,col_upper,obj,row_lower,row_upper) ;
  delete matrix ; matrix = 0 ;
  delete [] col_lower ; col_lower = 0 ;
  delete [] col_upper ; col_upper = 0 ;
  delete [] obj ; obj = 0 ;
  delete [] row_lower ; row_lower = 0 ;
  delete [] row_upper ; row_upper = 0 ;
}


/*!
  Constraints are described using a COIN packed matrix, constraint sense,
  right-hand-side, and (optional) range. The routine's parameters are
  destroyed once they have been copied to ODSI structures. Existing options
  and tolerances are preserved.
*/

inline void
ODSI::assignProblem (CoinPackedMatrix*& matrix, 
		     double*& lower, double*& upper, double*& obj, 
		     char*& sense, double*& rhs, double*& range)

{ loadProblem(*matrix, lower, upper, obj, sense, rhs, range) ;
  delete matrix ; matrix = 0 ;
  delete [] lower ; lower = 0 ;
  delete [] upper ; upper = 0 ;
  delete [] obj ; obj = 0 ;
  delete [] sense ; sense = 0 ;
  delete [] rhs ; rhs = 0 ;
  delete [] range ; range = 0 ;
}


/*!
  Constraints are described using an OSI packed matrix and upper and lower
  bounds for the left-hand-side. Existing options and tolerances are
  preserved if they exist; defaults are used otherwise.
*/

inline void
ODSI::loadProblem (const CoinPackedMatrix& matrix,
		   const double* col_lower, const double* col_upper,
		   const double* obj,
		   const double* row_lower, const double* row_upper)

{ int rowcnt = matrix.getNumRows() ;
  double *rhs = new double[rowcnt] ;
  double *rhslow = new double[rowcnt] ;
  contyp_enum *ctyp = new contyp_enum[rowcnt] ;

  gen_rowparms(rowcnt,rhs,rhslow,ctyp,row_lower,row_upper) ;
  load_problem(matrix,col_lower,col_upper,obj,ctyp,rhs,rhslow) ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ; }


/*!
  Constraints are described using an OSI packed matrix, constraint sense,
  right-hand-side, and (optional) range. Existing options and tolerances are
  preserved if they exist; defaults are used otherwise.
*/

inline void
ODSI::loadProblem (const CoinPackedMatrix& matrix, 
		   const double* col_lower, const double* col_upper,
		   const double* obj,
		   const char* sense, const double* rhsin, const double* range)

{ int rowcnt = matrix.getNumRows() ;
  double *rhs = new double[rowcnt] ;
  double *rhslow = new double[rowcnt] ;
  contyp_enum *ctyp = new contyp_enum[rowcnt] ;

  gen_rowparms(rowcnt,rhs,rhslow,ctyp,sense,rhsin,range) ;
  load_problem(matrix,col_lower,col_upper,obj,ctyp,rhs,rhslow) ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ; }

/*!
  Constraints are described using a standard column-major packed matrix
  structure and upper and lower bounds for the left-hand-side. Existing
  options and tolerances are preserved if they exist; defaults are used
  otherwise.
*/

inline void
ODSI::loadProblem (const int colcnt, const int rowcnt,
		   const int *start, const int *index, const double *value,
		   const double* col_lower, const double* col_upper,
		   const double* obj,
		   const double* row_lower, const double* row_upper)

{ double *rhs = new double[rowcnt] ;
  double *rhslow = new double[rowcnt] ;
  contyp_enum *ctyp = new contyp_enum[rowcnt] ;

  gen_rowparms(rowcnt,rhs,rhslow,ctyp,row_lower,row_upper) ;

  load_problem(colcnt,rowcnt,start,0,index,value,
	       col_lower,col_upper,obj,ctyp,rhs,rhslow) ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ; }


/*!
  Constraints are described using a standard column-major packed matrix
  structure, constraint sense, right-hand-side, and (optional) range.
  Existing options and tolerances are preserved if they exist; defaults are
  used otherwise.
*/

inline void
ODSI::loadProblem (const int colcnt, const int rowcnt,
		   const int *start, const int *index, const double *value,
		   const double* col_lower, const double* col_upper,
		   const double* obj,
		   const char* sense, const double* rhsin, const double* range)

{ double *rhs = new double[rowcnt] ;
  double *rhslow = new double[rowcnt] ;
  contyp_enum *ctyp = new contyp_enum[rowcnt] ;

  gen_rowparms(rowcnt,rhs,rhslow,ctyp,sense,rhsin,range) ;
  load_problem(colcnt,rowcnt,start,0,index,value,
	       col_lower,col_upper,obj,ctyp,rhs,rhslow) ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ; }

//@} // ODSILoadProb





/*! \defgroup SolverParms Methods to Set/Get Solver Parameters
    \brief Methods to set and retrive solver parameters.

  Dylp supports a limited set of the OSI parameters. Specifically,
  <dl>
    <dt>OsiDualTolerance</dt>
    <dd>(supported)
	This doesn't really correspond to a dylp parameter, but it can be
	faked (sort of).  Dylp's dual feasibility tolerance is dynamically
	scaled from an absolute dual zero tolerance. There is also an
	additional scaling factor (`lpcontrol dfeas' in the dylp
	documentation) that can be controlled by the user.  Given a value for
	the OsiDualTolerance parameter, the code calculates the appropriate
	value for dfeas such that OsiDualTolerance = dfeas*(absolute dual
	zero tolerance).
    </dd>

    <dt>OsiDualObjectiveLimit</dt>
    <dd>(facade)
	The dynamic simplex algorithm used in dylp works with a partial
	constraint system, adding and deleting constraints and variables
	as required and generally maintaining a minimal active constraint
	system.
	A side effect is that neither the dual or primal objective changes
	monotonically, hence it's not in general possible to terminate simplex
	iterations based on a limit on the dual objective.
	Currently, isDualObjectiveLimitReached simply bases its answer on the
	objective returned by dylp.
	This can (and will) be improved.
    </dd>

    <dt>OsiMaxNumIteration</dt>
    <dd>(supported)
	Limit on the number of simplex iterations.
    </dd>

    <dt>OsiMaxNumIterationHotStart</dt>
    <dd>(supported)
	Limit on the number of simplex iterations for an
	lp initiated by the solveFromHotStart function.
    </dd>

    <dt>OsiObjOffset</dt>
    <dd>(supported)
	A constant value added to the objective returned by the solver.
    </dd>
 
    <dt>OsiPrimalObjectiveLimit</dt>
    <dd>(facade)
	As for OsiDualObjectiveLimit.
    </dd>

    <dt>OsiPrimalTolerance</dt>
    <dd>(supported) Handled in the same manner as OsiDualTolerance.</dd>
  </dl>
*/

//@{

inline double ODSI::getInfinity () const { return odsiInfinity ; }


bool ODSI::setDblParam (OsiDblParam key, double value)
/*
  Aside from simply setting the relevant ODSI field, the work here is making
  sure the dylp tolerance values are consistent.

  OsiDualTolerance and OsiPrimalTolerance are both feasibility tolerances.
  dylp actually maintains feasibility tolerances scaled off the respective
  zero tolerances. dfeas_scale and pfeas_scale allow the scaled tolerances to
  be relaxed (or tightened, for that matter) with respect to the zero
  tolerances. Hence the best approach to consistency is to set *_scale based
  on the ratio of osi and dylp tolerances. This is, at best, imperfect, as it
  doesn't take into account dylp's internal scaling. Nor can it, really, as
  this information isn't guaranteed to be valid.

  Support for Osi[Dual,Primal]ObjectiveLimit is a facade at the moment. The
  values are stored, and used by is[Dual,Primal]ObjectiveLimitReached, but
  dylp knows nothing of them.
*/
{ if (key == OsiLastDblParam) return (false) ;

  bool retval = OSI::setDblParam(key,value) ;

  switch (key)
  { case OsiDualTolerance:
    { tolerances->dfeas_scale = value/tolerances->cost ;
      break ; }
    case OsiPrimalTolerance:
    { tolerances->pfeas_scale = value/tolerances->zero ;
      break ; }
    case OsiObjOffset:
    { break ; }
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    { break ; }
    default:
    { retval = false ;
      break ; } }
  
  return (retval) ; }


bool ODSI::getDblParam (OsiDblParam key, double& value) const
/*
  This simply duplicates OSI::getDblParam. I've kept it here until I decide
  what to do with the parameters Osi[Primal,Dual]ObjectiveLimit.
*/
{ if (key == OsiLastDblParam) return (false) ;

  bool retval ;

  switch (key)
  { case OsiDualTolerance:
    case OsiPrimalTolerance:
    case OsiObjOffset:
    { retval = OSI::getDblParam(key,value) ;
      break ; }
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    { retval = OSI::getDblParam(key,value) ;
      break ; }
    default:
    { OSI::getDblParam(key,value) ;
      retval = false ;
      break ; } }

  return (retval) ; }


bool ODSI::setIntParam (OsiIntParam key, int value)
/*
  Dylp's iteration limit is normally enforced on a per-phase basis, with an
  overall limit of 3*(per-phase limit). OsiMaxNumIterationHotStart is handled
  in solveFromHotStart.
*/
{ if (key == OsiLastIntParam) return (false) ;

  bool retval = OSI::setIntParam(key,value) ;

  switch (key)
  { case OsiMaxNumIteration:
    { initialSolveOptions->iterlim = value/3 ;
      resolveOptions->iterlim = initialSolveOptions->iterlim ;
      break ; }
    case OsiMaxNumIterationHotStart:
    { break ; }
    default:
    { retval = false ;
      break ; } }
  
  return (retval) ; }


inline bool ODSI::getIntParam (OsiIntParam key, int& value) const
/*
  This simply duplicates OSI::getIntParam. Retained here in anticipation of
  future extensions.
*/
{ bool retval ;

  if (key == OsiLastIntParam) return (false) ;

  switch (key)
  { case OsiMaxNumIteration:
    case OsiMaxNumIterationHotStart:
    { retval = OSI::getIntParam(key,value) ;
      break ; }
    default:
    { retval = false ;
      break ; } }
  
  return (retval) ; }


bool ODSI::setStrParam (OsiStrParam key, const string& value)
/*
  Set string paramaters. Note that OsiSolverName is, by definition, a
  constant. You can set it all you like, however :-).

  The problem name is pushed to the constraint system, if it exists.
*/

{ if (key == OsiLastStrParam) return (false) ;

  bool retval = OSI::setStrParam(key,value) ;

  switch (key)
  { case OsiProbName:
    { if (consys)
      { if (consys->nme) strfree(consys->nme) ;
	consys->nme = stralloc(value.c_str()) ; }
      break ; }
    case OsiSolverName:
    { break ; }
    default:
    { retval = false ;
      break ; } }
  
  return (retval) ; }

bool ODSI::getStrParam (OsiStrParam key, string& value) const
/*
  This mostly duplicates OSI::getStrParam. It enforces the notion that you
  can't set the solver name, which the default OSI routine does not.
*/
{ if (key == OsiLastStrParam) return (false) ;

  bool retval = true ;

  switch (key)
  { case OsiProbName:
    { retval = OSI::getStrParam(key,value) ;
      break ; }
    case OsiSolverName:
    { value = solvername ;
      break ; }
    default:
    { retval = false ;
      break ; } }

  return (retval) ; }


/*!
  A helper routine to deal with hints where dylp lacks any flexibility ---
  the facility is unimplemented, or always on. Failure is defined as
  OsiForceDo in the wrong direction; in this case an error is thrown. If the
  routine returns, then either the hint is compatible with dylp's capabilities,
  or it can be ignored.
*/

void ODSI::unimp_hint (bool dylpSense, bool hintSense,
		       OsiHintStrength hintStrength, const char *msgString)

{ if (dylpSense != hintSense)
  { string message("Dylp ") ;
    if (dylpSense == true)
    { message += "cannot disable " ; }
    else
    { message += "does not support " ; }
    message += msgString ;
    if (hintStrength == OsiForceDo)
    { handler_->message(ODSI_UNSUPFORCEDO,messages_)
        << message << CoinMessageEol ;
      throw CoinError(message,"setHintParam","OsiDylpSolverInterface") ; }
    else
    { handler_->message(ODSI_IGNOREDHINT,messages_)
	<< message << CoinMessageEol ; } }
  
  return ; }


/*!
  Dylp only looks at the options and tolerances structures for the solver
  instance --- it has no knowledge of hints. For hints that aren't
  implemented in the ODSI code, the approach is to stash them away and make
  the appropriate modifications in the options and tolerances structures.
  Hints are therefore persistent across all calls to the solver, until
  explicitly changed.

  Note that the client retains ownership of the blocks specified by the (void
  *) arguments. If ODSI took ownership, it would need to know how to allocate
  and free the block, and that way lies madness.
*/

bool ODSI::setHintParam (OsiHintParam key, bool sense,
			 OsiHintStrength strength, void *info)
/*
  OSI provides arrays for the sense and strength. ODSI provides the array for
  info.
*/

{ bool retval = false ;
/*
  Set the hint in the OSI structures. Unlike the other set*Param routines,
  setHintParam will return false for key == OsiLastHintParam. Unfortunately,
  it'll also throw for strength = OsiForceDo, without setting a return value.
  We need to catch that throw.
*/
  try
  { retval = OSI::setHintParam(key,sense,strength) ; }
  catch (CoinError)
  { retval = (strength == OsiForceDo) ; }
    
  if (retval == false) return (false) ;
  info_[key] = info ;
/*
  Did the client say `ignore this'? Who are we to argue.
*/
  if (strength == OsiHintIgnore) return (true) ;
/*
  We have a valid hint which would be impolite to simply ignore. Deal with
  it as best we can.
*/
  switch (key)
  {
/*
  In the process of implementing native presolve. initialSolve only at first.
  Nothing more to do here, as presolve is handled within ODSI.
*/
    case OsiDoPresolveInInitial:
    { retval = true ;
      break ; }
    case OsiDoPresolveInResolve:
    { unimp_hint(false,sense,strength,"presolve for resolve") ;
      retval = true ;
      break ; }
/*
  Dylp can suppress the dual, but cannot suppress the primal. Requesting
  OsiDoDual* with (true, OsiHint*) is ignored; (true,OsiForceDo) will throw
  an exception. On the other hand, suppressing the dual is generally not a
  good idea, so the hint is ignored at (false,OsiHintTry).
*/
    case OsiDoDualInInitial:
    case OsiDoDualInResolve:
    { unimp_hint(false,sense,strength,"exclusive use of dual simplex") ;
      retval = true ;
      lpopts_struct *opts =
	(key == OsiDoDualInInitial)?initialSolveOptions:resolveOptions ;
      if (sense == false)
      { if (strength >= OsiHintDo) opts->usedual = false ; }
      else
      { opts->usedual = true ; }
      break ; }
/*
  No issues here, dylp can go either way.
*/
    case OsiDoScale:
    { if (sense == false)
      { if (strength >= OsiHintTry) initialSolveOptions->scaling = 0 ; }
      else
      { initialSolveOptions->scaling = 2 ; }
      resolveOptions->scaling = initialSolveOptions->scaling ;
      break ; }
/*
  Dylp can construct a basis using only slacks and artificials (ibLOGICAL),
  or it can prefer architecturals over artificials (ibSLACK) or prefer
  architecturals over both slacks and artificials (ibARCH). The default is
  ibLOGICAL, so nothing happens if the hint sense is false. For true, change
  to ibSLACK, then ibARCH depending on strength.
*/
    case OsiDoCrash:
    { if (sense == true)
      { if (strength >= OsiHintDo)
	{ initialSolveOptions->coldbasis = ibARCH ; }
	else
	{ initialSolveOptions->coldbasis = ibSLACK ; } }
      else
      { initialSolveOptions->coldbasis = ibLOGICAL ; }
      retval = true ;
      break ; }
/*
  What to do with OsiDoInBranchAndCut? Try a context option for dylp, and see
  how it goes. If sense = true, BANDC is the right context for resolve,
  INITIALLP for initialSolve. If sense is false, go for INITIALLP for low
  strength, SINGLELP at ForceDo. In theory, resolve() should never see
  SINGLELP.
*/
    case OsiDoInBranchAndCut:
    { if (sense == true)
      { resolveOptions->context = cxBANDC ;
	initialSolveOptions->context = cxINITIALLP ; }
      else
      { if (strength < OsiForceDo)
	{ resolveOptions->context = cxINITIALLP ;
	  initialSolveOptions->context = cxINITIALLP ; }
	else
	{ resolveOptions->context = cxSINGLELP ;
	  initialSolveOptions->context = cxSINGLELP ; } }
      retval = true ;
      break ; }
/*
  The current absolute print settings are held in info_[OsiDoReducePrint].
  An unadorned hint changes the level by +/- 1, 2, or 3, depending on sense
  and strength (abusing the equivalence of enums and integers).  If info is
  non-zero, it is taken as a pointer to an integer which is the new absolute
  value for the print settings. In this case, sense is irrelevant.
  
  Note that CoinMessageHandler recognizes levels 0 -- 4 as generic levels.
  Conveniently, dylp has a similar set of generic levels. The larger powers
  of 2 (8, 16, 32, etc.) are left to trigger specific classes of debugging
  printout. These must be set using info. Currently:
    0x08	CoinMpsIO messages
    0x10	Echo to terminal
*/
    case OsiDoReducePrint:
    { int verbosity = reinterpret_cast<long>(info_[key]) ;
      mps_debug = false ;
      if (info)
      { verbosity = *reinterpret_cast<int *>(info) ;
	if (verbosity&0x8) mps_debug = true ; }
      else
      { if (sense == true)
	{ verbosity -= strength ;
	  verbosity = max(verbosity,0) ; }
	else
	{ verbosity += strength ;
	  verbosity = min(verbosity,4) ; } }
      info_[key] = reinterpret_cast<void *>(verbosity) ;

      dy_setprintopts(0,initialSolveOptions) ;
      dy_setprintopts(0,resolveOptions) ;
      if (verbosity&0x10)
      { initial_gtxecho = true ;
	resolve_gtxecho = true ; }
      else
      { initial_gtxecho = false ;
	resolve_gtxecho = false ; }
      messageHandler()->setLogLevel(verbosity) ;
      dy_setprintopts((verbosity&0x7),initialSolveOptions) ;
      dy_setprintopts((verbosity&0x7),resolveOptions) ;
      retval = true ;
      break ; }
/*
  The OSI spec says that unimplemented options (and, by implication, hints)
  should return false. In the case of a hint, however, we can ignore anything
  except OsiForceDo, so usability says we should anticipate new hints and set
  this up so the solver doesn't break. So return true.
*/
    default:
    { unimp_hint(!sense,sense,strength,"unrecognized hint") ;
      retval = true ;
      break ; } }

  return (retval) ; }


inline bool ODSI::getHintParam (OsiHintParam key, bool &sense,
			        OsiHintStrength &strength, void *&info) const
/*
  OSI provides the arrays for the sense and strength. ODSI provides the array
  for info.
*/

{ bool retval = OSI::getHintParam(key,sense,strength) ;

  if (retval == false) return (false) ;

  info = &info_[key] ;

  return (true) ; }


//@} // SolverParms



/*! \defgroup LPCommonMethods Common Solver Invocation Methods
    \brief Common core methods for invoking dylp
*/
//@{

lpret_enum ODSI::do_lp (ODSI_start_enum start)

/*
  This routine handles the job of running an lp and making a few attempts at
  recovery from errors. Ultimately, it's looking for one of lpOPTIMAL,
  lpINFEAS, or lpUNBOUNDED. In the event we turn up anything else, the
  strategy is to keep trying, forcing a cold start and gradually increasing
  the refactorisation frequency.

  dylp expects that >= constraint have been converted to <= constraints
  before it ever sees the system. So we do it here, and then flip back
  afterwards, because OSI (particularly the unitTest) isn't tolerant of the
  solver rearranging the constraint system.

  Note that the routines which dump the solution and statistics in human
  readable form also make use of the constraint system. They won't be put off
  by converting between >= and <=, but they'll be very lost if the constraint
  system they see has a different number of constraints or variables than the
  constraint system dylp was using.

  Parameter:
    start:	specifies cold, warm, or hot start
  
  Returns: whatever it gets from dylp
*/

{ int retries,ndx,cndx,flips ;
  basis_struct *basis ;
  lpret_enum lpret ;
  lpopts_struct lcl_opts ;
  lptols_struct lcl_tols ;
  flags persistent_flags ;
  const char *rtnnme = "do_lp" ;

  bool *flipped ;

# ifdef ODSI_INFOMSGS
  int print = 1 ;

  CoinMessageHandler *hdl = messageHandler() ; 
  hdl->message(ODSI_ALLDYLP,messages_)
    << startString(start) << (int) reinterpret_cast<CoinIntPtr>(this)
    << CoinMessageEol ;
# endif

/*
  Check for constraint system corruption, and set lp_retval to lpFATAL if
  we find it. If we currently own the solver, release it.
*/
  if (flgon(consys->opts,CONSYS_CORRUPT))
  { if (dylp_owner == this)
    { detach_dylp() ; }
    return (lpFATAL) ; }
/*
  Set up options and tolerances, make sure the phase isn't dyDONE, and
  reinitialise the statistics if we're collecting them.
*/
  lcl_tols = *tolerances ;
  switch (start)
  { case startCold:
    { lcl_opts = *initialSolveOptions ;
      break ; }
    case startWarm:
    { lcl_opts = *resolveOptions ;
      assert(lcl_opts.forcecold == false) ;
      break ; }
    case startHot:
    { lcl_opts = *resolveOptions ;
      assert(lcl_opts.forcecold == false) ;
      lcl_opts.forcewarm = false ;
      break ; }
    case startInvalid:
    { handler_->message(ODSI_CONFUSION,messages_)
        << __LINE__ << CoinMessageEol ;
      return (lpFATAL) ; } }
  dy_checkdefaults(consys,&lcl_opts,&lcl_tols) ;
  lpprob->phase = dyINV ;

# ifdef ODSI_STATISTICS
  if (statistics) dy_freestats(&statistics) ;
  dy_initstats(&statistics,lpprob->consys) ;
# endif

  retries = 0 ;
/*
  Skim off any persistent flags that we may need to reset after a failure.
  (Output requests, etc.)
*/
  persistent_flags = getflg(lpprob->ctlopts,lpctlACTVARSOUT) ;

/*
  Step through the constraints and replace ax >= b constraints with (-a)x <=
  -b constraints. consys_mulrow will take care of the necessary
  modifications. Normally this would be handled in mpsin, along with the
  deletion of empty constraints, but the OSI test suite isn't tolerant of
  changing the sense of constraints, let alone removing a few. Arguably a
  good thing.
*/
  flipped = (bool *) CALLOC(lpprob->consys->concnt+1,sizeof(bool)) ;
  flips = 0 ;
  for (ndx = lpprob->consys->concnt ; ndx > 0 ; ndx--)
  { if (lpprob->consys->ctyp[ndx] == contypGE)
    { if (consys_mulrow(lpprob->consys,ndx,-1) == false)
      { errmsg(112,rtnnme,lpprob->consys->nme,"scalar multiply","row",
	       consys_nme(lpprob->consys,'c',ndx,false,NULL),ndx) ;
	FREE(flipped) ;
	return (lpFATAL) ; }
      flipped[ndx] = true ;
      flips++ ; } }

/*
  Take an initial run at doing the lp as it comes in. If this doesn't work,
  we'll try harder. After the initial shot at the lp, we can clear the vector
  change flags (they're only relevant to a hot start).
*/
  lpret = dylp(lpprob,&lcl_opts,&lcl_tols,statistics) ;
# ifdef ODSI_INFOMSGS
  if (print >= 1)
  { if (lpret == lpOPTIMAL || lpret == lpINFEAS || lpret == lpUNBOUNDED)
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  success, status %s",
		  dy_prtlpret(lpprob->lpret)) ; }
    else
    if (lpret == lpITERLIM)
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n  premature termination, status %s",
		  dy_prtlpret(lpprob->lpret)) ; }
    else
    { dyio_outfmt(dy_logchn,dy_gtxecho,"\n  failed, status %s",
		  dy_prtlpret(lpprob->lpret)) ; } }
# endif
  clrflg(lpprob->ctlopts,lpctlUBNDCHG|lpctlLBNDCHG|lpctlRHSCHG|lpctlOBJCHG) ;

/*
  Did it work? If not, get down to trying to recover. The algorithm is
  to first try a cold start with the full system, then proceed to halve
  the refactor frequency.
  
  Since a cold start rebuilds all dylp data structures from scratch, we don't
  need to worry about the various *CHG flags. And if we're about to retry,
  dylp freed the data structure regardless of the NOFREE flag.
*/
  if (!(lpret == lpOPTIMAL || lpret == lpINFEAS || lpret == lpUNBOUNDED ||
	lpret == lpITERLIM))
  { 
#   ifdef DYLP_POSTMORTEM
/*
  Write an mps file with the failed problem.  This is the hard way to do
  this, given how close we are to dylp. But I wanted to try it to exercise
  the code. Arguably I should call destruct_*_cache in order that writeMps
  sees any constraint flips. But that's a consistency problem. Cached
  structural vectors (including row descriptions and matrices) are normally
  preserved across a call to a solve routine. If presolve is active, then
  they are saved in the presolve cache area, but if we're solving without
  presolve, they are still in the main cache area. Right now, seems more
  trouble to sort out than it's worth.
*/
    int saveInfo ;
    void *foo = &saveInfo ;
    bool saveSense,saveEcho[2] ;
    OsiHintStrength saveStrength ;
    getHintParam(OsiDoReducePrint,saveSense,saveStrength,foo) ;
    saveEcho[0] = initial_gtxecho ;
    saveEcho[1] = resolve_gtxecho ;

    int level = 4 ;
    level |= 0x10 ;
    setHintParam(OsiDoReducePrint,true,OsiForceDo,&level) ;
    lcl_opts.print = initialSolveOptions->print ;
    dy_gtxecho = initial_gtxecho ;
    std::cout << "Verbosity now maxed at " << level << ".\n" ;
    writeMps("dylpPostmortem","mps") ;
#   endif
    if (lcl_opts.forcecold == true) lcl_opts.factor /= 2 ;
    lcl_opts.forcecold = true ;   
    lcl_opts.fullsys = true ;
    lcl_tols.pfeas_scale *= 100 ;
    lcl_tols.dfeas_scale *= 100 ;
/*
  We're ready. Open a loop that'll call dylp, doing a cold start on each
  lp and halving the refactor frequency with each call. Don't forget to
  reset the persistent flags.
*/
    for ( ; lcl_opts.factor >= 10 ; lcl_opts.factor /= 2)
    { retries++ ;
#     ifdef ODSI_INFOMSGS
      if (print >= 1)
      { dyio_outfmt(dy_logchn,dy_gtxecho,".\n    retry %d: refactor = %d ...",
	            retries,lcl_opts.factor) ; }
#     endif
      setflg(lpprob->ctlopts,persistent_flags) ;
      lpprob->phase = dyINV ;
      lpret = dylp(lpprob,&lcl_opts,&lcl_tols,statistics) ;

      if (lpret == lpOPTIMAL || lpret == lpINFEAS || lpret == lpUNBOUNDED)
      {
#       ifdef ODSI_INFOMSGS
	if (print >= 1)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n  success, status %s",
		      dy_prtlpret(lpprob->lpret)) ; }
#       endif
	break ; }
#     ifdef ODSI_INFOMSGS
      else
      { if (print >= 1)
	{ dyio_outfmt(dy_logchn,dy_gtxecho,"\n  failed, status %s",
		      dy_prtlpret(lpprob->lpret)) ; } }
#     endif
    }
#   ifdef DYLP_POSTMORTEM
/*
  Note that sense and strength are unimportant, execpt that we need something
  stronger than OsiHintIgnore.
*/
    setHintParam(OsiDoReducePrint,saveSense,OsiHintTry,&saveInfo) ;
    initial_gtxecho = saveEcho[0] ;
    resolve_gtxecho = saveEcho[1] ;
#   endif
    }
/*
  Time to undo any constraint flips. We also have to tweak the corresponding
  duals --- flipping the sign of a row in the basis corresponds to flipping the
  sign of a column in the basis inverse, which means that the sign of the
  corresponding dual is flipped. (Or just look at it as yA = (-y)(-A).)  We
  need to walk the basis here. If dylp is in dynamic mode, there may be fewer
  active constraints than when we started.
*/
  if (flips > 0)
  { for (ndx = lpprob->consys->concnt ; ndx > 0 ; ndx--)
    { if (flipped[ndx] == true)
      { if (consys_mulrow(lpprob->consys,ndx,-1) == false)
	{ errmsg(112,rtnnme,lpprob->consys->nme,"scalar multiply","row",
		 consys_nme(lpprob->consys,'c',ndx,false,NULL),ndx) ;
	  FREE(flipped) ;
	  return (lpFATAL) ; } } }
    if (lpprob->y != NULL)
    { basis = lpprob->basis ;
      for (ndx = 1 ; ndx <= basis->len ; ndx++)
      { cndx = basis->el[ndx].cndx ;
	if (flipped[cndx] == true) lpprob->y[ndx] = -lpprob->y[ndx] ; } } }
  FREE(flipped) ;
  solnIsFresh = true ;
# if ODSI_TRACK_FRESH > 0
  std::cout
    << "ODSI(" << std::hex << this << std::dec
    << ")::solution refreshed." << std::endl ;
# endif
/*
  That's it, we've done our best. Do a little printing and return.
*/
# ifdef ODSI_INFOMSGS
  if (print >= 1)
  { if (lpprob->lpret == lpOPTIMAL)
      dyio_outfmt(dy_logchn,dy_gtxecho,"; objective %.8g",lpprob->obj) ;
    else
    if (lpprob->lpret == lpINFEAS)
      dyio_outfmt(dy_logchn,dy_gtxecho,"; infeasibility %.4g",lpprob->obj) ;
    if (lpprob->phase == dyDONE)
      dyio_outfmt(dy_logchn,dy_gtxecho," after %d pivots",lpprob->iters) ;
    dyio_outchr(dy_logchn,dy_gtxecho,'.') ;
    dyio_flushio(dy_logchn,dy_gtxecho) ; }
# endif

  return (lpret) ; }

//@} // LPCommonMethods




/*! \defgroup ColdStartMethods Cold Start Method
    \brief Method for solving an LP from a cold start.

    This method solves an LP from scratch. Any necessary initialisation of the
    dylp solver is performed, then dylp is called to solve the problem. Dylp
    uses the full constraint system for the initial solution, rather than
    dynamic simplex.
*/
//@{

/*!
  The first actions make sure all is in order: there must be a valid lpprob
  structure, the basis package must be initialised, and this ODSI object must
  own the solver. The routine then sets a few other dylp parameters to
  appropriate values for the initial solution of the lp and calls the
  solver.

  Forcing fullsys here is generally a good choice in terms of the performance
  of the code, but it does reduce flexibility. This should be made into a
  hint.
*/
void ODSI::initialSolve ()

{ CoinMessageHandler *hdl = messageHandler() ; 
  flags save_ctlopts = 0 ;
  bool save_finpurge_vars = false ;
  bool save_finpurge_cons = false ;
/*
  The minimum requirement is a constraint system. Options and tolerances should
  also be present --- they're established when the ODSI object is created.
*/
  assert(consys && initialSolveOptions && tolerances) ;
/*
  Do we have an lpprob structure yet?
*/
  if (!lpprob) construct_lpprob() ;
/*
  Is the basis package initialised? The second parameter controls how
  many basis updates the basis can hold before it requires refactoring.
  Adding 5 to dylp's refactor interval should give a safety margin.
*/
  if (basis_ready == false)
  { int count = static_cast<int>(1.5*getNumRows()) ;
    dy_initbasis(count,initialSolveOptions->factor+5,0) ;
    basis_ready = true ; }
/*
  Does some ODSI object (including this object) own the solver? If so, detach
  it. We're doing an initial solve and we'll want to start fresh.
*/
  if (dylp_owner != 0)
  { dylp_owner->detach_dylp() ; }
  clrflg(lpprob->ctlopts,lpctlDYVALID) ;
/*
  Are we going to do presolve and postsolve? If so, now's the time to presolve
  the constraint system.
*/
# ifdef ODSI_INFOMSGS
  double startTime = CoinCpuTime() ;
  double presolveTime = startTime ;
# endif
  bool presolving ;
  double bias = 0 ;
  { OsiHintStrength strength ;
    bool sense ;
    (void) getHintParam(OsiDoPresolveInInitial,sense,strength) ;
    if (sense == true)
    { presolving = true ; }
    else
    { presolving = false ; } }
  if (presolving == true)
  { preObj_ = initialisePresolve(false) ;
    doPresolve() ;
    if (evalPresolve() == true)
    { saveOriginalSys() ;
      installPresolve() ;
      bias = preObj_->dobias_ ; }
    else
    { destruct_presolve() ;
      presolving = false ; }
#   ifdef ODSI_INFOMSGS
    presolveTime = CoinCpuTime() ;
#   endif
  }
/*
  Remove any active basis and trash any cached solution. The calls to
  destruct_*_cache remove only the cached solution vectors. If we're
  presolving, saveOriginalSys has moved the cached structural vectors and
  matrices to a save place, otherwise they'll be unaffected.
*/
# if ODSI_TRACK_ACTIVE > 0
  std::cout
    << "ODSI(" << std::hex << this
    << ")::initialSolve(1): deleting active basis "
    << activeBasis << std::dec
    << "." << std::endl ;
# endif
  delete activeBasis ;
  activeBasis = 0 ;
  activeIsModified = false ;
  destruct_col_cache(false) ;
  destruct_row_cache(false) ;
/*
  Establish logging and echo values, and invoke the solver.
*/
  if (dyio_isactive(local_logchn)) dy_logchn = local_logchn ;
  dy_gtxecho = initial_gtxecho ;
  if (presolving == true)
  { save_ctlopts = getflg(lpprob->ctlopts,lpctlNOFREE|lpctlACTVARSOUT) ;
    clrflg(lpprob->ctlopts,lpctlNOFREE|lpctlACTVARSOUT) ;
    save_finpurge_vars = initialSolveOptions->finpurge.vars ;
    initialSolveOptions->finpurge.vars = false ;
    save_finpurge_cons = initialSolveOptions->finpurge.cons ;
    initialSolveOptions->finpurge.cons = false ; }
  lp_retval = do_lp(startCold) ;
  if (presolving == true)
  { lpprob->ctlopts = setflg(lpprob->ctlopts,save_ctlopts) ;
    initialSolveOptions->finpurge.vars = save_finpurge_vars ;
    initialSolveOptions->finpurge.cons = save_finpurge_cons ; }
# ifdef ODSI_INFOMSGS
  double firstLPTime = CoinCpuTime() ;

  hdl->message(ODSI_COLD,messages_) ;
  hdl->printing(presolving)
    << "presol" ;
  hdl->printing(true)
    << dy_prtlpret(lp_retval)
    << getObjSense()*(lpprob->obj+bias) << lpprob->iters
    << CoinMessageEol ;
# endif
/*
  Separate the failure cases from the successful cases. lpITERLIM is
  questionable in this context (initial solution to the lp) but it's
  considered ok in warm and hot start (in particular, for strong branching).
  Since we can't tell from here, include it in the successes.
*/
  bool lpOK ;
  if ((lp_retval == lpOPTIMAL || lp_retval == lpINFEAS ||
       lp_retval == lpUNBOUNDED || lp_retval == lpITERLIM))
  { lpOK = true ; }
  else
  { lpOK = false ; }

# ifdef ODSI_INFOMSGS
  dylp_printsoln(true,true) ;
# endif
/*
  If we did presolve, we now need to do a postsolve, install a warm start,
  and call dylp one more time to set up the optimal solution with the
  original constraint system. Note that initialisePostsolve destroys the
  presolve object (by converting it to the postsolve object).
  installPostsolve brings back the cached original constraint system, the
  cached structural vectors and matrices, and creates a warm start
  appropriate for the original system.

  When we're done, we again need to clean out the active basis that was
  established as we set up for the second call to dylp.
*/
  if (lpOK)
  {
#   ifdef ODSI_INFOMSGS
    double postsolveTime = firstLPTime ;
    double secondLPTime = firstLPTime ;
#   endif
    int presolIters = lpprob->iters ;
    if (presolving == true)
    { postObj_ = initialisePostsolve(preObj_) ;
      doPostsolve() ;
      installPostsolve() ;
#     ifdef ODSI_INFOMSGS
      postsolveTime = CoinCpuTime() ;
#     endif
      lp_retval = do_lp(startWarm) ;
#     ifdef ODSI_INFOMSGS
      secondLPTime = CoinCpuTime() ;
      hdl->message(ODSI_COLD,messages_) ;
      hdl->printing(true)
	<< "postsol" ;
      hdl->printing(true)
	<< dy_prtlpret(lp_retval) << getObjSense()*lpprob->obj << lpprob->iters
	<< CoinMessageEol ;
#     endif
      if (!(lp_retval == lpOPTIMAL || lp_retval == lpINFEAS ||
	    lp_retval == lpUNBOUNDED))
      { throw CoinError("Call to dylp failed (postsolve).",
			"initialSolve","OsiDylpSolverInterface") ; }
      lpprob->iters += presolIters ;
#     if ODSI_TRACK_ACTIVE > 0
      std::cout
	<< "ODSI(" << std::hex << this
	<< ")::initialSolve(2): deleting active basis "
	<< activeBasis << std::dec
	<< "." << std::endl ;
#     endif
      delete activeBasis ;
      activeBasis = 0 ;
      activeIsModified = false ; }
#   ifdef ODSI_INFOMSGS
    hdl->message(ODSI_SHORTSTATS,messages_)
      << consys->nme << secondLPTime-startTime
      << presolveTime-startTime
      << presolIters << firstLPTime-presolveTime
      << postsolveTime-firstLPTime
      << lpprob->iters-presolIters << secondLPTime-postsolveTime
      << CoinMessageEol ;
#   endif
  }
/*
  LP failure. Destroy the presolve object.
*/
  else
  { if (presolving == true)
    { destruct_presolve() ; } }
/*
  There should be no cached solution vectors at this point.
*/
  assert(_col_x == 0) ;
  assert(_col_cbar == 0) ;
  assert(_row_lhs == 0) ;
  assert(_row_price == 0) ;
/*
  Tidy up. If all went well, indicate this object owns the solver, set the
  objective, and set the active basis. If we've failed, do the opposite.

  Any of lpOPTIMAL, lpINFEAS, or lpUNBOUNDED can (in theory) be hot started
  after allowable problem modifications, hence can be flagged lpctlDYVALID.
  dylp overloads lpprob->obj with the index of the unbounded variable when
  returning lpUNBOUNDED, so we need to fake the objective.
*/
  if (lpOK && flgon(lpprob->ctlopts,lpctlDYVALID))
  { dylp_owner = this ;
#   ifdef ODSI_INFOMSGS
    hdl->message(ODSI_ATTACH,messages_)
      << "initialSolve" << (int) reinterpret_cast<CoinIntPtr>(this)
      << CoinMessageEol ;
#   endif
    if (lpprob->lpret == lpUNBOUNDED)
    { _objval = -getObjSense()*getInfinity() ; }
    else
    { _objval = getObjSense()*lpprob->obj ; }
    activeBasis = this->getWarmStart() ;
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::initialSolve: setting active basis "
      << activeBasis << std::dec
      << "." << std::endl ;
#   endif
  }
  else
  { dylp_owner = 0 ; }

  return ; }

//@} // ColdStartMethods



/*! \defgroup SolverTerm Methods Returning Solver Termination Status */
//@{


inline int ODSI::getIterationCount () const

{ return ((lpprob)?lpprob->iters:0) ; }


inline bool ODSI::isIterationLimitReached () const

{ return (lp_retval == lpITERLIM) ; }


inline bool ODSI::isProvenOptimal () const

{ return (lp_retval == lpOPTIMAL) ; }


inline bool ODSI::isProvenPrimalInfeasible () const

{ return (lp_retval == lpINFEAS) ; }


/*!
  Aka primal unbounded. There is some loss here, as dylp does not have
  a way to indicate both primal and dual infeasibility.
*/

inline bool ODSI::isProvenDualInfeasible () const

{ return (lp_retval == lpUNBOUNDED) ; }


/*!
  Returns true if dylp abandoned the problem. This could be due to numerical
  problems (accuracy check or singular basis), stalling, unexpected loss of
  feasibility, inability to allocate space, or other fatal internal error.
*/
/*
  As implemented here, we simply look for the standard codes that indicate
  optimality, unboundedness, or infeasibility. Anything else qualifies as
  abandoned.  Hitting the iteration limit is specifically not included here
  --- use isIterationLimitReached to check that. There are applications
  (e.g., strong branching) that deliberately use an artificially low
  iteration limit.
*/

inline bool ODSI::isAbandoned () const

{ if (lp_retval == lpOPTIMAL || lp_retval == lpUNBOUNDED ||
      lp_retval == lpINFEAS || lp_retval == lpITERLIM)
    return (false) ;
  else
    return (true) ; }


//@} // SolverTerm



/*! \defgroup SolnInfo Methods to Get Solution Information */
//@{

/*!
  Returns the objective function value.

  The default implementation calculates the dot product of the objective
  and the current primal solution, then subtracts any constant offset.
  Note that if there is no current solution, the objective will be zero in the
  absence of an offset.
*/

inline double ODSI::getObjValue () const
/*
  _objval is kept up-to-date and corrected for obj_sense. All we need to do
  here is deal with the offset.
*/
{ double objoffset ;

  getDblParam(OsiObjOffset,objoffset) ;

  return (_objval-objoffset) ; }
  


/*!
  Return a cached vector of primal variable values. If there is no cached copy,
  construct one and cache it.

  \note The values returned for superbasic variables are bogus (a limitation
	of the representation chosen to return the solution). Dylp ensures
	that an optimal or infeasible solution will not contain them. They
	may appear in an unbounded solution. The value returned is -inf; this
	ensures that you'll know if you inadvertently include it in a
	calculation.
*/

const double* ODSI::getColSolution () const
/*
  There are three possible sources of a primal solution: a call to dylp, a call
  to pessimal_primal(), or specified by the client (setColSolution).

  If we have a cached solution in _col_x, we'll use it. setColSolution writes
  the solution directly into _col_x, as does pessimal_primal().

  If the solver is empty (no consys), return null.
  
  If an lpprob exists and is fresh, then we have a solution from dylp.  Calls
  to the solver will invalidate cached data, including _col_x, but will not
  regenerate it, hence we may need to do that here.  Dylp returns a vector x
  of basic variables, in basis order, and a status vector. To assemble a
  complete primal solution, it's necessary to construct one vector that
  merges the two sources. The result is cached.

  NOTE the caution above about superbasic (vstatSB) variables. Nonbasic free
  variables (vstatNBFR) occur under the same termination conditions, but are
  legitimately zero. We really should return NaN for superbasics, using the
  function <limits>:std::numeric_limits::signaling_NaN(), but it's just not
  worth the hassle. Sun Workshop C++ (Rogue Wave Software) has it, but Gnu
  C++ prior to v3 doesn't even have <limits>, and the local v3 installation's
  <limits> points to another, non-existent header. I don't want to even
  discuss Microsquash.
*/
{ if (_col_x)
  { return _col_x ; }
  if (consys == 0)
  { return (0) ; }

  if (!solnIsFresh)
  { CoinMessageHandler *hdl = messageHandler() ; 
    hdl->message(ODSI_ACCESS_STALE,messages_)
      << "getColSolution"
      << CoinMessageEol ;
#   if ODSI_TRACK_SOLVERS > 0
    std::cout
      << "ODSI(" << std::hex << this << std::dec
      << ")::getColSolution: request for stale solution." << std::endl ;
#   endif
#   ifdef ODSI_STRICTLY_FRESH
    throw CoinError("Constraint system has changed since last call to solver.",
		    "getColSolution","OsiDylpSolverInterface") ;
    return (0) ;
#   endif
  }

  assert(lpprob->status && lpprob->x) ;

  int n = getNumCols() ;
  flags statj ;
  _col_x = new double[n] ;

/*
  Walk the status vector, taking values for basic variables from lpprob.x and
  values for nonbasic variables from the appropriate bounds vector.
*/
  for (int j = 0 ; j < n ; j++)
  { statj = lpprob->status[idx(j)] ;
    if (((int) statj) < 0)
    { int k = -((int) statj) ;
      _col_x[j] = lpprob->x[k] ; }
    else
    { switch (statj)
      { case vstatNBLB:
	case vstatNBFX:
	{ _col_x[j] = consys->vlb[idx(j)] ;
	  break ; }
	case vstatNBUB:
	{ _col_x[j] = consys->vub[idx(j)] ;
	  break ; }
	case vstatNBFR:
	{ _col_x[j] = 0 ;
	  break ; }
	case vstatSB:
	{ _col_x[j] = -odsiInfinity ;
	  break ; } } } }

  return (_col_x) ; }


/*!
  If necessary, a cached copy is built from the solution returned by the
  solver.
*/

const double* ODSI::getRowPrice () const
/*
  There are two possible sources of duals: a call to the solver, or the
  client (setRowPrice).

  If we have a cached solution in _row_price, we'll use it.  setRowPrice
  writes the duals directly into _row_price.

  If the solver is empty (no consys), return null.

  If an lpprob exists, we have a solution from dylp.  Dylp reports dual
  variables in basis order, so we need to write a vector with the duals
  dropped into position by row index.
*/

{ if (_row_price)
  { return (_row_price) ; }
  if (consys == 0)
  { return (0) ; }

  if (!solnIsFresh)
  { CoinMessageHandler *hdl = messageHandler() ; 
    hdl->message(ODSI_ACCESS_STALE,messages_)
      << "getRowPrice"
      << CoinMessageEol ;
#   ifdef ODSI_STRICTLY_FRESH
    throw CoinError("Constraint system has changed since last call to solver.",
		    "getRowPrice","OsiDylpSolverInterface") ;
    return (0) ;
#   endif
  }

  assert(lpprob->basis && lpprob->y) ;

  int m = getNumRows() ;
  _row_price = new double[m] ;
  basis_struct* basis = lpprob->basis ;

  memset(_row_price,0,m*sizeof(double)) ;

  for (int k = 1 ; k <= basis->len ; k++)
  { int i = inv(basis->el[k].cndx) ;
    _row_price[i] = lpprob->y[k]*getObjSense() ; }

  return (_row_price) ; }


/*!
  Return a cached vector of row activity. If there is no cached copy,
  calculate and cache Ax (i.e., the constraints evaluated at the primal
  solution).
*/

const double *ODSI::getRowActivity () const
/*
  This routine calculates the value of the left-hand-side of a constraint. The
  algorithm is straightforward: Retrieve the primal solution, then walk the
  variables, adding the contribution of each non-zero variable.
*/
{ 
/*
  If we have a cached copy, we're done.
*/
  if (_row_lhs) return (_row_lhs) ;
/*
  In order to go further, we'll need a primal solution and some rows.
*/
  int m = getNumRows() ;
  const double *x = getColSolution() ;

  if (m > 0 && x)
/*
  Create a vector to hold ax and clear it to 0.
*/
  { _row_lhs = new double[consys->concnt] ;
    memset(_row_lhs,0,m*sizeof(double)) ;
/*
  Walk the primal solution vector, adding in the contribution from non-zero
  variables. Take care that we don't scale infinity (sigh ... I hate this
  finite infinity business).
*/
    pkvec_struct *aj = pkvec_new(m) ;
    for (int j = 0 ; j < consys->varcnt ; j++)
    { double xj = x[j] ;
      if (xj != 0)
      { bool r = consys_getcol_pk(consys,idx(j),&aj) ;
	if (!r)
	{ delete[] _row_lhs ;
	  _row_lhs = 0 ;
	  if (aj) pkvec_free(aj) ;
	  return (0) ; }
	if (fabs(xj) >= odsiInfinity)
	{ for (int l = 0 ; l < aj->cnt ; l++)
	  { int i = inv(aj->coeffs[l].ndx) ;
	    if (fabs(_row_lhs[i]) < odsiInfinity)
	    { double aij = aj->coeffs[l].val ;
	      if (aij < 0)
	      { _row_lhs[i] = -xj ; }
	      else
	      if (aij > 0)
	      { _row_lhs[i] = xj ; } } } }
	else
	{ for (int l = 0 ; l < aj->cnt ; l++)
	  { int i = inv(aj->coeffs[l].ndx) ;
	    if (fabs(_row_lhs[i]) < odsiInfinity)
	    { _row_lhs[i] += xj*aj->coeffs[l].val ; } } } } }
    if (aj) pkvec_free(aj) ;
/*
  Groom the vector to eliminate tiny values.
*/
    for (int i = 0 ; i < consys->concnt ; i++)
      setcleanzero(_row_lhs[i],tolerances->zero) ; }
/*
  Cache the result and return.
*/

  return (_row_lhs) ; }


/*!
  The routine calculates cbar = (c - yA) by accumulating y[i]*a<i>.
*/

const double *ODSI::getReducedCost () const
/*
  Calculate the reduced cost as cbar = c - yA.

  It's tempting to want to dive into dylp for this calculation, but there are
  problems: First, some other ODSI object may own the solver. Second, we
  can't be sure where the duals are coming from --- they could be part of a
  solution from dylp, or they could be a vector supplied by the client.
  Third, both the objective coefficients and the duals are subject to sign
  change for maximisation. In the end, it seems best to rely on
  getObjCoefficients and getRowPrice to provide the objective and duals,
*/

{ 
/*
  Check for a cached value.
*/
  if (_col_cbar) return (_col_cbar) ;
/*
  Do we have columns to to work with? If yes, then we have a valid constraint
  system data structure. Initialize cbar with the objective coefficients.
*/
  int n = getNumCols() ;
  if (n == 0) return (0) ;
  _col_cbar = new double[n] ;
  COPY_VEC(double,getObjCoefficients(),_col_cbar,n) ;  
/*
  Do we have rows? Dual variables? If not, we're done. The presence of rows
  does not necessarily mean we have valid duals.
*/
  int m = getNumRows() ;
  const double *y = getRowPrice() ;
  if (!y) return (_col_cbar) ;
/*
  We have a constraint system, objective coefficients, and dual variables.
  Grind out the calculation, row by row.
*/
  pkvec_struct *ai = pkvec_new(n) ;
  for (int i = 0 ; i < m ; i++)
  { if (y[i] != 0)
    { bool r = consys_getrow_pk(consys,idx(i),&ai) ;
      if (!r)
      { delete[] _col_cbar ;
	_col_cbar = 0 ;
	if (ai) pkvec_free(ai) ;
	return (0) ; }
      for (int l = 0 ; l < ai->cnt ; l++)
      { int j = inv(ai->coeffs[l].ndx) ;
	_col_cbar[j] -= y[i]*ai->coeffs[l].val ; } } }
  if (ai) pkvec_free(ai) ;
/*
  Groom the vector to eliminate tiny values and we're done.
*/
  for (int j = 0 ; j < n ; j++) setcleanzero(_col_cbar[j],tolerances->cost) ;

  return (_col_cbar) ; }


//@} // SolnInfo



/*! \defgroup GetProbInfo Methods to Obtain Problem Information */
//@{

/*!
  Nonexistant variables return false. In the absence of type information, all
  variables are continuous.
*/

inline bool ODSI::isBinary (int i) const

{ if (!consys || i < 0 || i > consys->varcnt-1)
    return (false) ;
  else
  if (!consys->vtyp)
    return (false) ;
  else
    return (consys->vtyp[idx(i)] == vartypBIN) ; }


/*!
  Nonexistant variables return false. In the absence of type information, all
  variables are continuous.
*/

inline bool ODSI::isIntegerNonBinary (int i) const

{ if (!consys || i < 0 || i > consys->varcnt-1)
    return (false) ;
  else
  if (!consys->vtyp)
    return (false) ;
  else
    return (consys->vtyp[idx(i)] == vartypINT) ; }


/*!
  Nonexistant variables return false. In the absence of type information, all
  variables are continuous.
*/

inline bool ODSI::isInteger (int i) const

{ if (!consys || i < 0 || i > consys->varcnt-1)
    return (false) ;
  else
  if (!consys->vtyp)
    return (false) ;
  else
    return (consys->vtyp[idx(i)] == vartypINT ||
	    consys->vtyp[idx(i)] == vartypBIN) ; }



/*!
  Nonexistant variables return false. In the absence of type information, all
  variables are continuous.
*/

inline bool ODSI::isContinuous (int i) const

{ if (!consys || i < 0 || i > consys->varcnt-1)
    return (false) ;
  else
  if (!consys->vtyp)
    return (true) ;
  else
    return (consys->vtyp[idx(i)] == vartypCON) ; }


/*! Returns a pointer to the dylp data structure, with indexing shift. */

inline const double* ODSI::getColLower () const

{ if (!consys || !consys->vlb)
    return (0) ;
  else
    return (INV_VEC(double,consys->vlb)) ; }


/*! Returns a pointer to the dylp data structure, with indexing shift. */

inline const double* ODSI::getColUpper () const

{ if (!consys || !consys->vub)
    return (0) ;
  else
    return (INV_VEC(double,consys->vub)) ; }


inline int ODSI::getNumCols () const

{ if (!consys)
    return (0) ;
  else
    return (consys->varcnt) ; }


inline int ODSI::getNumIntegers () const

{ if (!consys)
    return (0) ;
  else
    return (consys->intvcnt+consys->binvcnt) ; }


inline int ODSI::getNumElements () const

{ if (!consys)
    return (0) ; 
  else
    return (consys->mtx.coeffcnt) ; }


inline int ODSI::getNumRows () const

{ if (!consys)
    return (0) ;
  else
    return (consys->concnt) ; }


/*! 
  Creates a cached copy, if it doesn't already exist, compensating for the
  objective sense. Returns a pointer to the cached copy.
*/

inline const double* ODSI::getObjCoefficients () const

{ if (!consys || !consys->obj) return (0) ;

  if (_col_obj) return _col_obj ;

  int n = getNumCols() ;
  _col_obj = new double[n] ;
  double *obj_0base = INV_VEC(double,consys->obj) ;

  if (getObjSense() < 0)
  { std::transform(obj_0base,obj_0base+n,_col_obj,std::negate<double>()) ; }
  else
  { COPY_VEC(double,obj_0base,_col_obj,n) ; }

  return (_col_obj) ; }


inline double ODSI::getObjSense () const

{ return (obj_sense) ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy. The row-major OSI matrix is generated from a column-major
  OSI matrix which is created if it doesn't already exist.
*/

const CoinPackedMatrix* ODSI::getMatrixByRow () const

{ if (!consys) return 0 ;
  if (_matrix_by_row) return _matrix_by_row ;

  _matrix_by_row = new CoinPackedMatrix ;
  _matrix_by_row->reverseOrderedCopyOf(*getMatrixByCol()) ;

  return _matrix_by_row ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  The routine generates the core structures of the OSI matrix: the paired
  arrays values (coefficients) and indices (row indices for coefficients),
  and the column description arrays start (column starting index in
  values/indices) and length (number of coefficients). These are filled from
  the consys_struct and then passed to assignMatrix to create the OSI matrix.

  This routine knows the internal structure of a consys_struct. This saves
  the overhead of a call to consys_getcol_pk and a packed vector intermediary.
  A debatable decision.
*/

const CoinPackedMatrix* ODSI::getMatrixByCol () const

{ if (!consys) return 0 ;
  if (_matrix_by_col) return _matrix_by_col ;

/*
  Get the column and coefficient counts and create the OSI core vectors.
*/
  int col_count = getNumCols() ;
  int coeff_count = consys->mtx.coeffcnt ;
  assert((col_count > 0 && coeff_count > 0) ||
	 (col_count == 0 && coeff_count == 0)) ;

  int* start = new int[col_count+1] ;
  int* len = new int[col_count] ;
  double* values = new double[coeff_count] ;
  int* indices = new int[coeff_count] ;
  CoinPackedMatrix* matrix = new CoinPackedMatrix ;
/*
  Scan out the coefficients from the consys_struct column by column.
*/
  colhdr_struct** cols = consys->mtx.cols ;
  int coeff_ndx = 0 ;
  for (int i = 0 ; i < col_count ; i++)
  { start[i] = coeff_ndx ;
    len[i] = cols[idx(i)]->len ;

    coeff_struct* c = cols[idx(i)]->coeffs ;
    for (int j = 0 ; j < len[i] ; j++, c = c->colnxt, coeff_ndx++)
    { values[coeff_ndx] = c->val ;
      indices[coeff_ndx] = inv(c->rowhdr->ndx) ; }
    assert(c == 0) ; }
  assert(coeff_ndx == coeff_count) ;
/*
  Create the OSI sparse matrix structure.
*/
  start[col_count] = coeff_ndx ;
  int row_count = getNumRows() ;
  matrix->assignMatrix(true, row_count, col_count, coeff_count, 
    values, indices, start, len) ;
/*
  Make a cached copy and return a pointer to it.
*/
  _matrix_by_col = matrix ;
  return matrix ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  A little work is required here to synthesize a vector that meets OSI's
  expectations. Dylp uses only rhs for all except range constraints and
  does not keep an explicit second bound. For non-binding (NB) constraints,
  it keeps no rhs values at all.
*/

const double* ODSI::getRowLower () const

{ if (!consys) return (0) ;
  if (_row_lower) return _row_lower ;

  int m = getNumRows() ;
  double* lower = new double[m] ;

  for (int i = 0 ; i < m ; i++)
  { contyp_enum ctypi = consys->ctyp[idx(i)] ;
    switch (ctypi)
    { case contypEQ:
      case contypGE:
      { lower[i] = consys->rhs[idx(i)] ;
	break ; }
      case contypRNG:
      { lower[i] = consys->rhslow[idx(i)] ;
	break ; }
      case contypLE:
      case contypNB:
      { lower[i] = -odsiInfinity ;
	break ; }
      default:
      { assert(0) ; } } }

  _row_lower = lower ;
  return (lower) ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  A little work is required here to synthesize a vector that meets OSI's
  expectations. Dylp uses only rhs for all except range constraints and
  does not keep an explicit second bound. For non-binding (NB) constraints,
  it keeps no rhs values at all.
*/

const double* ODSI::getRowUpper () const

{ if (!consys) return (0) ;
  if (_row_upper) return _row_upper ;

  int m = getNumRows() ;
  double* upper = new double[m] ;

  for (int i = 0 ; i < m ; i++)
  { if (consys->ctyp[idx(i)] == contypGE || consys->ctyp[idx(i)] == contypNB)
      upper[i] = odsiInfinity ;
    else
      upper[i] = consys->rhs[idx(i)] ; }

  _row_upper = upper ;
  return upper ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  This is a simple translation from dylp contyp_enum codes to the character
  codes used by OSI.
*/

const char* ODSI::getRowSense () const

{ if (!consys) return (0) ;
  if (_row_sense) return _row_sense ;
  
  int m = getNumRows() ;
  char* sense = new char[m] ;
  const contyp_enum* ctyp = INV_VEC(contyp_enum,consys->ctyp) ;

  std::transform(ctyp,ctyp+m,sense,type_to_sense) ;

  _row_sense = sense ;
  return sense ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  Dylp keeps explicit upper and lower row bounds only for range constraints.
  We need to construct a vector with information for all constraints. Note
  that by using getRowLower, getRowUpper, and getRowSense, we'll create cached
  copies of that information also.
*/

const double* ODSI::getRowRange () const

{ if (!consys) return (0) ;
  if (_row_range) return _row_range ;
  
  int m = getNumRows() ;  
  double* range = new double[m] ;
  const double* lower = getRowLower() ;
  const double* upper = getRowUpper() ;
  const char* sense = getRowSense() ;

  for (int i=0 ; i<m ; i++)
    if (sense[i] != 'R')
      range[i] = 0.0 ;
    else
      range[i] = upper[i] - lower[i] ;
  
  _row_range = range ;
  return range ; }


/*!
  Creates a cached copy if it doesn't already exist and returns a pointer to
  the cached copy.

  Some fiddling is required, as dylp's idea of how to store the
  right-hand-side information differs from OSI's.  The bulk of the translation
  is handled in getRowLower and getRowUpper. Here we just copy the correct
  value into rhs.
*/

const double* ODSI::getRightHandSide () const

{ if (!consys) return (0) ;
  if (_row_rhs) return _row_rhs ;
  
  int m = getNumRows() ;  
  double* rhs = new double[m] ;
  const double* lower = getRowLower() ;
  const double* upper = getRowUpper() ;
  const char* sense = getRowSense() ;

  for (int i = 0 ; i < m ; i++)
  { switch (sense[i])
    { case 'N':
      { rhs[i] = 0.0 ;
	break ; }
      case 'G':
      { rhs[i] = lower[i] ;
	break ; }
      case 'E':
      { rhs[i] = upper[i] ;
	break ; }
      case 'L':
      { rhs[i] = upper[i] ;
	break ; }
      case 'R':
      { rhs[i] = upper[i] ;
	break ; }
      default:
      { assert(0) ; } } }

  _row_rhs = rhs ;
  return (rhs) ; }


//@} // GetProbInfo



/*! \defgroup Unsupported Unsupported Methods */
//@{


/*!
  These functions are supported in dylp only when it's operating in the full
  constraint system mode.  Otherwise, the alternating addition and deletion of
  constraints and variables by the dynamic simplex algorithm means that the
  objective does not change monotonically.

  Right now, they're just a facade, until I train dylp to recognize the
  limits.
*/

bool ODSI::isDualObjectiveLimitReached () const

{ double objlim ;
  double objval = getObjValue() ;

  getDblParam(OsiDualObjectiveLimit,objlim) ;
  objlim *= getObjSense() ;

  if (getObjSense() > 0)
    return (objval > objlim) ;
  else
    return (objval < objlim) ; }


/*!
  Partially supported. See
  \link OsiDylpSolverInterface::isDualObjectiveLimitReached
	ODSI:isDualObjectiveLimitReached\endlink.
*/

bool ODSI::isPrimalObjectiveLimitReached () const

{ double objlim ;
  double objval = getObjValue() ;

  getDblParam(OsiDualObjectiveLimit,objlim) ;
  objlim *= getObjSense() ;

  if (getObjSense() > 0)
    return (objval < objlim) ;
  else
    return (objval > objlim) ; }


/*!
  \todo Will require modification to dylp --- currently it does not report
	this information.
*/

vector<double*> ODSI::getDualRays (int) const

{ throw CoinError("Unimplemented method.",
		  "getDualRays","OsiDylpSolverInterface") ;

  return (vector<double*>(0)) ; }


/*!
  \todo Will require modification to dylp --- currently it does not report
	this information.
*/

vector<double*> ODSI::getPrimalRays (int) const

{ throw CoinError("Unimplemented method.",
		  "getPrimalRays","OsiDylpSolverInterface") ;

  return (vector<double*>(0)) ; }


/*!
  dylp is strictly a simplex LP solver.
*/

void ODSI::branchAndBound ()

{ throw CoinError("Unimplemented method.",
		  "branchAndBound","OsiDylpSolverInterface") ;

  return ; }


//@}



/*! \defgroup WarmStart Warm Start Methods
    \brief Methods for solving an LP from a warm start.

  Dylp returns basis and status information as a matter of course when
  the solver termination status is optimal, infeasible, or unbounded. This
  information can be read back into dylp for a warm start.

  The \link OsiDylpSolverInterface::getWarmStart ODSI::getWarmStart \endlink
  routine simply copies the information into an OsiDylpWarmStartBasis object.
  This can be done after any call to a solver routine. getWarmStart will return
  an empty warm start object if no solution exists; this is occasionally useful
  in contexts where a basis will be generated from scratch or where it's
  necessary to establish the type of warm start object for later use.

  The \link OsiDylpSolverInterface::setWarmStart ODSI::setWarmStart \endlink
  routine rebuilds the required basis and status information from the
  OsiDylpWarmStartBasis object and sets ODSI options so that the solver will
  attempt a warm start at the next call to
  \link OsiDylpSolverInterface::resolve ODSI::resolve \endlink.

  The convention in OSI is that warm start objects derive from the base class
  CoinWarmStart. If you're using only getWarmStart and setWarmStart, that's
  really all you need to know. For additional details, see the documentation
  for the OsiDylpWarmStartBasis class.
*/
//@{

/*!
  This routine returns an empty OsiDylpWarmStartBasis object. Its purpose is
  to provide a way to give a client a warm start basis object of the
  appropriate type, which can resized and modified as desired.
*/

CoinWarmStart *ODSI::getEmptyWarmStart () const
{ return (dynamic_cast<CoinWarmStart *>(new OsiDylpWarmStartBasis())) ; }

/*!
  This routine constructs an OsiDylpWarmStartBasis structure from the basis
  and status information returned by the dylp solver. If no valid solution
  exists, an empty warm start object is returned.
*/

CoinWarmStart* ODSI::getWarmStart () const

/*
  This routine constructs an OsiDylpWarmStartBasis structure from the basis
  returned by dylp.

  Dylp takes the atttitude that in so far as it is possible, logicals should
  be invisible outside the solver. The status of nonbasic logicals is not
  reported, nor does dylp expect to receive it. For completeness, however,
  getWarmStart will synthesize status information for nonbasic logicals.
*/

{ int i,j,k,m,n ;
  flags statj ;

/*
  If we have an active basis, return a clone.
*/
  if (activeBasis) { return (activeBasis->clone()) ; }
/*
  Create an empty ODWSB object. If no solution exists, we can return
  immediately.
*/
  OsiDylpWarmStartBasis *wsb = new OsiDylpWarmStartBasis ;
  assert(wsb) ;
  if (!lpprob) return (wsb) ;
/*
  Size the ODWSB object, then grab pointers to the status vectors and the
  dylp basis vector.  setSize initialises structural and logical status
  vectors to isFree, but constraint status entries are initialised to
  indicate active (atLowerBound), for the benefit of the world at large.
  We'll need to fix that below.
*/
  m = consys->concnt ;
  n = consys->varcnt ;
  wsb->setSize(n,m) ;

  char *const strucStatus = wsb->getStructuralStatus() ;
  char *const artifStatus = wsb->getArtificialStatus() ;
  char *const conStatus = wsb->getConstraintStatus() ;
  basis_struct *basis = lpprob->basis ;

  if (lpprob->lpret == lpOPTIMAL)
    wsb->setPhase(dyPRIMAL2) ;
  else
    wsb->setPhase(dyPRIMAL1) ;
/*
  Clear the constraint status vector to isFree, then walk the basis and mark
  the active constraints and basic variables. Basic logical variables are
  entered as the negative of the constraint index.
*/
  for (i = 0 ; i < m ; i++)
  { setStatus(conStatus,i,CoinWarmStartBasis::isFree) ; }
  for (k = 1 ; k <= basis->len ; k++)
  { i = inv(basis->el[k].cndx) ; 
    setStatus(conStatus,i,CoinWarmStartBasis::atLowerBound) ;
    j = basis->el[k].vndx ;
    if (j < 0)
    { j = inv(-j) ;
      setStatus(artifStatus,j,CoinWarmStartBasis::basic) ; }
    else
    { j = inv(j) ;
      setStatus(strucStatus,j,CoinWarmStartBasis::basic) ; } }
/*
  Next, scan the conStatus array. For each inactive (loose => isFree)
  constraint, indicate that the logical is basic. For active (tight =>
  atLowerBound) constraints, if the corresponding logical isn't already
  marked as basic, use the value of the dual to select one of atLowerBound
  (negative dual) or (for range constraints) atUpperBound (positive dual).
*/
  const double *y = getRowPrice() ;
  for (i = 0 ; i < m ; i++)
  { if (getStatus(conStatus,i) == CoinWarmStartBasis::isFree)
    { setStatus(artifStatus,i,CoinWarmStartBasis::basic) ; }
    else
    if (getStatus(artifStatus,i) == CoinWarmStartBasis::isFree)
    { if (y[i] > 0)
      { setStatus(artifStatus,i,CoinWarmStartBasis::atUpperBound) ; }
      else
      { setStatus(artifStatus,i,CoinWarmStartBasis::atLowerBound) ; } } }
/*
  Now scan the status vector and record the status of nonbasic structural
  variables. Some information is lost here --- CoinWarmStartBasis::Status
  doesn't encode NBFX.
*/
  for (j = 1 ; j <= n ; j++)
  { statj = lpprob->status[j] ;
    if (((int) statj) > 0)
    { switch (statj)
      { case vstatNBLB:
	case vstatNBFX:
	{ setStatus(strucStatus,inv(j), CoinWarmStartBasis::atLowerBound) ;
	  break ; }
	case vstatNBUB:
	{ setStatus(strucStatus,inv(j), CoinWarmStartBasis::atUpperBound) ;
	  break ; }
	case vstatNBFR:
	{ setStatus(strucStatus,inv(j), CoinWarmStartBasis::isFree) ;
	  break ; }
	default:
	{ delete wsb ;
	  wsb = 0 ;
	  throw CoinError("Warm start construction failed --- "
			  "invalid status in dylp basis.",
			  "getWarmStart","OsiDylpSolverInterface") ; } } } }

# if ODSI_PARANOIA >= 2
  if (wsb)
  { wsb->checkBasis() ; }
# endif

  return (wsb) ; }


/*!
  This routine installs the basis snapshot from an OsiDylpWarmStartBasis
  object and sets ODSI options so that dylp will attempt a warm start on the
  next call to \link OsiDylpSolverInterface::resolve ODSI::resolve \endlink.
  A basis with 0 rows and 0 columns, or a null parameter, is taken as a
  request to delete the existing warm start information held in activeBasis.

  Note that the size (rows x columns) of the CoinWarmStart information should
  match the size of the constraint system. The final basis built for dylp can
  be smaller, as inactive constraints will not be included. A basis with 0
  active constraints is legal. (In fact, fairly common. When a B&C code fixes
  all integer variables to confirm or regenerate a solution, a problem with
  only integer variables may have all variables fixed and no tight
  architectural constraints.)

  It can happen that the client will, for one reason or another, fix a basic
  variable. The symptom here is that we end up short a few basic variables.
  Unfortunately, the COIN data structures do not contain enough information
  for setWarmStart to tell whether this is intentional or an error. In the
  case of cbc, this does occur intentionally. So we compensate, and promote
  some random nonbasic variable to basic status. If the log level is
  sufficiently high, you'll get a message.
*/

bool ODSI::setWarmStart (const CoinWarmStart *ws)

{ int i,j,k ;
  CoinWarmStartBasis::Status osi_statlogi,osi_statconi,osi_statvarj ;

/*
  A null parameter says delete the current active basis and return.
*/
  if (!ws)
  { 
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::setWarmStart: deleting active basis "
      << activeBasis << std::dec
      << "(null param)." << std::endl ;
#   endif
    delete activeBasis ;
    activeBasis = 0 ;
    activeIsModified = false ;
    return (true) ; }
/*
  Use a dynamic cast to make sure we have a CoinWarmStartBasis. Then check
  the size --- 0 x 0 is just a request to remove the active basis.
*/
  const CoinWarmStartBasis *cwsb =
			dynamic_cast<const CoinWarmStartBasis *>(ws) ;
  if (!cwsb)
  { handler_->message(ODSI_NOTODWSB,messages_) << "Coin" ;
    return (false) ; }
  int varcnt = cwsb->getNumStructural() ;
  int concnt = cwsb->getNumArtificial() ;
  if (varcnt == 0 && concnt == 0)
  { 
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::setWarmStart: deleting active basis "
      << activeBasis << std::dec
      << " (0x0 CWSB)." << std::endl ;
#   endif
    delete activeBasis ;
    activeBasis = 0 ;
    activeIsModified = false ;
    return (true) ; }
/*
  Use a dynamic cast to see if we have an OsiDylpWarmStartBasis. If not,
  make one. 
*/
  const OsiDylpWarmStartBasis *wsb ;
  bool ourBasis = false ;
  wsb = dynamic_cast<const OsiDylpWarmStartBasis *>(ws) ;
  if (!wsb)
  { wsb = new OsiDylpWarmStartBasis(*cwsb) ;
    ourBasis = true ; }
# if ODSI_PARANOIA >= 2
  wsb->checkBasis() ;
# endif
  varcnt = wsb->getNumStructural() ;
  concnt = wsb->getNumArtificial() ;
  if (varcnt == 0 && concnt == 0)
  { 
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::setWarmStart: deleting active basis "
      << activeBasis << std::dec
      << " (0x0 ODWSB)." << std::endl ;
#   endif
    delete activeBasis ;
    activeBasis = 0 ;
    activeIsModified = false ;
    if (ourBasis == true) delete wsb ;
    return (true) ; }
/*
  Do we have an lpprob structure yet? If not, construct one.
*/
  assert(consys && consys->vlb && consys->vub) ;
  if (!lpprob) construct_lpprob() ;
  assert(resolveOptions) ;

/*
  Extract the info in the warm start object --- size and status vectors.  The
  number of variables and constraints in the warm start object should match
  the full size of the constraint system. Note that getWarmStart can create
  an empty ODWSB object (see comments with getWarmStart).

  Create a dylp basis_struct and status vector of sufficient size to hold the
  information in the OsiDylpWarmStartBasis. This space may well be freed or
  realloc'd by dylp, so use standard calloc to acquire it.  We'll only use as
  much of the basis as is needed for the active constraints.
*/
  if (!(varcnt == getNumCols() && concnt == getNumRows()))
  { handler_->message(ODSI_ODWSBBADSIZE,messages_)
      << concnt << varcnt << getNumRows() << getNumCols() ;
    if (ourBasis == true) delete wsb ;
    return (false) ; }

  const char *const strucStatus = wsb->getStructuralStatus() ;
  const char *const artifStatus = wsb->getArtificialStatus() ;
  const char *const conStatus = wsb->getConstraintStatus() ;

  flags *status = new flags[idx(varcnt)] ;
  basis_struct basis ;
  basis.el = new basisel_struct[idx(concnt)] ;
/*
  We never use these entries, but we need to initialise them to squash a
  `read from uninitialised memory' error during block copies.
*/
  basis.el[0].cndx = 0 ;
  basis.el[0].vndx = 0 ;
  status[0] = 0 ;
/*
  Walk the constraint status vector and build the set of active constraints.
*/
  int actcons = 0 ;
  for (i = 1 ; i <= concnt ; i++)
  { osi_statconi = getStatus(conStatus,inv(i)) ;
    if (osi_statconi == CoinWarmStartBasis::atLowerBound)
    { actcons++ ;
      basis.el[actcons].cndx = i ;
      basis.el[actcons].vndx = 0 ; } }
  basis.len = actcons ;
/*
  Now walk the structural status vector. For each basic variable, drop it into
  a basis entry and note the position in the status vector. For each nonbasic
  variable, set the proper status flag. We need to check bounds to see if the
  variable should be fixed.
*/
  k = 0 ;
  for (j = 1 ; j <= varcnt ; j++)
  { osi_statvarj = getStatus(strucStatus,inv(j)) ;
    switch (osi_statvarj)
    { case CoinWarmStartBasis::basic:
      { k++ ;
	assert(k <= actcons) ;
	basis.el[k].vndx = j ;
	status[j] = (flags) (-k) ;
	break ; }
      case CoinWarmStartBasis::atLowerBound:
      { if (consys->vlb[j] > -odsiInfinity)
	{ if (consys->vlb[j] == consys->vub[j])
	  { status[j] = vstatNBFX ; }
	  else
	  { status[j] = vstatNBLB ; } }
	else
	{ if (consys->vub[j] < odsiInfinity)
	  { status[j] = vstatNBUB ; }
	  else
	  { status[j] = vstatNBFR ; }
	  handler_->message(ODSI_ODWSBBADSTATUS,messages_)
	    << consys_nme(consys,'v',j,false,0) << j
	    << dy_prtvstat(vstatNBLB) << dy_prtvstat(status[j])
	    << CoinMessageEol ; }
	break ; }
      case CoinWarmStartBasis::atUpperBound:
      { if (consys->vub[j] < odsiInfinity)
	{ if (consys->vlb[j] == consys->vub[j])
	  { status[j] = vstatNBFX ; }
	  else
	  { status[j] = vstatNBUB ; } }
	else
	{ if (consys->vlb[j] > -odsiInfinity)
	  { status[j] = vstatNBLB ; }
	  else
	  { status[j] = vstatNBFR ; }
	  handler_->message(ODSI_ODWSBBADSTATUS,messages_)
	    << consys_nme(consys,'v',j,false,0) << j
	    << dy_prtvstat(vstatNBUB) << dy_prtvstat(status[j])
	    << CoinMessageEol ; }
	break ; }
      case CoinWarmStartBasis::isFree:
      { status[j] = vstatNBFR ;
	break ; } } }
/*
  Now we need to finish out the basis by adding the basic logicals.  The
  convention is to represent this with the negative of the index of the
  constraint which spawns the logical. Note that we do this only if the
  constraint is active. When running in fullsys mode, this is the common
  case. If we're running dynamic simplex, then all we're picking up here
  is the occasional logical that's basic at bound.
*/
  for (i = 1 ; i <= concnt ; i++)
  { osi_statlogi = getStatus(artifStatus,inv(i)) ;
    osi_statconi = getStatus(conStatus,inv(i)) ;
    if (osi_statlogi == CoinWarmStartBasis::basic &&
	osi_statconi == CoinWarmStartBasis::atLowerBound)
    { k++ ;
      assert(k <= actcons) ;
      basis.el[k].vndx = -i ; } }
/*
  Up to now we've worried about too many basic variables. Time to check to see
  if we have enough. If we're short, it could be an error, or it could be
  intentional. We can't tell from here (see note at head of routine). Fill
  out the basis by grabbing nonbasic logicals and declaring them basic. Such
  variables must exist. It'd be nice to correct the warm start, but we've
  been handed a const object, so we can't do that. If I was being picky here,
  I'd look for inequalities (slacks) in preference to equalities (artificials).
  But that would require keeping an artificial in reserve in case we didn't
  find a slack. My initial reaction is ``More trouble than it's worth.''
*/
  if (k < actcons)
  { handler_->message(ODSI_ODWSBSHORTBASIS,messages_)
      << consys->nme << k << actcons << CoinMessageEol ;
    for (i = 1 ; i <= concnt ; i++)
    { osi_statlogi = getStatus(artifStatus,inv(i)) ;
      if (osi_statlogi != CoinWarmStartBasis::basic)
      { k++ ;
	basis.el[k].vndx = -i ;
	if (k >= actcons) break ; } } }
  assert(k == actcons) ;
/*
  Now install the new basis in the lpprob_struct. If the present allocated
  capacity is sufficient, we can just copy over the existing information. If
  the capacity is insufficient, we'll free the existing space and allocate
  new. There's also the possibility that this WarmStartBasis came from some
  other object, and there are no vectors here at all. Because dylp relies on
  lpprob->colsze and lpprob->rowsze for the allocated capacity of actvars, x,
  and y, we need to reallocate them too.
*/
  if (lpprob->colsze < varcnt)
  { if (lpprob->status)
    { FREE(lpprob->status) ;
      lpprob->status = 0 ; }
    if (lpprob->actvars)
    { lpprob->actvars =
	(bool *) REALLOC(lpprob->actvars,idx(varcnt)*sizeof(bool)) ; }
    lpprob->colsze = varcnt ; }
  if (!lpprob->status)
  { lpprob->status = (flags *) CALLOC(idx(lpprob->colsze),sizeof(flags)) ; }
  if (!lpprob->actvars)
  { lpprob->actvars = (bool *) CALLOC(idx(lpprob->colsze),sizeof(bool)) ; }

  if (lpprob->rowsze < actcons)
  { if (lpprob->x)
    { lpprob->x =
	(double *) REALLOC(lpprob->x,idx(actcons)*sizeof(double)) ; }
    if (lpprob->y)
    { lpprob->y =
	(double *) REALLOC(lpprob->y,idx(actcons)*sizeof(double)) ; }
    if (lpprob->basis && lpprob->basis->el)
    { FREE(lpprob->basis->el) ;
      lpprob->basis->el = 0 ; }
    lpprob->rowsze = actcons ; }
  if (!lpprob->x)
  { lpprob->x = (double *) CALLOC(idx(lpprob->rowsze),sizeof(double)) ; }
  if (!lpprob->y)
  { lpprob->y = (double *) CALLOC(idx(lpprob->rowsze),sizeof(double)) ; }
  if (!lpprob->basis)
  { lpprob->basis =  (basis_struct *) CALLOC(1,sizeof(basis_struct)) ; }
  if (!lpprob->basis->el)
  { lpprob->basis->el =
      (basisel_struct *) CALLOC(idx(lpprob->rowsze),sizeof(basisel_struct)) ; }
/*
  Whew. The actual copy is dead easy.
*/
  copy_basis(&basis,lpprob->basis) ;
  COPY_VEC(flags,status,lpprob->status,idx(varcnt)) ;
  delete[] basis.el ;
  delete[] status ;
/*
  And we're done. The new status and basis are installed in lpprob, and the
  x, y, and actvars arrays have been resized if needed. Set the phase, signal
  a warm start, and then we're done.
*/
  lpprob->phase = wsb->getPhase() ;
  resolveOptions->forcecold = false ;
  resolveOptions->forcewarm = true ;
/*
  One last thing --- make a copy of wsb to be the active basis (that is, unless
  wsb is the active basis, and we're simply installing it in lpprob). If we
  already have a copy, so much the better.
*/
  if (wsb != activeBasis)
  { 
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::setWarmStart: replacing active basis "
      << activeBasis << " with " ;
#   endif
    delete activeBasis ;
    if (ourBasis == false)
    { activeBasis = wsb->clone() ; }
    else
    { activeBasis = const_cast<OsiDylpWarmStartBasis *>(wsb) ;
      ourBasis = false ; }
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << activeBasis << std::dec << "." << std::endl ;
#   endif
    activeIsModified = false ; }
  
  if (ourBasis == true) delete wsb ;

  return (true) ; }


void ODSI::resolve ()
/*
  Grossly oversimplified, this routine calls dylp after making sure that the
  forcecold option is turned off.  But we do need to make sure the solver is
  ours to use, and there's some state to maintain.

  If we're reoptimising, then the basis should be ready and we should have
  warm start information. For our purposes here, that boils down to the
  presence of an activeBasis object (created by setWarmStart() or a recent call
  to the solver).

  Note that we don't actually force a warm start unless we install the active
  basis here. Similarly, setWarmStart will force a warm start. But if we have
  an unmodified active basis, dylp should be able to manage a hot start.
*/
{ CoinMessageHandler *hdl = messageHandler() ; 
  assert(lpprob && lpprob->basis && lpprob->status &&
	 consys && resolveOptions && tolerances) ;

/*
  Does some other ODSI object own the solver? If so, detach it. If we don't
  own the solver, we surely don't have retained data structures.
*/
  if (dylp_owner != this)
  { if (dylp_owner != 0) dylp_owner->detach_dylp() ;
    clrflg(lpprob->ctlopts,lpctlDYVALID) ; }
/*
  We're reoptimising, so you might ask ``Why should we need to initialise the
  basis?''. Consider this scenario: we've optimised with ODSI object1, and
  captured a warm start. Then we destroy object1, create object2, install the
  warm start, and call resolve(). That's why we need to do this.  The second
  parameter controls how many basis updates the basis can hold before it
  requires refactoring.  Adding 5 to dylp's refactor interval should give a
  safety margin.
*/
  if (basis_ready == false)
  { int count = static_cast<int>((1.5*getNumRows()) + 2*getNumCols()) ;
    dy_initbasis(count,initialSolveOptions->factor+5,0) ;
    basis_ready = true ; }
/*
  We can hope that activeBasis is already installed in the lpprob, but there
  are several circumstances to trap for here:
    * The user has been playing with cuts since the last call to dylp, but
      kept within the rules for modifying the active basis. In this case
      activeIsModified will be true, and we need to install it.
    * The solver has been used by some other ODSI object since the last
      call by this object. In this case, we've just done a detach_dylp()
      and dylp_owner will be null. Again, we need to install activeBasis.
  If we have an active basis, install it and force dylp to warm start. If the
  installation fails, remove activeBasis and throw.
*/
  if (!activeBasis)
  { throw CoinError("Warm start failed --- empty active basis.",
		    "resolve","OsiDylpSolverInterface") ; }
  else
  if (activeBasis && (activeIsModified == true || dylp_owner == 0))
  { if (setWarmStart(activeBasis) == false)
    { 
#     if ODSI_TRACK_ACTIVE > 0
      std::cout
	<< "ODSI(" << std::hex << this
	<< ")::resolve: deleting basis "
	<< activeBasis << std::dec
	<< "." << std::endl ;
#     endif
      delete activeBasis ;
      activeBasis = 0 ;
      activeIsModified = false ;
      
      throw CoinError("Warm start failed --- invalid active basis.",
		      "resolve","OsiDylpSolverInterface") ; }
    resolveOptions->forcewarm = true ; }
/*
  Choose options appropriate for reoptimising and go to it. Phase is not
  crucial here --- dylp will sort it out as it starts.
*/
  assert(resolveOptions->forcecold == false) ;
  if (dyio_isactive(local_logchn)) dy_logchn = local_logchn ;
  dy_gtxecho = resolve_gtxecho ;
  dyphase_enum phase = lpprob->phase ;
  if (!(phase == dyPRIMAL1 || phase == dyPRIMAL2 || phase == dyDUAL))
    lpprob->phase = dyINV ;
  lp_retval = do_lp(startWarm) ;
/*
  Separate the failure cases from the successful cases. lpITERLIM is
  considered ok in warm and hot start (in particular, for strong branching).
  Since we can't tell from here, include it in the successes.
*/
# ifdef ODSI_INFOMSGS
  hdl->message(ODSI_WARM,messages_)
    << dy_prtlpret(lp_retval)
    << getObjSense()*(lpprob->obj) << lpprob->iters
    << CoinMessageEol ;
# endif
  bool lpOK ;
  if ((lp_retval == lpOPTIMAL || lp_retval == lpINFEAS ||
       lp_retval == lpUNBOUNDED || lp_retval == lpITERLIM))
  { lpOK = true ; }
  else
  { lpOK = false ; }
/*
  Trash the cached solution. This needs to happen regardless of how the lp
  turned out. It needs to be done ahead of getWarmStart, which will try to
  access reduced costs. We should only need to empty the solution vectors, not
  the structural vectors.
*/
  destruct_col_cache(false) ;
  destruct_row_cache(false) ;
/*
  Tidy up. If all went well, indicate this object owns the solver, set the
  objective, and set the active basis. If we've failed, do the opposite.
  Note that we have to remove any current active basis first, or getWarmStart
  will hand it back to us.

  Any of lpOPTIMAL, lpINFEAS, or lpUNBOUNDED can (in theory) be hot started
  after allowable problem modifications, hence can be flagged lpctlDYVALID.
  dylp overloads lpprob->obj with the index of the unbounded variable when
  returning lpUNBOUNDED, so we need to fake the objective.
*/
# if ODSI_TRACK_ACTIVE > 0
  std::cout
    << "ODSI(" << std::hex << this
    << ")::resolve: deleting active basis "
    << activeBasis << std::dec
    << "." << std::endl ;
# endif
  delete activeBasis ;
  activeBasis = 0 ;
  activeIsModified = false ;
  if (lpOK && flgon(lpprob->ctlopts,lpctlDYVALID))
  { dylp_owner = this ;
#   ifdef ODSI_INFOMSGS
    hdl->message(ODSI_ATTACH,messages_)
      << "resolve" << (int) reinterpret_cast<CoinIntPtr>(this)
      << CoinMessageEol ;
#   endif
    if (lpprob->lpret == lpUNBOUNDED)
    { _objval = -getObjSense()*getInfinity() ; }
    else
    { _objval = getObjSense()*lpprob->obj ; }
    activeBasis = this->getWarmStart() ;
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::resolve: setting active basis "
      << activeBasis << std::dec
      << "." << std::endl ;
#   endif
    resolveOptions->forcewarm = false ; }
  else
  { dylp_owner = 0 ; }
  
  return ; }

//@} // WarmStart



/*! \defgroup HotStartMethods Hot Start Methods
    \brief Methods for solving an LP from a hot start.

  The restrictions on problem modification outlined in the OSI documentation
  are somewhat relaxed in dylp: you can change upper or lower bounds on
  variables, upper or lower bounds on constraints, or objective function
  coefficients. No changes to the coefficient matrix are allowed --- for
  this, you need to drop back to a warm start. Philosophically, the guiding
  principle is to avoid the need to refactor the basis before resuming
  pivoting.

  To dylp, a hot start is simply a resumption of pivoting under the
  assumption that the data structures left from the previous call are
  unchanged, modulo the modifications permitted in the previous paragraph.
  There is no native hot start object.  This is not sufficient to provide the
  semantics required of markHotStart() and solveFromHotStart(). In
  particular, if several ODSI objects are sharing the dylp solver, they will
  interfere with one another.

  To hide this behaviour from the client, ODSI provides a fall back method
  along the lines described in the OSI documentation. markHotStart() captures
  a warm start object. If a hot start is possible, according to dylp's notion
  of hot start, this object will not be used.  If solveFromHotStart()
  determines that intervening use of the solver by some other ODSI object has
  made it impossible to perform a hot start for this object, it will
  automatically fall back to a warm start using the object captured by
  markHotStart().

  unmarkHotStart() simply deletes the warm start object captured by
  markHotStart().
*/
//@{

/*!
  This routine sets options in the ODSI object so that dylp will attempt a
  hot start at each call of
  \link OsiDylpSolverInterface::solveFromHotStart
	ODSI::solveFromHotStart \endlink.
*/

inline void ODSI::markHotStart ()
{ 
/*
  If we're attempting this, then it should be true that this OSI object
  owns the solver, has successfully solved an lp, and has left dylp's data
  structures intact. Should be, but who knows where dylp's been in the
  meantime. Applications written to a generic OSI interface don't tend to worry
  about these niceties. If we don't own the solver, try for a resolve.
*/
  if (dylp_owner != this)
  { resolve() ; }
  assert(lpprob && resolveOptions) ;
  assert(flgon(lpprob->ctlopts,lpctlDYVALID)) ;
/*
  Turn off the forcecold and forcewarm options, allowing a hot start.
*/
  resolveOptions->forcecold = false ;
  resolveOptions->forcewarm = false ;
/*
  Capture a warm start object, in case some other ODSI object uses the solver
  between hot starts by this ODSI object.
*/
  if (hotstart_fallback) delete hotstart_fallback ;
  hotstart_fallback = getWarmStart() ;

  return ; }

/*!
  Usually, this routine simply calls the solver. For dylp's model of hot
  start, that's all that's required. If some other ODSI object has grabbed
  the solver in the meantime, fall back to a warm start using the object
  captured by markHotStart().
*/

void ODSI::solveFromHotStart ()

{ 
/*
  If some other ODSI object has used the solver, then fall back to a warm
  start. Similarly, if the last call to dylp failed and the solver does not
  hold retained data structures, we need to fall back to a warm start.  Throw
  an error if there's no warm start object to fall back on.  (The throw
  message lies a little, for the benefit of users who haven't checked the
  details of the code.)
*/

  if (dylp_owner != this || flgoff(lpprob->ctlopts,lpctlDYVALID))
  { if (hotstart_fallback && setWarmStart(hotstart_fallback))
    { resolve() ; }
    else
    { throw CoinError("Hot start failed --- invalid/missing hot start object.",
		      "solveFromHotStart","OsiDylpSolverInterface") ; }
    return ; }
/*
  If no other ODSI object has used the solver, all we need to do is check the
  iteration limit. The basis should be ready.  Note that dylp does need to
  know what's changed: any of bounds, rhs & rhslow, or objective. The various
  routines that make these changes set the flags.
*/
  int tmp_iterlim = -1 ;
  int hotlim ;

  assert(lpprob && lpprob->basis && lpprob->status && basis_ready &&
	 consys && resolveOptions && tolerances) ;

  if (dyio_isactive(local_logchn)) dy_logchn = local_logchn ;
  dy_gtxecho = resolve_gtxecho ;
/*
  Phase can be anything except dyDONE, which will cause dylp to free data
  structures and return.
*/
  lpprob->phase = dyINV ;
  getIntParam(OsiMaxNumIterationHotStart,hotlim) ;
  if (hotlim > 0)
  { tmp_iterlim = resolveOptions->iterlim ;
    resolveOptions->iterlim = (hotlim/3 > 0)?hotlim/3:1 ; }

  lp_retval = do_lp(startHot) ;
# ifdef ODSI_INFOMSGS
  CoinMessageHandler *hdl = messageHandler() ; 
  hdl->message(ODSI_HOT,messages_)
    << dy_prtlpret(lp_retval)
    << getObjSense()*(lpprob->obj) << lpprob->iters
    << CoinMessageEol ;
# endif
/*
  Separate the failure cases from the successful cases. lpITERLIM is
  considered ok in warm and hot start (in particular, for strong branching).
  Since we can't tell from here, include it in the successes.
*/
  bool lpOK ;
  if ((lp_retval == lpOPTIMAL || lp_retval == lpINFEAS ||
       lp_retval == lpUNBOUNDED || lp_retval == lpITERLIM))
  { lpOK = true ; }
  else
  { lpOK = false ; }
/*
  Trash the cached solution. This needs to happen regardless of how the lp
  turned out. It needs to be done ahead of getWarmStart, which will try to
  access reduced costs. We should only need to empty the solution vectors, not
  the structural vectors.
*/
  destruct_col_cache(false) ;
  destruct_row_cache(false) ;
/*
  Tidy up. If all went well, set the objective and set the active basis (we
  already own the solver). If we've failed, do the opposite.  Note that we
  have to remove any current active basis first, or getWarmStart will hand it
  back to us.

  dylp overloads lpprob->obj with the index of the unbounded variable when
  returning lpUNBOUNDED, so we need to fake the objective.
*/
# if ODSI_TRACK_ACTIVE > 0
  std::cout
    << "ODSI(" << std::hex << this
    << ")::solveFromHotStart: deleting active basis "
    << activeBasis << std::dec
    << "." << std::endl ;
# endif
  delete activeBasis ;
  activeBasis = 0 ;
  activeIsModified = false ;
  if (lpOK && flgon(lpprob->ctlopts,lpctlDYVALID))
  { if (lpprob->lpret == lpUNBOUNDED)
    { _objval = -getObjSense()*getInfinity() ; }
    else
    { _objval = getObjSense()*lpprob->obj ; }
    activeBasis = this->getWarmStart() ;
#   if ODSI_TRACK_ACTIVE > 0
    std::cout
      << "ODSI(" << std::hex << this
      << ")::solveFromHotStart: setting active basis "
      << activeBasis << std::dec
      << "." << std::endl ;
#   endif
  }
  else
  { dylp_owner = 0 ; }
  if (tmp_iterlim > 0) resolveOptions->iterlim = tmp_iterlim ;

  return ; }


inline void ODSI::unmarkHotStart ()
/*
  Nothing to do here but to destroy the fall back warm start object.
*/

{ if (hotstart_fallback)
  { delete hotstart_fallback ;
    hotstart_fallback = 0 ; }

  return ; }

//@} // HotStartMethods



/*! \defgroup DebugMethods Debugging Methods */
//@{

/*!
  Because dylp retains information in the solver, the current solver state
  may need to be preserved around the activation of the row cut debugger.
  The problem arises because the row cut debugger clones the solver object
  and uses it to solve an lp; this destroys the retained state and transfers
  ownership of the solver to the clone. To preserve the solver state around
  debugger activation, we need to create a warm start object and then restore
  it after activating the debugger.

  activateRowCutDebugger tries to make an informed guess as to whether to
  preserve state. The condition is that some object owns the solver (not
  necessarily this object) and that the lpprob for that object exists and
  indicates preserved state in the solver. This is conservative --- it's more
  likely that the previous owner is no longer interested but was never
  detached.

  It's also conservative in the sense that unless the row cut debugger
  recognizes the model name, it won't activate.
  See activateRowCutDebugger(const double*) if you have a model the debugger
  doesn't already know about.
*/

void ODSI::activateRowCutDebugger (const char *modelName)

{ delete rowCutDebugger_ ;

  if (dylp_owner && dylp_owner->lpprob &&
      flgon(dylp_owner->lpprob->ctlopts,lpctlDYVALID))
  { CoinWarmStart *ws = dylp_owner->getWarmStart() ;
    ODSI *prev_owner = dylp_owner ;
    prev_owner->detach_dylp() ;
    rowCutDebugger_ = new OsiRowCutDebugger(*this,modelName) ;
    prev_owner->setWarmStart(ws) ;
    prev_owner->resolve() ;
    delete ws ; }
  else
  { rowCutDebugger_ = new OsiRowCutDebugger(*this,modelName) ; }

  return ; }

/*!
  Activate the row cut debugger for a model that's not one of the models known
  to the debugger. You must provide a full solution vector, but only the
  integer variables will be checked.

  See also the comments for activateRowCutDebugger(const char*).
*/

void ODSI::activateRowCutDebugger (const double *solution)

{ delete rowCutDebugger_ ;

  if (dylp_owner && dylp_owner->lpprob &&
      flgon(dylp_owner->lpprob->ctlopts,lpctlDYVALID))
  { CoinWarmStart *ws = dylp_owner->getWarmStart() ;
    ODSI *prev_owner = dylp_owner ;
    prev_owner->detach_dylp() ;
    rowCutDebugger_ = new OsiRowCutDebugger(*this,solution) ;
    prev_owner->setWarmStart(ws) ;
    prev_owner->resolve() ;
    delete ws ; }
  else
  { rowCutDebugger_ = new OsiRowCutDebugger(*this,solution) ; }

  return ; }


# if ODSI_PARANOIA >= 1
/*!
  This routine will check that the index passed as a parameter is a valid
  variable (isCol == true) or constraint (isCol == false) index and throw an
  error if it's out of bounds. The index is assumed to be zero-based, which is
  appropriate for indices specified by an ODSI client.

  This routine must be active (i.e., ODSI_PARANOIA must be 1 or greater) in
  order for OsiCbc(dylp) to pass the OsiCbc unit test.
*/
inline void ODSI::indexCheck (int k, bool isCol, std::string rtnnme)

{ std::string message ;

  if (!consys)
  { message = "No constraint system!" ;
    throw CoinError(message,rtnnme,"OsiDylpSolverInterface") ; }

  int m = getNumRows() ;
  int n = getNumCols() ;

  if (isCol)
  { if (0 > k || k > n)
    { message = "Column index out of range!" ;
      throw CoinError(message,rtnnme,"OsiDylpSolverInterface") ; } }
  else
  { if (0 > k || k > m)
    { message = "Row index out of range!" ;
      throw CoinError(message,rtnnme,"OsiDylpSolverInterface") ; } }
  
  return ; }

#endif

//@} // DebugMethods



/*! \defgroup DylpMethods Dylp-Specific Methods */
//@{

/*!
  Process a dylp options (.spc) file.

  \param name file name
  \param silent true to suppress command echo
  \param mustexist true (default) if the file must exist

  This is imperfect (to say the least) as one often wants different settings
  for the initial call to the solver vs. subsequent calls to reoptimize.
*/

void ODSI::dylp_controlfile (const char *name,
			     const bool silent, const bool mustexist)

{ if (name == 0 || *name == 0) return ;
  string mode = (mustexist)?"r":"q" ;
  dy_cmdchn = dyio_openfile(name,mode.c_str()) ;
  if (!(dy_cmdchn == IOID_INV || dy_cmdchn == IOID_NOSTRM))
  { dyio_setmode (dy_cmdchn, 'l') ;  
    main_lpopts = initialSolveOptions ;
    main_lptols = tolerances ;
    bool r UNUSED = (process_cmds(silent) != 0) ;
    (void) dyio_closefile(dy_cmdchn) ;
    dy_cmdchn = IOID_NOSTRM ;
    assert(r == cmdOK) ;
/*
  Copy user settings into resolveOptions, except for forcecold and fullsys.
*/
    lpopts_struct saveOptions = *resolveOptions ;
    memcpy(resolveOptions,initialSolveOptions,sizeof(lpopts_struct));
    resolveOptions->forcecold = saveOptions.forcecold ;
    resolveOptions->fullsys = saveOptions.fullsys ; }

  dy_cmdchn = IOID_NOSTRM ;

  return ; }
  
/*!
  Establish a log (.log) file.

  \param name file name
  \param echo true to echo log messages to the terminal
*/

void ODSI::dylp_logfile (const char *name, bool echo)

{ if (name == 0 || *name == 0) return ;

  string lognme = make_filename(name,".mps",".log") ;
  local_logchn = dyio_openfile(lognme.c_str(),"w") ;
  if (local_logchn == IOID_INV)
  { local_logchn = IOID_NOSTRM ; }
  else
  { (void) dyio_chgerrlog(lognme.c_str(),true) ; }
  initial_gtxecho = echo ;
  resolve_gtxecho = echo ;

  return ; }

/*!
  Establish an output (.out) file.

  \param name file name
*/

void ODSI::dylp_outfile (const char *name)

{ if (name == 0 || *name == 0) return ;

  string outnme = make_filename(name,".mps",".out") ;
  local_outchn = dyio_openfile(outnme.c_str(),"w") ;
  if (local_outchn == IOID_INV) local_outchn = IOID_NOSTRM ;
  
  return ; }

/*!
  Print the solution and/or statistics to the output file.

  \param wantSoln true to print the lp solution
  \param wantStats true to print execution statistics for the lp
*/

void ODSI::dylp_printsoln (bool wantSoln, bool wantStats)

{ if (dyio_isactive(local_outchn))
  { if (wantStats) dy_dumpstats(local_outchn,false,statistics,consys) ;
    if (wantSoln) dy_dumpcompact(local_outchn,false,lpprob,false) ; }

  return ; }

//@} // DylpMethods
