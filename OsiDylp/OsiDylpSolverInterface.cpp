/*! \legal
  Copyright (C) 2002, 2003
  Lou Hafer, Stephen Tse, International Business Machines Corporation and
  others. All Rights Reserved.
*/

#ifdef COIN_USE_DYLP

#ifdef _MSC_VER

/* Turn off compiler warning about long names */
#  pragma warning(disable:4786)

/*
  MS C++ doesn't want to cope with this set of templates. Since they're just
  paranoid checks, it's easiest to simply disable them.
*/

#  define assert_same

/*
  In general, templates do not seem to be a strong point for MS C++. To get
  around this, there are a number of macros defined below.  For Sun Workshop
  and Gnu programming environments, these will expand to template functions.
  For MS, they expand to some code plus calls to STL functions. Apparently MS
  C++ doesn't have trouble with the definitions of the templates, just
  instantiations.
*/

#endif

/* Cut name lengths for readability. */

#define ODSI OsiDylpSolverInterface
#define OSI OsiSolverInterface

/*!
   \file OsiDylpSolverInterface.cpp

   \brief Implementation of COIN OSI layer for dylp.

   This file contains the implementation of a COIN OSI layer for dylp, the lp
   solver for the bonsaiG MILP code.

   More information on the COIN/OR project and the OSI layer specification
   can be found at http://www.coin-or.org.

   -- Implementation Principles --

   <strong>Option and Tolerance Existence</strong>:
   These structures are created and initialised by the ODSI constructor,
   so they're valid from the moment the ODSI object is created. Option and
   tolerance settings are preserved when a new problem is loaded or assigned.

   <strong>Constraint System Existence</strong>:
   An invariant in the interface is that if getNumCols or getNumRows returns
   a value greater than 0, then a valid constraint system structure (consys)
   exists, and attached to it are valid vectors for the objective
   coefficients (obj), variable type (vtyp) and upper and lower bounds (vub,
   vlb), and constraint type (ctyp) and right-hand side (rhs, rhslow).

   <strong>Construction of a Constraint System</strong>:
   A constraint system can be batch loaded from a file (MPS format) or from
   a data structure, or it can be built incrementally. When building a
   constraint system incrementally, keep in mind that you must create a row
   or column (addRow or addCol, respectively) before you can adjust other
   properties (row or column bounds, objective, variable values, etc.)

   <strong>LP Problem Existence</strong>
   The LP problem structure is not created until it's needed. A valid LP
   problem structure will exist after the first call to the solver
   (initialSolve) or after a warm start object is loaded (setWarmStart).

   <strong>Caching</strong>:
   Since vectors returned by OsiSolverInterface "get" functions are constant
   (not modifiable by users) and it is convenient for users to make repeated
   calls for the same information, a cache mechanism is used when necessary
   to avoid repeated expensive conversions in the interface. All cache
   variables are prefixed with underscore (_col*, _row*, _matrix*).  In
   general, modifying the problem causes the cache to be invalidated.

   <strong>Index Base</strong>:
   OsiSolverInterface indexes variables and constraints from 0 (the standard
   approach for C/C++) while dylp indexes them from 1 (for robustness; see
   consys.h). Some caution is needed in the interface to avoid off-by-one
   errors. Similarly, arrays (vectors) for k elements in bonsai occupy k+1
   space with vector[0] unused.  ODSI::idx, ODSI::inv, and ODSI::inv_vec are
   used to make these trivial conversions explicit.

   <strong>Copy and Assert</strong>:
   To clone, each part of dylp's data structures (lpprob, consys, basis) is
   copied (ODSI::copy*) and then (optionally) verified (ODSI::assert_same*
   and assert). Helpers for copying and verifying primitive types (integer
   and double) and array types are implemented with C++ templates. (If you
   skipped it, this is a good moment to go back and read the note above about
   MS C++ and templates.)
*/

namespace {
  char sccsid[] = "%W%	%G%" ;
  char cvsid[] = "$Id$" ;
}

#include <string>
#include <cassert>
#include <OsiColCut.hpp>
#include <OsiRowCut.hpp>
#include "OsiDylpSolverInterface.hpp"
#include "OsiDylpMessages.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinMpsIO.hpp"

using std::string ;
using std::vector ;


extern "C"
{

#include "bonsai.h"

/*
  Dylp uses HUGE_VAL for infinity; generally this seems to resolve to IEEE
  infinity but if it doesn't, you'll need to have a look at vector.h for the
  appropriate compile-time symbol definition. A macro seems to be necessary
  here because I can't figure out how to call ODSI.getInfinity() from inside
  a static member function without explicitly passing in the value of `this'
  from the context of the call.  That seems like overkill.
*/

#define DYLP_INFINITY HUGE_VAL


#ifndef ERRMSGDIR
# define ERRMSGDIR "."
#endif


extern bool dy_mpsin(const char *filename, consys_struct **consys) ;
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
  mechanisms are not a good fit for C++ objects and the OSI framework.
  Four global variables, cmdchn, cmdecho, logchn, and gtxecho, control
  command input/echoing and log message output/echoing.

  Gtxecho can be indirectly controlled using the OsiDoReducePrint hint; it
  is set to true whenever messageHandler()->logLevel() > 0.

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
/*! \var ioid cmdchn
    \brief ioid used to read option (.spc) files
*/
/*! \var ioid logchn
    \brief ioid used for logging to a file
*/
/*! \var bool cmdecho
    \brief controls echoing of commands from option files to the terminal
*/
/*! \var bool gtxecho
    \brief controls echoing of generated text to the terminal
*/

ioid cmdchn = IOID_NOSTRM,
     logchn = IOID_NOSTRM ;

bool cmdecho = false,
     gtxecho = false ;

//@}

/*!
  \defgroup DylpResidual dylp residual control variables
  \brief Variables controlling persistent dylp subsystems

  These variables are used to control initialisation and shutdown of the i/o
  and basis maintenance subsystems.
  
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


/*! \brief Specialised copy function for a dylp basis_struct */
  
inline void ODSI::copy_basis (const basis_struct* src, basis_struct* dst)
/*
  We use CALLOC here to minimize mixing of memory allocation with new & alloc.
  This structure stands a good chance of being reallocated within dylp.
*/

{ if (!src) return ;

  dst->len = src->len ;
  if (dst->el != 0) FREE(dst->el) ;
  dst->el = (basisel_struct *) CALLOC(idx(src->len),sizeof(basisel_struct)) ;
  memcpy(dst->el,src->el,idx(src->len)*sizeof(double)) ;

# ifndef NDEBUG
  assert_same(*dst, *src, true) ;
# endif

  return ; }


/*! \brief Specialised copy function for a dylp basis_struct */

inline basis_struct* ODSI::copy_basis (const basis_struct* src)

{ if (!src) return 0 ;
  basis_struct* dst = new basis_struct ;
  dst->el = 0 ;
  copy_basis(src, dst) ;
  return dst ; }


/*! \brief Specialised copy function for the constraint matrix of a
	   consys_struct

  The algorithm creates an empty row for each constraint, then adds the
  coefficients column by column.

  \todo This could be done more efficiently as a primitive routine in
	consys_utils.c. The column additions aren't so bad, but the row
	additions are incurring unnecessary overhead.
*/


/*! \brief Specialised copy function for a dylp lpprob_struct */

lpprob_struct* ODSI::copy_lpprob (const lpprob_struct* src)

{ if (!src) return 0 ;

  int col_count = idx(src->colsze) ;
  int row_count = idx(src->rowsze) ;

  lpprob_struct* dst = NULL;
  CLONE(lpprob_struct,src,dst);

  dst->basis = copy_basis(src->basis) ;
  CLONE_VEC(flags,src->status,dst->status,col_count);
  CLONE_VEC(double,src->x,dst->x,row_count);
  CLONE_VEC(double,src->y,dst->y,row_count);
  CLONE_VEC(bool,src->actvars,dst->actvars,col_count);

# ifndef NDEBUG
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

  Destroy cached copies of computed values for columns.
  See the \link OsiDylpSolverInterface.cpp comments \endlink at the head of
  the file for a brief description of the ODSI cache mechanism.
*/

inline void ODSI::destruct_col_cache ()

{ delete [] _col_x ; _col_x = 0 ;
  delete [] _col_obj ; _col_obj = 0 ;
  delete [] _col_cbar ; _col_cbar = 0 ;
}


/*! \brief Destroy cached values
    \ingroup DestructorHelpers

  Destroy cached copies of computed values for rows.
  See the \link OsiDylpSolverInterface.cpp comments \endlink at the head of
  the file for a brief description of the ODSI cache mechanism.
*/

void ODSI::destruct_row_cache ()

{ delete [] _row_lhs ; _row_lhs = 0 ;
  delete [] _row_lower ; _row_lower = 0 ;
  delete [] _row_price ; _row_price = 0 ;
  delete [] _row_range ; _row_range = 0 ;
  delete [] _row_rhs ; _row_rhs = 0 ;
  delete [] _row_sense ; _row_sense = 0 ;
  delete [] _row_upper ; _row_upper = 0 ;
}

/*! \brief Destroy cached values
    \ingroup DestructorHelpers

  Destroy cached copies of computed values.
  See the \link OsiDylpSolverInterface.cpp file comments \endlink
  for a brief description of the ODSI cache mechanism.
*/
  
inline void ODSI::destruct_cache ()

{ destruct_row_cache() ;
  destruct_col_cache() ;

  delete _matrix_by_row ; _matrix_by_row = 0 ;
  delete _matrix_by_col ; _matrix_by_col = 0 ;

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

{ double inf = DYLP_INFINITY ;

  if (upper == lower) return contypEQ ;
  bool finite_low = (lower > -inf) ;
  bool finite_up = (upper < inf) ;
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
  Add the column.
*/
  bool r = consys_addcol_pk(consys,vtypj,pk_colj,objj,vlbj,vubj) ;
  pkvec_free(pk_colj) ;
  assert(r) ;
/*
  After adding a column, the best we can do is a warm start.
*/
  resolveOptions->forcewarm = true ;

  destruct_cache() ; }


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
  Add the row.
*/
  bool r = consys_addrow_pk(consys,clazzi,ctypi,pk_rowi,rhsi,rhslowi,0,0) ;
  pkvec_free(pk_rowi) ;
  assert(r) ;
/*
  After adding a constraint, the best we can do is a warm start.
*/
  resolveOptions->forcewarm = true ;

  destruct_cache() ; }


/*! \brief Establish a worst-case primal solution.

  This routine sets up a `worst case' primal solution, so that the objective
  reported before invoking the solver will be a worst case bound. Each primal
  variable is set to the bound which gives the worst objective value.

  The interface design guarantees that if a consys_struct exists, then
  the objective and bounds vectors are already attached to it.
*/

void ODSI::worst_case_primal ()

{ 

/*
  No columns, no solution. A non-zero value guarantees a consys structure.
*/
  int n = getNumCols() ;
  if (n == 0) return ;
  assert(consys && consys->obj && consys->vub && consys->vlb) ;

  double *obj = consys->obj ;
  double *vlb = consys->vlb ;
  double *vub = consys->vub ;
/*
  We have columns. Allocate a cached primal solution vector, if it doesn't
  already exist.
*/
  if (!_col_x) _col_x = new double[n] ;
/*
  Walk the objective vector, taking values for variables from the appropriate
  bounds vector. The objective attached to consys is always a minimisation
  objective.
*/
  for (int j = 1 ; j <= n ; j++)
  { if (obj[j] > 0)
    { _col_x[inv(j)] = vub[j] ; }
    else
    { _col_x[inv(j)] = vlb[j] ; } }
  
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

/*! \brief Construct a dylp lpprob_struct (LP problem).

    \todo Remove the call to dy_setprintopts, once I've figured out a better
	  scheme for handling dylp parameters.
*/

inline void ODSI::construct_lpprob ()

{
  dy_checkdefaults(consys,initialSolveOptions,tolerances) ;
  dy_checkdefaults(consys,resolveOptions,tolerances) ;
  lpprob = new lpprob_struct ;
  memset(lpprob, 0, sizeof(lpprob_struct)) ;
  setflg(lpprob->ctlopts,lpctlNOFREE) ;
  lpprob->phase = dyINV ;
  lpprob->consys = consys ;
  lpprob->rowsze = consys->rowsze ;
  lpprob->colsze = consys->colsze ;
}


/*! \brief Construct and initialize default options, tolerances, and
	   statistics structures
*/

inline void ODSI::construct_options ()

{ 
/*
  Acquire the default options and tolerances from dylp. Set the OSI options
  and tolerances to match. Turn dylp's printing options down to `catatonic'.
*/
  initialSolveOptions = new lpopts_struct ;
  tolerances = new lptols_struct ;
  dy_defaults(&initialSolveOptions,&tolerances) ;
  CLONE(lpopts_struct,initialSolveOptions,resolveOptions);
  dy_setprintopts(0,initialSolveOptions) ;
  dy_setprintopts(0,resolveOptions) ;

  setIntParam(OsiMaxNumIteration,3*initialSolveOptions->iterlim) ;
  setIntParam(OsiMaxNumIterationHotStart,3*initialSolveOptions->iterlim) ;
  setDblParam(OsiDualTolerance,tolerances->cost) ;
  setDblParam(OsiPrimalTolerance,tolerances->zero) ;
/*
  Differentiate options for initial solution vs. reoptimising.
*/
  initialSolveOptions->forcecold = true ;
  initialSolveOptions->fullsys = true ;
/*
  Acquire and clear a statistics structure.
*/
  statistics = new lpstats_struct ;
  memset(statistics,0,sizeof(lpstats_struct)) ;
}

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
  consys = consys_create(0,parts,opts,rows,cols) ;
  assert(consys) ;

  return ; }


/*! \brief Initialise dylp i/o subsystems

  Initialise dylp i/o subsystems for normal and error output. This should be
  done only once.
*/

void ODSI::dylp_ioinit ()

{ if (reference_count > 1) return ;

  std::string errfile = std::string(ERRMSGDIR)+std::string("/bonsaierrs.txt") ;
# ifndef NDEBUG
  errinit(const_cast<char *>(errfile.c_str()),0,true) ;
# else
  errinit(const_cast<char *>(errfile.c_str()),0,false) ;
# endif
  bool r1 = ioinit() ;
  assert(r1) ;
}

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
  Free the existing problem structures, preserving options and tolerances
  and resetting statistics.
*/
  destruct_problem(true) ;
/*
  Create an empty consys_struct.
*/
  int m = matrix.getNumCols() ;
  int n = matrix.getNumRows() ;

  construct_consys(n,m) ;
/*
  First loop: Insert empty constraints into the new constraint system.
*/
  pkvec_struct* rowi = pkvec_new(0) ;
  assert(rowi) ;

  for (int i = 0 ; i < n ; i++)
  { rowi->nme = 0 ;
    bool r = consys_addrow_pk(consys,'a',ctyp[i],rowi,rhs[i],rhslow[i],0,0) ;
    assert(r) ;
  }

  if (rowi) pkvec_free(rowi) ;
/*
  Second loop. Insert the coefficients by column. If we need to create a
  column-ordered copy, take advantage and cache it.  The size of colj is a
  gross overestimate, but it saves the trouble of constantly reallocating it.
*/
  const CoinPackedMatrix *matrix2 ;
  if (matrix.isColOrdered())
  { matrix2 = &matrix ; }
  else
  { _matrix_by_col = new CoinPackedMatrix ;
    _matrix_by_col->reverseOrderedCopyOf(matrix) ;
    matrix2 = _matrix_by_col ; }
  
  pkvec_struct* colj = pkvec_new(n) ;

  for (int j = 0 ; j < m ; j++)
  { const CoinShallowPackedVector coin_col = matrix2->getVector(j) ;
    packed_vector(coin_col,n,colj) ;
    double objj = obj?obj[j]:0 ;
    double vlbj = col_lower?col_lower[j]:0 ;
    double vubj = col_upper?col_upper[j]:DYLP_INFINITY ;
    colj->nme = 0 ;
    bool r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    assert(r) ; }

  pkvec_free(colj) ;
  assert(matrix2->isEquivalent(*getMatrixByCol())) ;

/*
  Finish up. Establish a primal solution and an objective. The primal solution
  is `worst-case', in the sense that it's constructed by setting each primal
  variable to the bound that maximises the objective.
*/
  worst_case_primal() ;
  calc_objval() ;

  return ; }


/*! \brief Load a problem description

  This routine expects a constraint system described using a standard column-
  major packed matrix structure and dylp constraint sense and right-hand-side
  (rhs, rhslow) vectors.

  The matrix description consists of row and column size and three arrays.
  Coefficients and corresponding row indices are given in value and index,
  respectively. The vector start contains the starting index in (value,
  index) for each column.

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
	   const int *start, const int *index, const double *value,
	   const double* col_lower, const double* col_upper, const double* obj,
	   const contyp_enum *ctyp, const double* rhs, const double* rhslow)

{ 
/*
  Free the existing problem structures, preserving options and tolerances
  and resetting statistics.
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

  for (int i = 0 ; i < rowcnt ; i++)
  { rowi->nme = 0 ;
    bool r = consys_addrow_pk(consys,'a',ctyp[i],rowi,rhs[i],rhslow[i],0,0) ;
    assert(r) ;
  }

  if (rowi) pkvec_free(rowi) ;
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
    int lenj = start[j+1]-startj ;

    for (int ndx = 0 ; ndx < lenj ; ndx++)
    { coeffs[ndx].ndx = idx(index[startj+ndx]) ;
      coeffs[ndx].val = value[startj+ndx] ; }
    colj->cnt = lenj ;
  
    double objj = obj?obj[j]:0 ;
    double vlbj = col_lower?col_lower[j]:0 ;
    double vubj = col_upper?col_upper[j]:DYLP_INFINITY ;
    colj->nme = 0 ;
    bool r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    assert(r) ;
  }

  if (colj) pkvec_free(colj) ;
/*
  Finish up. Establish a primal solution and an objective.
*/
  worst_case_primal() ;
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

{ double inf = DYLP_INFINITY ;
  double rowlbi,rowubi ;

  for (int i = 0 ; i < rowcnt ; i++)
  { rowlbi = rowlb?rowlb[i]:-inf ;
    rowubi = rowub?rowub[i]:inf ;
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
      { rhs[i] = inf ;
	rhslow[i] = -inf ;
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
  will be instantiated and initialized, but not the lp problem or constraint
  system structures.

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

    local_logchn(IOID_NOSTRM),
    initial_gtxecho(false),
    resolve_gtxecho(false),
    lp_retval(lpINV),
    obj_sense(1.0),

    solvername("dylp"),
    mps_debug(0),

    _col_x(0),
    _col_obj(0),
    _col_cbar(0),
    _row_lhs(0),
    _row_lower(0),
    _row_price(0),
    _row_range(0),
    _row_rhs(0),
    _row_sense(0),
    _row_upper(0),
    _matrix_by_row(0),
    _matrix_by_col(0)

{
/*
  Replace the OSI default messages with ODSI messages.
*/
  setOsiDylpMessages(CoinMessages::us_en) ;
/*
  Clear the hint info_ array.
*/
  for (int i = 0 ; i < OsiLastHintParam ; i++) info_[i] = 0 ;
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
    assert(eq(DYLP_INFINITY, DYLP_INFINITY)) ; }
  
  return ; }


/*!
  A true copy --- no data structures are shared with the original. Cached
  information is also replicated.
*/

ODSI::OsiDylpSolverInterface (const OsiDylpSolverInterface& src)

  : OsiSolverInterface(src),
    initial_gtxecho(src.initial_gtxecho),
    resolve_gtxecho(src.resolve_gtxecho),
    lp_retval(src.lp_retval),
    obj_sense(src.obj_sense),

    solvername(src.solvername),
    mps_debug(src.mps_debug),
  
    _matrix_by_row(0),
    _matrix_by_col(0)

{ if (src.consys)
  { bool r = consys_dupsys(src.consys,&consys,src.consys->parts) ;
    assert(r) ; }
  else
  { consys = 0 ; }
  if (src.lpprob)
  { lpprob = copy_lpprob(src.lpprob) ;
    lpprob->consys = consys ; }
  else
  { lpprob = 0 ; }

  CLONE(lpopts_struct,src.initialSolveOptions,initialSolveOptions) ;
  CLONE(lpopts_struct,src.resolveOptions,resolveOptions) ;
  CLONE(lpstats_struct,src.statistics,statistics) ;
  CLONE(lptols_struct,src.tolerances,tolerances) ;

  if (src.local_logchn != IOID_NOSTRM)
  { local_logchn = openfile(idtopath(src.local_logchn),
			    const_cast<char *>("w")) ; }
  else
  { local_logchn = IOID_NOSTRM ; }
  assert(local_logchn == src.local_logchn) ;

  int col_count = src.getNumCols() ;
  int row_count = src.getNumRows() ;

  CLONE_VEC(double,src._col_x,_col_x,col_count) ;
  CLONE_VEC(double,src._col_obj,_col_obj,col_count) ;
  CLONE_VEC(double,src._col_cbar,_col_cbar,col_count) ;
  CLONE_VEC(double,src._row_lhs,_row_lhs,row_count) ;
  CLONE_VEC(double,src._row_lower,_row_lower,row_count) ;
  CLONE_VEC(double,src._row_price,_row_price,row_count) ;
  CLONE_VEC(double,src._row_range,_row_range,row_count) ;
  CLONE_VEC(double,src._row_rhs,_row_rhs,row_count) ;
  CLONE_VEC(char,src._row_sense,_row_sense,row_count) ;
  CLONE_VEC(double,src._row_upper,_row_upper,row_count) ;

  if (src._matrix_by_row) 
    _matrix_by_row = new CoinPackedMatrix(*src._matrix_by_row) ;
  if (src._matrix_by_col) 
    _matrix_by_col = new CoinPackedMatrix(*src._matrix_by_col) ;

  COPY_VEC(void *,&src.info_[0],&info_[0],OsiLastHintParam) ;

  reference_count++ ;

# ifndef NDEBUG
  assert_same(*this, src, true) ;
# endif

}


/*!
  An alternate copy constructor; the parameter is ignored (and assumed true).
*/

inline OsiSolverInterface* ODSI::clone (bool copyData) const

{ if (copyData)
  { return new OsiDylpSolverInterface(*this) ; }
  else
  { return new OsiDylpSolverInterface() ; }
}


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

  If the base ODSI object will be retained, set preserve_interface to true.
  This will retain the current options and tolerance settings, and reinitialise
  the statistics structure.
*/

void ODSI::destruct_problem (bool preserve_interface)

{ 
/*
  If this object claims ownership of the solver, it should possess an lpprob
  structure created when the solver was called.
*/
  assert((dylp_owner != this) || (dylp_owner == this && lpprob)) ;

  if (dylp_owner == this) detach_dylp() ;

  if (lpprob)
  { assert(lpprob->consys == consys) ;
    consys_free(consys) ;
    consys = 0 ;
    dy_freesoln(lpprob) ;
    delete lpprob ;
    lpprob = 0 ; }
  else
  if (consys)
  { consys_free(consys) ;
    consys = 0 ; }

  destruct_cache() ;

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
  
  if (statistics)
  { if (preserve_interface == true)
    { memset(statistics,0,sizeof(lpstats_struct)) ; }
    else
    { delete statistics ;
      statistics = 0 ; } }
  
  return ; }


/*! \brief Detach an ODSI instance from the dylp solver

  Dylp retains static data structures from one call to the next, in
  anticipation of reoptimizing following problem modification. When an ODSI
  instance is destroyed or another instance needs to use dylp, the current
  instance must be detached. Setting lpctlONLYFREE and phase == dyDONE tells
  dylp the call is solely to free data structures.
*/

void ODSI::detach_dylp ()

{ assert(dylp_owner == this && lpprob && lpprob->consys) ;
  assert(flgon(lpprob->ctlopts,lpctlDYVALID)) ;
  clrflg(lpprob->ctlopts,lpctlNOFREE) ;
  setflg(lpprob->ctlopts,lpctlONLYFREE) ;
  lpprob->phase = dyDONE ;
/*
  Either of the initialSolve or resolve option blocks are ok here; dylp does
  not look.
*/
  dylp(lpprob,initialSolveOptions,tolerances,statistics) ;
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
  file, if open.
*/
  destruct_problem(false) ;
  if (local_logchn != IOID_INV && local_logchn != IOID_NOSTRM)
    (void) closefile(local_logchn) ;

  reference_count-- ;
  if (reference_count == 0)
  { if (basis_ready == true)
    { dy_freebasis() ;
      basis_ready = false ; }
    ioterm() ;
    errterm() ; }
  
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
  destruct_problem(false) ;
  if (local_logchn != IOID_INV && local_logchn != IOID_NOSTRM)
  { (void) closefile(local_logchn) ;
    local_logchn = IOID_NOSTRM ; }
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
  obj_sense = 1.0 ;
  mps_debug = false ;
  construct_options() ;
  setOsiDylpMessages(CoinMessages::us_en) ;
  for (int i = 0 ; i < OsiLastHintParam ; i++) info_[i] = 0 ;

  return ; }



/*! \defgroup ProbAdjust Functions for Problem Modification
    \brief Public functions for problem modification

  This group of functions comprise the public interface for adjusting the
  optimization problem seen by the dylp solver.
*/
//@{

inline void ODSI::setContinuous (int j)

{ if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    assert(r) ; }

  consys->vtyp[idx(j)] = vartypCON ; }


inline void ODSI::setContinuous (const int* indices, int len)

{ if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    assert(r) ; }

  for (int i = 0 ; i < len ; i++)
    consys->vtyp[idx(indices[i])] = vartypCON ; }


inline void ODSI::setInteger (int j)

{ if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    assert(r) ; }

  if (getColLower()[j] == 0.0 && getColUpper()[j] == 1.0)
    consys->vtyp[idx(j)] = vartypBIN ;
  else
    consys->vtyp[idx(j)] = vartypINT ; }


inline void ODSI::setInteger (const int* indices, int len)

{ if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    assert(r) ; }

  for (int i = 0 ; i < len ; i++) setInteger(indices[i]) ; }



inline void ODSI::setColLower (int i, double val)

{ if (!consys->vlb)
  { bool r = consys_attach(consys,CONSYS_VLB,sizeof(double),
			   reinterpret_cast<void **>(&consys->vlb)) ;
    assert(r) ; }

  consys->vlb[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlLBNDCHG) ;
  if (isInteger(i)) setInteger(i) ; }


inline void ODSI::setColUpper (int i, double val)

{ if (!consys->vub)
  { bool r = consys_attach(consys,CONSYS_VUB,sizeof(double),
			   reinterpret_cast<void **>(&consys->vub)) ;
    assert(r) ; }

  consys->vub[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlUBNDCHG) ;
  if (isInteger(i)) setInteger(i) ; }


/*!
  A call to this routine destroys all cached row values.
*/

void ODSI::setRowType (int i, char sense, double rhs, double range)

{ int k = idx(i) ;

  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		sense,rhs,range) ;
  if (resolveOptions) resolveOptions->forcewarm = true ;

  destruct_row_cache() ; }


/*!
  A call to this routine destroys all cached row values.
*/

void ODSI::setRowUpper (int i, double val)

{ int k = idx(i) ;
  double clbi ;

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
    { clbi = -DYLP_INFINITY ;
      break ; }
    default:
    { assert(false) ; } }
  
  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		clbi,val) ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlRHSCHG) ;

  destruct_row_cache() ; }


/*!
  A call to this routine destroys all cached row values.
*/

void ODSI::setRowLower (int i, double val)

{ int k = idx(i) ;
  double cubi ;
  contyp_enum ctypi = consys->ctyp[k] ;

  if (ctypi == contypGE || ctypi == contypNB)
    cubi = DYLP_INFINITY ;
  else
    cubi = consys->rhs[k] ;

  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		val,cubi) ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlRHSCHG) ;

  destruct_row_cache() ; }


/*!
  Add a row to the constraint system given the coefficients and upper and lower
  bounds on the left-hand-side.

  A call to this routine destroys all cached values.
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

  A call to this routine destroys all cached values.
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

  A call to this routine destroys all cached values.
*/

void ODSI::applyRowCut (const OsiRowCut &cut)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,cut.lb(),cut.ub()) ;
  add_row(cut.row(),'c',ctypi,rhsi,rhslowi) ;
  
  return ; }


/*!
  A call to this routine destroys all cached values.
*/

void ODSI::deleteRows (int count, const int* rows)

{ for (int k = 0 ; k < count ; k++)
  { int i = idx(rows[k]) ;
    bool r = consys_delrow_stable(consys,i) ;
    assert(r) ; }
  if (resolveOptions) resolveOptions->forcecold = true ;

  destruct_cache() ; }



/*!
  Change a coefficient in the objective function.
*/

void ODSI::setObjCoeff (int j, double objj)

{ if (j >= getNumCols()) return ;

  consys->obj[idx(j)] = obj_sense*objj ;
  if (_col_obj) _col_obj[j] = objj ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ;
  calc_objval() ; }
  

/*!
  Change the sense of the objective; use 1.0 to request the objective be
  minimised, -1.0 to request it be maximised.
*/

void ODSI::setObjSense (double val)
/*
  The `natural' action of OSI (and dylp) is minimisation. Maximisation is
  accomplished as min -cx. So using -1 for maximisation is more natural
  than you'd think at first glance.
*/

{ int n = getNumCols() ;

  if (n > 0 && val != obj_sense)
  { double *tmpobj = INV_VEC(double,consys->obj) ;
    std::transform(tmpobj,tmpobj+n,tmpobj,std::negate<double>()) ;
    if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ; }
  
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
*/

void ODSI::deleteCols (int count, const int* cols)

{ for (int i = 0 ; i < count ; i++)
  { int j = idx(cols[i]) ;
    bool r = consys_delcol(consys, j) ;
    assert(r) ; }
  if (resolveOptions) resolveOptions->forcecold = true ;

  destruct_cache() ; }


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
  when comparing floating point values. The inexact test looks for equality
  within a tolerance of 1.0e-10, scaled by the maximum of the two values.  To
  do a toleranced test, both values must be finite.
*/

void ODSI::assert_same (double d1, double d2, bool exact)

{ if (d1 == d2) return ;

  assert(!exact && finite(d1) && finite(d2)) ;

  static const double epsilon = 1.e-10 ;
  double tol = std::max(fabs(d1),fabs(d2))+1 ;
  double diff = fabs(d1 - d2) ;

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

  int size = b1.len*sizeof(basisel_struct) ;
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

  Compare two bonsai consys_struct's for equality. There are two levels of
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
  assert(c1.maxcolndx == c2.maxcolndx) ;
  assert(c1.concnt == c2.concnt) ;
  assert(c1.archccnt == c2.archccnt) ;
  assert(c1.cutccnt == c2.cutccnt) ;
  assert(c1.maxrowlen == c2.maxrowlen) ;
  assert(c1.maxrowndx == c2.maxrowndx) ;
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

  Verify equivalence by checking each component in turn. Rely on the
  CoinPackedMatrix equivalence check to verify equivalence of the cached
  copies.
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
  assert(o1.local_logchn == o2.local_logchn) ;
  assert(o1.initial_gtxecho == o2.initial_gtxecho) ;
  assert(o1.resolve_gtxecho == o2.resolve_gtxecho) ;
  assert(o1.lp_retval == o2.lp_retval) ;
  assert(o1.obj_sense == o2.obj_sense) ;
  assert(o1.solvername == o2.solvername) ;
  assert(o1.mps_debug == o2.mps_debug) ;
/*
  Options, tolerances, and statistics should be byte-for-byte identical.
*/
  assert_same(*o1.initialSolveOptions, *o2.initialSolveOptions, true) ;
  assert_same(*o1.resolveOptions, *o2.resolveOptions, true) ;
  assert_same(*o1.statistics, *o2.statistics, true) ;
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
  const CoinPackedMatrix* m2 = o2.getMatrixByCol() ;
  if (m1)
  { assert(m2 || m1->isEquivalent(*m2)) ; }
  else
  { assert(!m2) ; }
}

//@} // CopyVerifiers

#endif /* !_MSC_VER */



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

std::string ODSI::make_filename (const char *filename,
				 const char *ext1, const char *ext2)

{ std::string basename(filename) ;
  std::string ext1str(ext1), ext2str(ext2) ;

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

   The readMps, writeMps, and writeMpsNative routines extend the default OSI
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

  Before reading the MPS problem file, readMps looks for a dylp control file
  (extension "spc"). The name is constructed by stripping the extension 
  (if present) from the given file name and adding the extension "spc".

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
    ext:	file extension; defaults to "mps"

  Returns: -1 if the MPS file could not be opened, otherwise the number of
	   errors encountered while reading the file.
*/

{ int errcnt ;
  CoinMpsIO mps ;
  CoinMessageHandler *mps_handler = mps.messageHandler() ;

  if (mps_debug)
  { mps_handler->setLogLevel(handler_->logLevel()) ; }
  else
  { mps_handler->setLogLevel(0) ; }

/*
  See if there's an associated .spc file and read it first, if present.
*/
  std::string filename = make_filename(basename,ext,"spc") ;
  dylp_controlfile(filename.c_str(),true,false) ;
/*
  Make sure the MPS reader has the same idea of infinity as dylp, then
  attempt to read the MPS file.
*/
  mps.setInfinity(DYLP_INFINITY) ;
  filename = make_filename(basename,ext,ext) ;
  errcnt = mps.readMps(filename.c_str(),0) ;
  handler_->message(ODSI_MPSFILEIO,messages_)
    << filename << "read" << errcnt << CoinMessageEol ;
  if (errcnt != 0) return (errcnt) ;
/*
  Transfer the major data vectors into a consys_struct and build a worst
  case solution.
*/
  loadProblem(*mps.getMatrixByCol(),mps.getColLower(),mps.getColUpper(),
	       mps.getObjCoefficients(),mps.getRowSense(),
	       mps.getRightHandSide(),mps.getRowRange()) ;
/*
  All those little bits and pieces. First, the problem and objective names.
  Then acquire any constant component for the objective.
*/
  setStrParam(OsiProbName,mps.getProblemName()) ;
  if (consys->objnme) strfree(consys->objnme) ;
  consys->objnme = stralloc(mps.getObjectiveName()) ;
  setDblParam(OsiObjOffset,mps.objectiveOffset()) ;
/*
  Next, the integer variables. It's a bit of a pain that CoinMpsIO doesn't
  distinguish binary from general integer, but we can do it here by checking
  bounds. The test is drawn from OSI::isBinary and sucks in variables fixed
  at 0 or 1. There's a separate routine, OSI::isFreeBinary, which excludes
  variables fixed at 0 or 1.
*/
  const char *const intvars = mps.integerColumns() ;
  if (intvars)
  { int j ;
    int n = mps.getNumCols() ;
    vartyp_enum *const vtyp = consys->vtyp ;
    const double *const vlb = mps.getColLower() ;
    const double *const vub = mps.getColUpper() ;

    for (j = 0 ; j < n ; j++)
    { if (intvars[j])
      { if ((vlb[j] == 0.0 || vlb[j] == 1.0) &&
	    (vub[j] == 0.0 || vub[j] == 1.0))
	{ vtyp[idx(j)] = vartypBIN ;
	  consys->binvcnt++ ; }
	else
	{ vtyp[idx(j)] = vartypINT ;
	  consys->intvcnt++ ; } } } }

/*
  Success!
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

  if (objsense == getObjSense())
  { outputobj = const_cast<double *>(solverobj) ; }
  else
  { outputobj = new double[n] ;
    std::transform(solverobj,solverobj+n,outputobj,std::negate<double>()) ; }

  char *vartyp = new char[n] ;
  typedef char *charp ;
  char **colnames = new charp[n],
       **rownames = new charp[m] ;
  int i,j ;

  for (j = 0 ; j < n ; j++) vartyp[j] = isInteger(j) ;

  for (i = 0 ; i < m ; i++)
    rownames[i] = consys_nme(consys,'c',idx(i),false,0) ;
  
  for (j = 0 ; j < n ; j++)
    colnames[j] = consys_nme(consys,'v',idx(j),false,0) ;

  mps.setMpsData(*getMatrixByRow(),DYLP_INFINITY,
		 getColLower(),getColUpper(),outputobj,vartyp,
		 getRowLower(),getRowUpper(),colnames,rownames) ;
/*
  We really need to work on symbolic names for these magic numbers.
*/
  int errcnt = mps.writeMps(filename.c_str(),0,4,2) ;
  handler_->message(ODSI_MPSFILEIO,messages_)
    << filename << "read" << errcnt << CoinMessageEol ;
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
  load_problem(colcnt,rowcnt,start,index,value,
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
  load_problem(colcnt,rowcnt,start,index,value,
	       col_lower,col_upper,obj,ctyp,rhs,rhslow) ;

  delete[] rhs ;
  delete[] rhslow ;
  delete[] ctyp ; }

//@} // ODSILoadProb





/*! \defgroup SolverParms Methods to Set/Get Solver Parameters

  Dylp supports a limited set of the OSI parameters. Specifically,
  <dl>
    <dt>OsiDualTolerance</dt>
    <dd>(supported)
	This doesn't really correspond to a dylp parameter, but
	it can be faked (sort of).
	Dylp's dual feasibility tolerance is dynamically
	scaled from an absolute dual zero tolerance. There is also an
	additional scaling factor (`lpcontrol dfeas' in the dylp
	documentation, normally defaulted to 1.0) that can be controlled by
	the user.  Given a value for the OsiDualTolerance parameter, the code
	calculates the appropriate value for dfeas such that
	OsiDualTolerance = dfeas*(absolute dual zero tolerance).
    </dd>

    <dt>OsiDualObjectiveLimit</dt>
    <dd>(unsupported)
	The dynamic simplex algorithm used in dylp works with a partial
	constraint system, adding and deleting constraints and variables
	as required and generally maintaining a minimal active constraint
	system.
	A side effect is that neither the dual or primal objective changes
	monotonically, hence it's not really possible to support a limit
	on the dual objective.
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
    <dd>(unsupported)
	As for OsiDualObjectiveLimit.
    </dd>

    <dt>OsiPrimalTolerance</dt>
    <dd>(supported) Handled in the same manner as OsiDualTolerance.</dd>
  </dl>
*/

//@{

inline double ODSI::getInfinity () const { return DYLP_INFINITY ; }


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
    { retval = false ;
      break ; }
    default:
    { retval = false ;
      break ; } }
  
  return (retval) ; }


bool ODSI::getDblParam (OsiDblParam key, double& value) const
/*
  This simply duplicates OSI::getDblParam. I've kept it here until I decide
  what to do with the unsupported parameters Osi[Primal,Dual]ObjectiveLimit.
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
    { OSI::getDblParam(key,value) ;
      retval = false ;
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


bool ODSI::setStrParam (OsiStrParam key, const std::string& value)
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

bool ODSI::getStrParam (OsiStrParam key, std::string& value) const
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
  { std::string message("Dylp ") ;
    if (dylpSense = true)
    { message += "cannot disable " ; }
    else
    { message += "does not support " ; }
    message += msgString ;
    if (hintStrength == OsiForceDo)
    { handler_->message(ODSI_UNSUPFORCEDO,messages_)
        << message << CoinMessageEol ;
      throw CoinError(message,"setHintParam","OsiDylpSolverInterface") ; }
    else
    { handler_->message(ODSI_IGNORED,messages_)
	<< message << CoinMessageEol ; } }
  
  return ; }


/*!
  Since dylp really only looks at the options and tolerances structures for
  the solver instance, the approach taken for hints is to stash them away and
  make the appropriate modifications in the options and tolerances
  structures.  Hints are therefore persistent across all calls to the solver,
  until explicitly changed.

  Note that the client retains ownership of the blocks specified by the
  (void *) arguments. If ODSI took ownership, it would need to know how to
  allocate and free the block, and that way lies madness.
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
  Dylp has no presolve capabilities.
*/
    case OsiDoPresolveInInitial:
    case OsiDoPresolveInResolve:
    { unimp_hint(false,sense,strength,"presolve") ;
      retval = true ;
      break ; }
/*
  Dylp can suppress the dual, but cannot suppress the primal. Requesting
  OsiDoDual* with (true, OsiHint*) is ignored; (true,OsiForceDo) will throw
  an exception. On the other hand, suppressing the dual is generally not a
  good idea, so the hint is ignored at (true,OsiHintTry).
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
  Dylp knows only one crash, which preferably uses slacks, then architecturals,
  and artificials as a last resort.
*/
    case OsiDoCrash:
    { unimp_hint(true,sense,strength,"basis crash") ;
      retval = true ;
      break ; }
/*
  An unadorned hint changes the level by +/- 1, 2, or 3, depending on sense
  and strength (abusing the equivalence of enums and integers).  If info is
  non-zero, it is taken as a pointer to an integer which is the absolute
  value for the log level. In this case, sense is irrelevant. Note that
  CoinMessageHandler recognizes levels 0 -- 4 as generic levels.
  Conveniently, dylp has a similar set of generic levels. The larger powers
  of 2 (8, 16, 32, etc.) left to trigger specific classes of debugging
  printout. These must be set using info. Currently:
    0x08	CoinMpsIO messages
*/
    case OsiDoReducePrint:
    { int verbosity = messageHandler()->logLevel() ;
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

      dy_setprintopts(0,initialSolveOptions) ;
      dy_setprintopts(0,resolveOptions) ;
      if (verbosity == 0)
      { initial_gtxecho = false ;
	resolve_gtxecho = false ; }
      else
      { initial_gtxecho = true ;
	resolve_gtxecho = true ; }
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

  info = info_[key] ;

  return (true) ; }


//@} // SolverParms



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
  appropriate values for the initial solution of the lp and calls the solver.

  Forcing fullsys here is a good choice in terms of the performance of the
  code, but it does reduce flexibility. This should be made into a hint.
*/
void ODSI::initialSolve ()

{ 
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
  Is the basis package initialised?
*/
  if (basis_ready == false)
  { int count = static_cast<int>((1.5*getNumRows()) + 2*getNumCols()) ;
    dy_initbasis(count,initialSolveOptions->factor,0) ;
    basis_ready = true ; }
/*
  Does some other ODSI object own the solver? If so, detach it.
*/
  if (dylp_owner != 0 && dylp_owner != this) dylp_owner->detach_dylp() ;
/*
  Choose options appropriate for the initial solution of the lp and go to it.
*/
  if (isactive(local_logchn)) logchn = local_logchn ;
  gtxecho = initial_gtxecho ;
  lpprob->phase = dyINV ;
  lp_retval = dylp_dolp(lpprob,initialSolveOptions,tolerances,statistics) ;
  if (flgon(lpprob->ctlopts,lpctlDYVALID))
  { dylp_owner = this ;
    _objval = obj_sense*lpprob->obj ; }
  destruct_row_cache() ;
  destruct_col_cache() ;

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
  Aka primal unbounded.
*/

inline bool ODSI::isProvenDualInfeasible () const

{ return (lp_retval == lpUNBOUNDED) ; }


/*!
  Returns true if dylp abandoned the problem. This could be due to numerical
  problems (accuracy check or singular basis), stalling, unexpected loss of
  feasibility, inability to allocate space, or other fatal internal error.
*/

inline bool ODSI::isAbandoned () const

{ if (lp_retval == lpACCCHK || lp_retval == lpSINGULAR ||
      lp_retval == lpSTALLED || lp_retval == lpNOSPACE ||
      lp_retval == lpLOSTFEAS || lp_retval == lpPUNT || lp_retval == lpFATAL)
    return (true) ;
  else
    return (false) ; }


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
  We could have a cached solution, either supplied by the client or a
  previous solution from dylp. If not, we have to build one.  Dylp returns a
  vector x of basic variables, in basis order, and a status vector. To
  assemble a complete primal solution, it's necessary to construct one vector
  that merges the two sources. The result is cached.

  NOTE the caution above about superbasic (vstatSB) variables. Nonbasic free
  variables (vstatNBFR) occur under the same termination conditions, but are
  legitimately zero. We really should return NaN for superbasics, using the
  function <limits>:std::numeric_limits::signaling_NaN(), but it's just not
  worth the hassle. Sun Workshop C++ (Rogue Wave Software) has it, but Gnu C++
  prior to v3 doesn't even have <limits>, and the local v3 installation's
  <limits> points to another, non-existent header.
*/
{ if (_col_x) return _col_x ;

  if (lpprob)
  { assert(lpprob->status && lpprob->x) ;

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
	  { _col_x[j] = -DYLP_INFINITY ;
	    break ; } } } } }

  return (_col_x) ; }


/*!
  If necessary, a cached copy is built from the solution returned by the
  solver.
*/

const double* ODSI::getRowPrice () const
/*
  Dylp reports dual variables in basis order, so we need to write a vector
  with the duals dropped into position by row index.
*/

{ if (_row_price) return (_row_price) ;

  if (lpprob)
  { assert(lpprob->basis && lpprob->y) ;

    int m = getNumRows() ;
    _row_price = new double[m] ;
    basis_struct* basis = lpprob->basis ;

    memset(_row_price,0,m*sizeof(double)) ;

    for (int k = 1 ; k <= basis->len ; k++)
    { int i = inv(basis->el[k].cndx) ;
      _row_price[i] = lpprob->y[k]*obj_sense ; } }

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
  variables.
*/
    pkvec_struct *aj = pkvec_new(m) ;
    for (int j = 0 ; j < consys->varcnt ; j++)
    { if (x[j] != 0)
      { bool r = consys_getcol_pk(consys,idx(j),&aj) ;
	if (!r)
	{ delete[] _row_lhs ;
	  if (aj) pkvec_free(aj) ;
	  return (0) ; }
	for (int l = 0 ; l < aj->cnt ; l++)
	{ int i = inv(aj->coeffs[l].ndx) ;
	  _row_lhs[i] += x[j]*aj->coeffs[l].val ; } } }
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

  if (obj_sense < 0)
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
  assert(col_count > 0) ;
  int coeff_count = consys->mtx.coeffcnt ;
  assert(coeff_count > 0) ;

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

  int n = getNumRows() ;
  double* lower = new double[n] ;

  for (int i = 0 ; i < n ; i++)
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
      { lower[i] = -DYLP_INFINITY ;
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

  int n = getNumRows() ;
  double* upper = new double[n] ;

  for (int i = 0 ; i < n ; i++)
  { if (consys->ctyp[idx(i)] == contypGE || consys->ctyp[idx(i)] == contypNB)
      upper[i] = DYLP_INFINITY ;
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
  
  int n = getNumRows() ;
  char* sense = new char[n] ;
  const contyp_enum* ctyp = INV_VEC(contyp_enum,consys->ctyp) ;

  std::transform(ctyp,ctyp+n,sense,type_to_sense) ;

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
  
  int n = getNumRows() ;  
  double* range = new double[n] ;
  const double* lower = getRowLower() ;
  const double* upper = getRowUpper() ;
  const char* sense = getRowSense() ;

  for (int i=0 ; i<n ; i++)
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
  
  int n = getNumRows() ;  
  double* rhs = new double[n] ;
  const double* lower = getRowLower() ;
  const double* upper = getRowUpper() ;
  const char* sense = getRowSense() ;

  for (int i = 0 ; i < n ; i++)
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
  Effectively unsupportable in dylp, due to the use of the dynamic simplex
  algorithm. Alternating addition/deletion of constraints and variables
  means that the objective does not change monotonically.
*/

bool ODSI::isDualObjectiveLimitReached () const

{ throw CoinError("Unsupported: Objective does not change monotonically "
		  "in a dynamic simplex algorithm.",
		  "isDualObjectiveLimitReached",
		  "OsiDylpSolverInterface") ;

  return false ; }


/*!
  Effectively unsupportable. See
  \link OsiDylpSolverInterface::isDualObjectiveLimitReached
	ODSI:isDualObjectiveLimitReached\endlink.
*/

bool ODSI::isPrimalObjectiveLimitReached () const

{ throw CoinError("Unsupported: Objective does not change monotonically "
		  "in a dynamic simplex algorithm.",
		  "isPrimalObjectiveLimitReached",
		  "OsiDylpSolverInterface") ;

  return false ; }


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
  This can be done after any call to a solver routine.

  The \link OsiDylpSolverInterface::setWarmStart ODSI::setWarmStart \endlink
  routine rebuilds the required basis and status information from the
  OsiDylpWarmStartBasis object and sets ODSI options so that the solver will
  attempt a warm start at the next call to
  \link OsiDylpSolverInterface::resolve ODSI::resolve \endlink.

  The convention in OSI is that warm start objects derive from the base class
  CoinWarmStart. If you're using only getWarmStart and setWarmStart, that's
  really all you need to know.
  
  The dylp warm start object OsiDylpWarmStartBasis derives from
  CoinWarmStartBasis, which in turn derives from CoinWarmStart. A new derived
  class is needed because dylp does not always work with the full constraint
  system. This means the warm start object must specify the active
  constraints.  For ease of use, constraint status is handled just like
  variable status. There's an array, one entry per constraint, coded as
  CoinWarmStartBasis::Status.  Inactive constraints have status isFree, active
  constraints use atLowerBound.
*/
//@{

/*!
  This routine constructs an OsiDylpWarmStartBasis structure from the basis
  and status information returned as a matter of course by the dylp solver.
*/

CoinWarmStart* ODSI::getWarmStart () const

/*
  This routine constructs a OsiDylpWarmStartBasis structure from the basis
  returned by dylp.

  Dylp takes the atttitude that in so far as it is possible, logicals should
  be invisible outside the solver. The status of nonbasic logicals is not
  reported, nor does dylp expect to receive it.
*/

{ int i,j,k ;
  flags statj ;

  if (!lpprob) return (0) ;
/*
  Create and size an OsiDylpWarmStartBasis object. setSize initialises all
  status entries to isFree. Then grab pointers to the status vectors and the
  dylp basis vector.
*/
  OsiDylpWarmStartBasis *wsb = new OsiDylpWarmStartBasis ;
  wsb->setSize(consys->varcnt,consys->concnt) ;

  char *const strucStatus = wsb->getStructuralStatus() ;
  char *const artifStatus = wsb->getArtificialStatus() ;
  char *const conStatus = wsb->getConstraintStatus() ;
  basis_struct *basis = lpprob->basis ;

  if (lpprob->lpret == lpOPTIMAL)
    wsb->setPhase(dyPRIMAL2) ;
  else
    wsb->setPhase(dyPRIMAL1) ;
/*
  Walk the basis and mark the active constraints and basic variables. Basic
  logical variables are entered as the negative of the constraint index.
*/
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
  Now scan the status vector and record the status of nonbasic structural
  variables. Some information is lost here --- CoinWarmStartBasis::Status
  doesn't encode NBFX.
*/
  for (j = 1 ; j <= consys->varcnt ; j++)
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
	  wsb = 0 ; } } } }

  return (wsb) ; }


/*!
  This routine installs the basis snapshot from an OsiDylpWarmStartBasis
  structure and sets ODSI options so that dylp will attempt a warm start on
  the next call to
  \link OsiDylpSolverInterface::resolve ODSI::resolve \endlink.
*/

bool ODSI::setWarmStart (const CoinWarmStart *ws)

{ int i,j,k ;
  CoinWarmStartBasis::Status osi_stati ;

/*
  Use a dynamic cast to make sure we have an OsiDylpWarmStartBasis.
*/
  const OsiDylpWarmStartBasis *wsb =
			dynamic_cast<const OsiDylpWarmStartBasis *>(ws) ;
  if (!wsb) return (false) ;

/*
  Do we have an lpprob structure yet? If not, construct one.
*/
  if (!lpprob) construct_lpprob() ;

  assert(resolveOptions && consys && consys->vlb && consys->vub) ;

/*
  Extract the info in the warm start object --- size and status vectors.  The
  number of variables and constraints in the warm start object should match
  the full size of the constraint system.

  Create a dylp basis_struct and status vector of sufficient size to hold the
  information in the OsiDylpWarmStartBasis. This space may well be freed or
  realloc'd by dylp, so use standard calloc to acquire it.  We'll only use as
  much of the basis as is needed for the active constraints.
*/
  int varcnt = wsb->getNumStructural() ;
  int concnt = wsb->getNumArtificial() ;
  assert(varcnt == getNumCols() && concnt == getNumRows()) ;

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
  { osi_stati = getStatus(conStatus,inv(i)) ;
    if (osi_stati == CoinWarmStartBasis::atLowerBound)
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
  { osi_stati = getStatus(strucStatus,inv(j)) ;
    switch (osi_stati)
    { case CoinWarmStartBasis::basic:
      { k++ ;
	assert(k <= actcons) ;
	basis.el[k].vndx = j ;
	status[j] = (flags) (-k) ;
	break ; }
      case CoinWarmStartBasis::atLowerBound:
      { if (consys->vlb[j] == consys->vub[j])
	{ status[j] = vstatNBFX ; }
	else
	{ status[j] = vstatNBLB ; }
	break ; }
      case CoinWarmStartBasis::atUpperBound:
      { status[j] = vstatNBUB ;
	break ; }
      case CoinWarmStartBasis::isFree:
      { status[j] = vstatNBFR ;
	break ; } } }
/*
  Now we need to finish out the basis by adding the basic logicals.  The
  convention to represent this with the negative of the index of the
  constraint which spawns the logical.
*/
  for (i = 1 ; i <= concnt ; i++)
  { osi_stati = getStatus(artifStatus,inv(i)) ;
    if (osi_stati == CoinWarmStartBasis::basic)
    { k++ ;
      assert(k <= actcons) ;
      basis.el[k].vndx = -i ; } }
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
  
  return (true) ; }


void ODSI::resolve ()
/*
  This routine simply calls the solver, after making sure that the forcecold
  option is turned off.  We do need to make sure the solver is ours to use.

  If we're reoptimising, then the basis should be ready and we should have
  warm start information. But the client could be doing something clever, like
  using a warm start object to jump start a problem in a virgin solver.
*/
{ assert(lpprob && lpprob->basis && lpprob->status && consys &&
	 resolveOptions && tolerances && statistics) ;

  if (dylp_owner != 0 && dylp_owner != this) dylp_owner->detach_dylp() ;

/*
  Is the basis package initialised?
*/
  if (basis_ready == false)
  { int count = static_cast<int>((1.5*getNumRows()) + 2*getNumCols()) ;
    dy_initbasis(count,initialSolveOptions->factor,0) ;
    basis_ready = true ; }

  if (isactive(local_logchn)) logchn = local_logchn ;
  gtxecho = resolve_gtxecho ;

  dyphase_enum phase = lpprob->phase ;
  if (!(phase == dyPRIMAL1 || phase == dyPRIMAL2 || phase == dyDUAL))
    lpprob->phase = dyINV ;

  resolveOptions->forcecold = false ;

  lp_retval = dylp_dolp(lpprob,resolveOptions,tolerances,statistics) ;  
  if (flgon(lpprob->ctlopts,lpctlDYVALID))
  { dylp_owner = this ;
    _objval = obj_sense*lpprob->obj ; }

  destruct_col_cache() ;
  destruct_row_cache() ; }

//@} // WarmStart



/*! \defgroup HotStartMethods Hot Start Methods
    \brief Methods for solving an LP from a hot start.

  In the absence of instructions to the contrary, dylp will always go for a
  hot start. The obligation of the caller is simply to insure that dylp's
  internal static data structures are intact from the previous call to the
  solver. In the context of OSI, this means no intervening calls to the
  solver by any ODSI object.

  The restrictions outlined in the OSI documentation are somewhat relaxed in
  dylp: you can change upper or lower bounds on variables, upper or lower
  bounds on constraints, or objective function coefficients. No changes to
  the coefficient matrix are allowed --- for this, you need to drop back to
  a warm start.

  There is, however, no hot start object, and the semantics of a hot start
  for dylp are probably not quite as you would expect (or as implied by the
  OSI specification). dylp simply assumes that its data structures are intact
  and resumes pivoting after making any necessary adjustments for allowable
  changes.  If you want the ability to repeatedly reset the solver to a
  particular basis, you need to use a warm start.

  Philosophically, a hot start tries to avoid having to refactor the basis
  as part of the restart. Resetting to a particular basis or changing the
  coefficient matrix would force a refactor before pivoting could resume.
*/
//@{

/*!
  This routine sets options in the ODSI object so that dylp will attempt a
  hot start at each call of
  \link OsiDylpSolverInterface::solveFromHotStart
	ODSI::solveFromHotStart \endlink.
*/

inline void ODSI::markHotStart ()
/*
  We should check that dylp has valid data structures, but there's no
  provision to return failure ... Otherwise, all we do is turn off the
  forcecold and forcewarm options.
*/

{ assert(lpprob && resolveOptions) ;

  assert(flgon(lpprob->ctlopts,lpctlDYVALID)) ;

  resolveOptions->forcecold = false ;
  resolveOptions->forcewarm = false ;

  return ; }

/*!
  This routine simply calls the solver. For dylp's model of hot start, that's
  all that's required.
*/

void ODSI::solveFromHotStart ()
/*
  This routine simply calls the solver, after making sure that the forcecold
  and forcewarm options are turned off. If we're going for a hot start, the
  basis should be ready.

  Note that dylp does need to know what's changed: any of bounds, rhs &
  rhslow, or objective. The various routines that make these changes should
  set the flags.
*/

{ int tmp_iterlim = -1 ;
  int hotlim ;

  assert(lpprob && lpprob->basis && lpprob->status && basis_ready &&
	 consys && resolveOptions && tolerances && statistics) ;

  if (isactive(local_logchn)) logchn = local_logchn ;
  gtxecho = resolve_gtxecho ;
/*
  This is overly cautious, but the only practical solution. dylp actually
  expects the user to set the proper phase after tweaking the problem. That's
  not practical in ODSI --- exposing this to the client would violate the
  generic nature of the interface.
*/
  lpprob->phase = dyPRIMAL1 ;
  resolveOptions->forcecold = false ;
  resolveOptions->forcewarm = false ;
  getIntParam(OsiMaxNumIterationHotStart,hotlim) ;
  if (hotlim > 0)
  { tmp_iterlim = resolveOptions->iterlim ;
    resolveOptions->iterlim = (hotlim/3 > 0)?hotlim/3:1 ; }

  lp_retval = dylp_dolp(lpprob,resolveOptions,tolerances,statistics) ;
  _objval = obj_sense*lpprob->obj ;

  if (tmp_iterlim > 0) resolveOptions->iterlim = tmp_iterlim ;

  destruct_col_cache() ;
  destruct_row_cache() ; }


inline void ODSI::unmarkHotStart ()
/*
  I can't think of anything to do here. There's no point in forcing one of
  the other start modes.
*/

{ return ; }

//@} // HotStartMethods



/*! \defgroup DylpMethods Dylp-Specific Methods */
//@{

/*!
  Process a dylp options (.spc) file.

  \param name file name
  \param silent true to suppress command echo
  \param mustexist true (default) if the file must exist
*/

void ODSI::dylp_controlfile (const char *name,
			     const bool silent, const bool mustexist)

{ if (name == 0 || *name == 0) return ;
  string mode = (mustexist)?"r":"q" ;
  cmdchn = openfile(const_cast<char *>(name),const_cast<char *>(mode.c_str())) ;
  if (!(cmdchn == IOID_INV || cmdchn == IOID_NOSTRM))
  { setmode (cmdchn, 'l') ;  
    main_lpopts = initialSolveOptions ;
    main_lptols = tolerances ;
    bool r = (process_cmds(silent) != 0) ;
    (void) closefile(cmdchn) ;
    assert(r == cmdOK) ; }
/*
  Differentiate options for initial solution vs. reoptimising.
*/
  memcpy(resolveOptions,initialSolveOptions,sizeof(lpopts_struct));
  initialSolveOptions->forcecold = true ;
  initialSolveOptions->fullsys = true ;

  cmdchn = IOID_NOSTRM ; }
  

/*!
  Establish a log (.log) file.

  \param name file name
  \param silent true to echo log messages to the terminal
*/

void ODSI::dylp_logfile (const char *name, bool echo)

{ if (name == 0 || *name == 0) return ;

  string log = make_filename(name,".mps",".log") ;
  local_logchn = openfile(const_cast<char *>(log.c_str()),
			  const_cast<char *>("w")) ;
  if (local_logchn == IOID_INV) local_logchn = IOID_NOSTRM ;
  initial_gtxecho = echo ;
  resolve_gtxecho = echo ; }

//@} // DylpMethods

#endif /* COIN_USE_DYLP */
