/*! \legal
  Copyright (C) 2002, Lou Hafer, Stephen Tse, International Business Machines
  Corporation and others. All Rights Reserved.
*/

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
  around this, there are a number of macros defined below: INV_VEC, CLONE,
  CLONE_VEC, and COPY_VEC. For Sun Workshop and Gnu programming environments,
  these will expand to template functions. For MS, they expand to some code
  plus calls to STL functions. Apparently MS C++ doesn't have trouble with the
  definitions of the templates, just instantiations.
*/

#endif

/* Cut name length for readability. */

#define ODSI OsiDylpSolverInterface


#ifndef ERRMSGDIR
# define ERRMSGDIR "."
#endif

/*!
   \file OsiDylpSolverInterface.cpp

   \brief Implementation of COIN OSI layer for dylp.

   This file contains the implementation of a COIN OSI layer for dylp, the lp
   solver for the bonsaiG MILP code.

   More information on the COIN/OR project and the OSI layer specification
   can be found at http://www.coin-or.org.

   -- Implementation Details --

   CACHE MECHANISM: Since vectors returned by OsiSolverInterface "get"
   functions are constant (not modifiable by users) and it is convenient for
   users to make repeated calls for the same information, a cache mechanism
   is used to avoid repeated expensive conversions in the interface. All
   cache variables are prefixed with underscore (_col*, _row*, _matrix*).
   Note, however, that we need to invalidate or update the cache whenever the
   problem is modified.

   INDEX BASE: OsiSolverInterface indexes variables and constraints from 0
   (the standard approach for C/C++) while dylp indexes them from 1 (for
   robustness; see consys.h). Some caution is needed in the interface to
   avoid off-by-one errors. Similarly, arrays (vectors) for k elements in
   bonsai occupy k+1 space with vector[0] unused.  ODSI::idx, ODSI::inv, and
   ODSI::inv_vec are used to make these trivial conversions explicit.

   COPY AND ASSERT: Since OsiSolverInterface requires cloning while dylp does
   not provide that internally, a good part of this code is devoted to doing
   the job. To clone, each part of dylp's data structures (lpprob, consys,
   basis) is copied (with ODSI::copy*) and then (optionally) verified (with
   ODSI::assert_same* and assert). Helpers for copying and verifying
   primitive types (integer and double) and array types are implemented with
   C++ templates.
*/

static char sccsid[] = "@(#)OsiDylpSolverInterface.cpp	1.9	11/26/02" ;

#include <string>
#include <cassert>
#include <OsiColCut.hpp>
#include <OsiRowCut.hpp>
#include "OsiDylpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinMpsIO.hpp"

using std::string ;
using std::vector ;


extern "C"
{

#include "bonsai.h"

extern bool dy_mpsin(const char* filename, consys_struct** consys) ;
extern void dy_initbasis(int concnt, int factor_freq, double zero_tol) ;
extern void dy_freebasis() ;
extern char *stralloc(char *str) ;

}

#define UNSUPPORTED assert(0)

/*
  Dylp uses HUGE_VAL for infinity; generally this seems to resolve to IEEE
  infinity.
*/

#define DYLP_INFINITY HUGE_VAL

/*!
  \defgroup DylpIO dylp i/o control variables
  \brief Variables controlling dylp file and terminal i/o.

  Dylp is capable of generating a great deal of output, but the control
  mechanisms are not a good fit for C++ objects and the OSI framework.
  Four global variables, cmdchn, cmdecho, logchn, and gtxecho, control
  command input/echoing and log message output/echoing.

  By default, solver command (option) files are handled in the following
  manner: When an MPS file example.mps is read, the solver looks for an
  options file example.spc in the same directory. If it exists, it is
  processed. This is a local file open/close, well-defined and clean, but
  with little user control. Commands are not echoed during processing.

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
#ifndef NDEBUG
bool cmdecho = false,
     gtxecho = false ;
#else
bool cmdecho = false,
     gtxecho = false ;
#endif

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
  an lp (initialSolve, basis_ready == false) and shut down when the last ODSI
  object is destroyed (destructor, reference_count == 0, basis_ready == true).
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
  inv, and inv_vec make this conversion explicit.
  See the \link OsiDylpSolverInterface.cpp file comments \endlink for a brief
  description.
  This group also provides a function for conversion between the OSI notion of
  a packed vector and the dylp notion of a packed vector.
*/
//@{

/*! \brief Convert 0-based index to 1-based index */

inline int ODSI::idx (int i) { return i + 1 ; }

/*! \brief Convert 1-based index to 0-based index */

inline int ODSI::inv (int i) { return i - 1 ; }


/*! \brief Convert 1-based vector pointer to 0-based vector pointer

  For cases where it's inconvenient to adjust indices, the alternative is to
  adjust the pointer to the vector so it points to vector[1].
*/
/*
  This function used to be called just inv, but ... Sun CC can't figure out
  the right type given inv(some_vector), and Gnu g++ claims inv<some_type> is
  a syntax error! Both are satisfied if the function is named inv_vec.
*/

template<class T> inline T* ODSI::inv_vec (T* vec) { return vec + 1 ; }

#ifdef _MSC_VER
#  define INV_VEC(zz_type_zz,zz_vec_zz) (zz_vec_zz+1)
#else
#  define INV_VEC(zz_type_zz,zz_vec_zz) inv_vec<zz_type_zz>(zz_vec_zz)
#endif

/*! \brief Convert an CoinShallowPackedVector to a dylp pkvec.

  pkvec's created with this routine must eventually be freed with pkvec_free.
  Note that 0-based indexing is used for the entries of the coeff vector
  (i.e., coeffs[0] is valid).

  \param dimension The length of the vector if it were to be unpacked.
*/

pkvec_struct* ODSI::packed_vector (const CoinShallowPackedVector src,
				   int dimension)

{ int n = src.getNumElements() ;
  const int* indices = src.getIndices() ;
  const double* elements = src.getElements() ;
  assert(indices && elements) ;
 
  pkvec_struct* dst = pkvec_new(n) ;

  dst->cnt = n ;
  dst->dim = dimension ;
  pkcoeff_struct* coeffs = dst->coeffs ;

  for (int i = 0 ; i < n ; i++)
  { coeffs[i].ndx = idx(indices[i]) ;
    coeffs[i].val = elements[i] ; }

  return dst ; }


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
	    { /*zz_dst_zz = (zz_type_zz *)malloc(sizeof(zz_type_zz)) ;*/ \
	      zz_dst_zz = new zz_type_zz ; \
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

void ODSI::copy_consysmtx (const consys_struct* src, consys_struct* dst)

{
/*
  Check that dst has a properly formed but empty constraint matrix.
*/
  assert(dst && dst->mtx.cols && dst->mtx.rows && dst->mtx.coeffcnt == 0) ;
/*
  Add an empty row to dst for each row in src.
*/
  pkvec_struct* drow = pkvec_new(0) ;
  int i ;

  for (i = idx(0) ; i < idx(src->concnt) ; i++)
  { rowhdr_struct* srow = src->mtx.rows[i] ;
    drow->ndx = srow->ndx ;
    drow->nme = srow->nme ;
    contyp_enum type = src->ctyp[i] ;
    bool r = consys_addrow_pk(dst, (i <= src->archccnt)?'a':'c',
			      type, drow, 0.0, 0.0, 0, 0) ;
    assert(r) ; }
  if (drow) pkvec_free(drow) ;
/*
  Add the coefficients, column by column.
*/
  pkvec_struct* col = 0 ;
  for (i = idx(0) ; i < idx(src->varcnt) ; i++)
  { consys_struct* src2 = const_cast<consys_struct*>(src) ;
    bool r1 = consys_getcol_pk(src2, i, &col) ;
    assert(r1) ;
    vartyp_enum type = src->vtyp[i] ;
    bool r2 = consys_addcol_pk(dst, type, col, 0.0, 0.0, 0.0) ;
    assert(r2) ; }
  if (col) pkvec_free(col) ;

  return ; }


/*! \brief Specialised copy function for a consys_struct

  This routine calls consys_create to create and initialise an empty
  consys_struct, then calls copy_consysmtx to copy the coefficient matrix.
  Finally, any associated vectors (objective, bounds, right-hand-side) are
  copied.

  \note Any additional vectors that have been attached to the constraint
	system are not copied. Arguably, they are not really part of the
	constraint system.
*/

consys_struct* ODSI::copy_consys (const consys_struct* src)

{ if (!src) return 0 ;

  flags parts = src->parts ;
  flags options = src->opts ;
/*
  Call consys_utils:consys_create to create a new header structure and attach
  (empty) vectors as specified in src.parts. ODSI::copy_consysmtx clones the
  coefficient matrix.
*/
  consys_struct* dst = consys_create(src->nme, parts, options, 
				     src->colsze, src->rowsze) ;
  assert(dst) ;
  assert(dst->nme == src->nme) ;
  assert(dst->parts == src->parts) ;
  assert(dst->opts == src->opts) ;
  assert(dst->colsze == src->colsze) ;
  assert(dst->rowsze == src->rowsze) ;

/*
  Copy the constraint matrix, then verify (sort of) by checking the various
  counts that indicate the contents of the coefficient matrix along with the
  variable and constraint type vectors. Note that we check the allocated size
  a second time. It should still be the same.
*/
  int var_count = idx(src->varcnt) ;
  int con_count = idx(src->concnt) ;
  copy_consysmtx(src, dst) ;

  assert(dst->varcnt == src->varcnt) ;
  assert(dst->archvcnt == src->archvcnt) ;
  assert(dst->logvcnt == src->logvcnt) ;
  assert(dst->intvcnt == src->intvcnt) ;
  assert(dst->binvcnt == src->binvcnt) ;
  assert(dst->maxcollen == src->maxcollen) ;
  assert(dst->maxcolndx == src->maxcolndx) ;
  assert(dst->concnt == src->concnt) ;
  assert(dst->archccnt == src->archccnt) ;
  assert(dst->cutccnt == src->cutccnt) ;
  assert(dst->maxrowlen == src->maxrowlen) ;
  assert(dst->maxrowndx == src->maxrowndx) ;
  assert(dst->colsze == src->colsze) ;
  assert(dst->rowsze == src->rowsze) ;
# ifndef NDEBUG
  assert_same(dst->vtyp, src->vtyp,var_count,true) ;
  assert_same(dst->ctyp, src->ctyp,con_count,true) ;
# endif

/*
  Copy fields related to the objective function.
*/
  dst->objnme = stralloc(src->objnme) ;
  dst->objndx = src->objndx ;
  dst->xzndx = src->xzndx ;
/*
  Copy the contents of attached vectors. In theory this could have been
  done in copy_consysmtx (as parameters to addrow_pk and addcol_pk), but a
  block copy should be more efficient here.
*/
  bool do_obj = (parts & CONSYS_OBJ) != 0 ;
  bool do_vub = (parts & CONSYS_VUB) != 0 ;
  bool do_vlb = (parts & CONSYS_VLB) != 0 ;
  bool do_rhs = (parts & CONSYS_RHS) != 0 ;
  bool do_rhslow = (parts & CONSYS_RHSLOW) != 0 ;
  bool do_cub = (parts & CONSYS_CUB)!= 0 ;
  bool do_clb = (parts & CONSYS_CLB) != 0 ;

  if (do_obj) COPY_VEC(double,src->obj,dst->obj,var_count) ;
  if (do_vub) COPY_VEC(double,src->vub,dst->vub,var_count) ;
  if (do_vlb) COPY_VEC(double,src->vlb,dst->vlb,var_count) ;
  if (do_rhs) COPY_VEC(double,src->rhs,dst->rhs,con_count) ;
  if (do_rhslow) COPY_VEC(double,src->rhslow,dst->rhslow,con_count) ;
  if (do_cub) COPY_VEC(conbnd_struct,src->cub,dst->cub,con_count) ;
  if (do_clb) COPY_VEC(conbnd_struct,src->clb,dst->clb,con_count) ;

# ifndef NDEBUG
  assert_same(*dst, *src, true) ;
# endif

  return dst ;
}


/*! \brief Specialised copy function for a dylp lpprob_struct */

lpprob_struct* ODSI::copy_lpprob (const lpprob_struct* src)

{ if (!src) return 0 ;

  int col_count = idx(src->colsze) ;
  int row_count = idx(src->rowsze) ;

  lpprob_struct* dst = NULL;
  CLONE(lpprob_struct,src,dst);

  dst->consys = copy_consys(src->consys) ;
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

    This group of functions implement the conversions between the OSI notion of
    a constraint system and the dylp notion of a constraint system.
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
*/

void ODSI::add_col (const CoinPackedVectorBase& osi_coli,
		    vartyp_enum vtypi, double vlbi, double vubi, double obji)

{ pkvec_struct *pk_coli = packed_vector(osi_coli,getNumRows()) ;

  bool r = consys_addcol_pk(consys,vtypi,pk_coli,obji,vlbi,vubi) ;
  if (pk_coli) pkvec_free(pk_coli) ;
  assert(r) ;
  if (options) options->forcewarm = true ;

  destruct_cache() ; }


/*! \brief Add a row to the constraint matrix.

  Convert the OSI packed vector to a dylp packed vector and install the row.
  The remaining parameters are expected to be in the correct dylp form.
*/

void ODSI::add_row (const CoinPackedVectorBase& osi_rowi, char clazzi,
		    contyp_enum ctypi, double rhsi, double rhslowi)

{ pkvec_struct *pk_rowi = packed_vector(osi_rowi,getNumCols()) ;

  bool r = consys_addrow_pk(consys,clazzi,ctypi,pk_rowi,rhsi,rhslowi,0,0) ;
  if (pk_rowi) pkvec_free(pk_rowi) ;
  assert(r) ;
  if (options) options->forcewarm = true ;

  destruct_cache() ; }


//@} // ConsysHelpers



/* OsiDylpSolverInterface constructors and related subroutines. */


/*! \defgroup ConstructorHelpers Helper functions for problem construction */
//@{

/*! \brief Construct a dylp lpprob_struct (LP problem).

    \todo Remove the call to dy_setprintopts, once I've figured out a better
	  scheme for handling dylp parameters.
*/

void ODSI::construct_lpprob ()

{
  dy_checkdefaults(consys, options,tolerances) ;
  lpprob = new lpprob_struct ;
  memset(lpprob, 0, sizeof(lpprob_struct)) ;
  setflg(lpprob->ctlopts,lpctlNOFREE) ;
  lpprob->phase = dyINV ;
  lpprob->consys = consys ;
  lpprob->rowsze = consys->rowsze ;
  lpprob->colsze = consys->colsze ;
  dy_setprintopts(1, options) ;
}


/*! \brief Construct and initialize default options, tolerances, and
	   statistics structures

  If any OSI parameters have been set, they will be incorporated here.
*/

void ODSI::construct_options ()

{
/*
  Acquire the default options and tolerances, then modify them to reflect any
  OSI parameters.
*/
  dy_defaults(&options,&tolerances) ;

  if (osi_dual_tolerance > 0.0)
  { tolerances->dfeas_scale = osi_dual_tolerance/tolerances->cost ; }
  if (osi_primal_tolerance > 0.0)
  { tolerances->pfeas_scale = osi_primal_tolerance/tolerances->zero ; }
  if (osi_iterlim > -1)
  { options->iterlim = osi_iterlim/3 ; }
/*
  Acquire and clear a statistics structure.
*/
  statistics = new lpstats_struct ;
  memset(statistics,0,sizeof(lpstats_struct)) ;
}


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
  packed matrix and dylp constraint sense and right-hand-side (rhs, rhslow)
  vectors.

  The first action is to free the existing problem.
  Existing options and tolerances settings are saved, if they exist, and will
  be restored after the new problem is loaded.

  Next the routine builds an empty consys_struct of the proper size with a
  call to consys_create. The main body of the routine fills the constraint
  matrix in two loops. In the first, calls to consys_addrow_pk insert empty
  constraints and establish the constraint type and rhs. Then calls to
  consys_addcol_pk insert the variables, along with their bounds and
  objective coefficient.

  Finally, an lpprob_struct is created to hold the new problem, and options
  and tolerances are reattached or initialised from defaults.
*/

void ODSI::load_problem (const CoinPackedMatrix& matrix,
	   const double* col_lower, const double* col_upper, const double* obj,
	   const contyp_enum *ctyp, const double* rhs, const double* rhslow)

{ lpopts_struct *preserve_options = 0 ;
  lptols_struct *preserve_tolerances = 0 ;

/*
  Preserve existing options and tolerances, then free the existing problem
  structures.
*/
  if (options)
  { preserve_options = options ;
    options = 0 ; }
  if (tolerances)
  { preserve_tolerances = tolerances ;
    tolerances = 0 ; }
  destruct_problem() ;
/*
  Create an empty consys_struct. Request the standard attached vectors:
  objective, variable upper & lower bounds, variable and constraint types,
  and constraint right-hand-side and range (rhslow). On the off chance that
  we might duplicate a matrix with 0 rows or columns, quietly bump the value
  passed to consys_create to 1, so that the user doesn't have to be aware that
  consys_create interprets 0 as `use the (rather large) default value'.
*/
  int colcnt = matrix.getNumCols() ;
  int rowcnt = matrix.getNumRows() ;
  flags consys_parts = CONSYS_OBJ | CONSYS_VUB | CONSYS_VLB |
		CONSYS_RHS | CONSYS_RHSLOW | CONSYS_VTYP | CONSYS_CTYP ;
  flags consys_options = CONSYS_WRNATT ;

  consys = consys_create(0,consys_parts,consys_options,
			 (rowcnt <= 0)?1:rowcnt,(colcnt <= 0)?1:colcnt) ;
  assert(consys) ;
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
  Second loop. Insert the coefficients by column.
*/
  CoinPackedMatrix matrix2 = matrix ;
  if (!matrix2.isColOrdered()) matrix2.reverseOrdering() ;

  for (int j = 0 ; j < colcnt ; j++)
  { const CoinShallowPackedVector osi_col = matrix2.getVector(j) ;
    pkvec_struct* colj = packed_vector(osi_col,rowcnt) ;
    colj->nme = 0 ;
    double objj = obj?obj[j]:0 ;
    double vlbj = col_lower?col_lower[j]:0 ;
    double vubj = col_upper?col_upper[j]:DYLP_INFINITY ;
    bool r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    pkvec_free(colj) ;
    assert(r) ;
  }
  assert(matrix2.isEquivalent(*getMatrixByCol())) ;
/*
  Get a default set of options and tolerances, then replace them with the
  preserved parameters, if any. Finally, construct an lpprob_struct.
*/
  construct_options() ;
  if (preserve_options)
  { /* delete options ;removed cause sometimes malloc'ed and should then be freed to avoid crashes */
    options = preserve_options ; }
  if (preserve_tolerances)
  { /* delete tolerances ; removed cause sometimes malloc'ed and should then be freed to avoid crashes */
    tolerances = preserve_tolerances ; }
  construct_lpprob() ;

}

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
  call to consys_create. The main body of the routine fills the constraint
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

{ lpopts_struct *preserve_options = 0 ;
  lptols_struct *preserve_tolerances = 0 ;

/*
  Preserve existing options and tolerances, then free the existing problem
  structures.
*/
  if (options)
  { preserve_options = options ;
    options = 0 ; }
  if (tolerances)
  { preserve_tolerances = tolerances ;
    tolerances = 0 ; }
  destruct_problem() ;
/*
  Create an empty consys_struct. Request the standard attached vectors:
  objective, variable upper & lower bounds, variable and constraint types,
  and constraint right-hand-side and range (rhslow). The user may ask for 0
  rows or columns (intending to add them later). In this case, we quietly
  bump the value passed to consys_create to 1, so that the user doesn't have
  to be aware that consys_create interprets 0 as `use the (rather large)
  default value'.
*/
  flags consys_parts = CONSYS_OBJ | CONSYS_VUB | CONSYS_VLB |
		CONSYS_RHS | CONSYS_RHSLOW | CONSYS_VTYP | CONSYS_CTYP ;
  flags consys_options = CONSYS_WRNATT ;

  consys = consys_create(0,consys_parts,consys_options,
			 (rowcnt <= 0)?1:rowcnt,(colcnt <= 0)?1:colcnt) ;
  assert(consys) ;
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
  overestimate, but it saves us the trouble of constantly reallocating it.
*/
  pkvec_struct *colj = pkvec_new(rowcnt) ;
  assert(colj) ;
  colj->dim = rowcnt ;
  colj->nme = 0 ;
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
    bool r = consys_addcol_pk(consys,vartypCON,colj,objj,vlbj,vubj) ;
    assert(r) ;
  }

  if (colj) pkvec_free(colj) ;
/*
  Get a default set of options and tolerances, then replace them with the
  preserved parameters, if any. Finally, construct an lpprob_struct.
*/
  construct_options() ;
  if (preserve_options)
  { delete options ;
    options = preserve_options ; }
  if (preserve_tolerances)
  { delete tolerances ;
    tolerances = preserve_tolerances ; }
  construct_lpprob() ;

}



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
  Creates the shell of an ODSI object. The options and tolerances substructures
  are instantiated, but there are no lp problem or constraint system
  structures.

  Creation of the first ODSI object triggers initialisation of the i/o
  subsystem.
*/

ODSI::OsiDylpSolverInterface ()

  : options(0),
    tolerances(0),
    consys(0),
    lpprob(0),

    local_logchn(IOID_NOSTRM),
    local_gtxecho(false),
    lp_retval(lpINV),
    obj_sense(1.0),

    osi_dual_tolerance(0.0),
    osi_primal_tolerance(0.0),
    osi_obj_offset(0.0),
    osi_iterlim(-1),
    osi_hotiterlim(-1),
    osi_probname(),

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
  Establish structures for parameters, tolerances, and statistics.
*/
  ODSI::construct_options() ;
  reference_count++ ;

  if (reference_count == 1)
  { dylp_ioinit() ;
    CoinRelFltEq eq ;
    assert(eq(DYLP_INFINITY, DYLP_INFINITY)) ; }

}


/*!
  A true copy --- no data structures are shared with the original. Cached
  information is also replicated.
*/

ODSI::OsiDylpSolverInterface (const OsiDylpSolverInterface& src)

  : local_gtxecho(src.local_gtxecho),
    lp_retval(src.lp_retval),
    obj_sense(src.obj_sense),

    osi_dual_tolerance(src.osi_dual_tolerance),
    osi_primal_tolerance(src.osi_primal_tolerance),
    osi_obj_offset(src.osi_obj_offset),
    osi_iterlim(src.osi_iterlim),
    osi_hotiterlim(src.osi_hotiterlim),
    osi_probname(src.osi_probname),
  
    _matrix_by_row(0),
    _matrix_by_col(0)

{ lpprob = copy_lpprob(src.lpprob) ;
  consys = lpprob ? lpprob->consys : 0 ;

  CLONE(lpopts_struct,src.options,options) ;
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

  reference_count++ ;

# ifndef NDEBUG
  assert_same(*this, src, true) ;
# endif

}


/*!
  An alternate copy constructor; the parameter is ignored (and assumed true).
*/

inline OsiSolverInterface* ODSI::clone (bool copyData) const

{
  return new OsiDylpSolverInterface(*this) ;
}


//@}



/* OsiDylpSolverInterface destructor and related subroutines. */

/*! \defgroup DestructorHelpers Destructor Helpers */
//@{

/*! \brief Free a dylp lpprob_struct

  Free the problem-specific structures in an ODSI object.  First the dylp
  lpprob_struct, the main structure passed to dylp. consys_free will free the
  constraint system, and dy_freesoln will free the remaining data structures
  associated with lpprob. There remains only to destroy lpprob itself. To
  finish the job, free the options, tolerances, and statistics structures and
  call destruct_cache to free the cache structures.
*/

void ODSI::destruct_problem ()

{ if (lpprob)
  { assert(lpprob->consys == consys) ;
    consys_free(consys) ;
    consys = 0 ;
    dy_freesoln(lpprob) ;
    delete lpprob ;
    lpprob = 0 ; }
  
  if (options)
  { /*delete options ;removed cause sometimes malloc'ed and should then be freed to avoid crashes */
    options = 0 ; }
  if (tolerances)
  { /* delete tolerances ;removed cause sometimes malloc'ed and should then be freed to avoid crashes */
    tolerances = 0 ; }
  if (statistics)
  { delete statistics ;
    statistics = 0 ; }

  destruct_cache() ;
}


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
  dylp(lpprob,options,tolerances,statistics) ;
  dylp_owner = 0 ; }

//@}

/*! \ingroup ODSIConstructorsDestructors
    \brief Destructor for OsiDylpSolverInterface

    Destruction of the last ODSI object triggers shutdown of the basis
    maintenance and i/o subsystems and deallocation of all dylp internal
    data structures.  The call to detach_dylp is solely to free data
    structures.

*/

ODSI::~OsiDylpSolverInterface ()

{ if (dylp_owner == this) detach_dylp() ;

  destruct_problem() ;
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

  consys->vtyp[idx(j)] = vartypINT ; }


inline void ODSI::setInteger (const int* indices, int len)

{ if (!consys->vtyp)
  { bool r = consys_attach(consys,CONSYS_VTYP,sizeof(vartyp_enum),
			   reinterpret_cast<void **>(&consys->vtyp)) ;
    assert(r) ; }

  for (int i = 0 ; i < len ; i++)
    consys->vtyp[idx(indices[i])] = vartypCON ; }


inline void ODSI::setColLower (int i, double val)

{ if (!consys->vlb)
  { bool r = consys_attach(consys,CONSYS_VLB,sizeof(double),
			   reinterpret_cast<void **>(&consys->vlb)) ;
    assert(r) ; }

  consys->vlb[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlLBNDCHG) ; }


inline void ODSI::setColUpper (int i, double val)

{ if (!consys->vub)
  { bool r = consys_attach(consys,CONSYS_VUB,sizeof(double),
			   reinterpret_cast<void **>(&consys->vub)) ;
    assert(r) ; }

  consys->vub[idx(i)] = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlUBNDCHG) ; }


/*!
  A call to this routine destroys all cached row values.
*/

void ODSI::setRowType (int i, char sense, double rhs, double range)

{ int k = idx(i) ;

  gen_rowiparms(&consys->ctyp[k],&consys->rhs[k],&consys->rhslow[k],
		sense,rhs,range) ;
  if (options) options->forcewarm = true ;

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

inline void ODSI::addRow (const CoinPackedVectorBase& row,
			  double lower, double upper)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,lower,upper) ;
  add_row(row,'a',ctypi,rhsi,rhslowi) ;

  return ; }


/*!
  Add a row to the constraint system given the coefficients, constraint sense,
  right-hand-side value, and range.

  A call to this routine destroys all cached values.
*/

inline void ODSI::addRow (const CoinPackedVectorBase& osi_rowi, 
			  char sense, double rhs, double range)

{ contyp_enum ctypi ;
  double rhsi,rhslowi ;

  gen_rowiparms(&ctypi,&rhsi,&rhslowi,sense,rhs,range) ;
  add_row(osi_rowi,'a',ctypi,rhsi,rhslowi) ;

  return ; }


/*!
  One cut at a time, please, when it comes to rows.

  A call to this routine destroys all cached values.
*/

void ODSI::applyRowCut (const OsiRowCut& cut)

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
  if (options) options->forcecold = true ;

  destruct_cache() ; }



/*!
  Create an objective function vector, if it doesn't already exist, then
  insert the new value. Update the cached copy, if it exists.
*/

void ODSI::setObjCoeff (int j, double objj)

{ if (!consys->obj)
  { bool r = consys_attach(consys,CONSYS_OBJ,sizeof(double),
			   reinterpret_cast<void **>(&consys->obj)) ;
    assert(r) ; }
  
  consys->obj[idx(j)] = obj_sense*objj ;
  if (_col_obj) _col_obj[j] = objj ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ; }
  

/*!
  If the objective exists and the desired sense disagrees with the current
  sense, multiply the objective by -1. Note that the cached copy is not
  affected.
*/

void ODSI::setObjSense (double val)

{ if (consys->obj && val != obj_sense)
  { for (int j = 1 ; j <= consys->varcnt ; j++) consys->obj[j] *= -1.0 ; }
  
  obj_sense = val ;
  if (lpprob) setflg(lpprob->ctlopts,lpctlOBJCHG) ; }


/*!
  Add a column to the constraint system given the coefficients, upper and lower
  bounds, and objective function coefficient. The variable type defaults to
  continuous.

  A call to this routine destroys all cached values.
*/

inline void ODSI::addCol (const CoinPackedVectorBase& osi_coli,
			  double vlbi, double vubi, double obji)

{ add_col(osi_coli,vartypCON,vlbi,vubi,obji) ; }


/*!
  Perhaps better named `applyColCuts', as an OsiColCut object can contain
  any number of adjustments to variable bounds. The routine simply steps
  through the new bounds, tightening existing bounds as necessary. Note that
  a bound is never loosened by this routine.
*/

void ODSI::applyColCut (const OsiColCut& cut)

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
  if (options) options->forcecold = true ;

  destruct_cache() ; }


//@} // ProbAdjust



#ifndef _MSC_VER

/*! \defgroup CopyVerifiers Copy Verification Routines
    \brief Private helper functions to verify correctness of copies

  This group of helper functions verify the correctness of copies of the
  data structures used in ODSI objects. For uniformity, all functions take
  a parameter, `exact', which influences the type of checks for some of the
  more complex structures.
*/
//@{


/*! \brief Compare two double values for equality.

 Exact requires exact bit-for-bit equality, which is usually not what you want.
 Inexact tests for equality within a tolerances of 1.0e-10, scaled by the
 maximum of the two doubles.
*/
void ODSI::assert_same (double d1, double d2, bool exact)

{ if (d1 == d2) return ;
  assert(!exact && finite(d1) && finite(d2)) ;

  //-- from OsiFloatEqual.hpp
  static const double epsilon = 1.e-10 ;
  double tol = std::max(fabs(d1), fabs(d2)) + 1 ;
  double diff = fabs(d1 - d2) ;

  assert(diff <= tol * epsilon) ;
}


/*! \brief Compare two integer values for equality. */

void ODSI::assert_same (int i1, int i2, bool exact)

{ assert(exact && i1 == i2) ; }


/*! \brief Byte-wise equality comparison of a structure

  Byte-for-byte comparison of a structure.
*/

template<class T>
  void ODSI::assert_same (const T& t1, const T& t2, bool exact)


{ assert(exact) ;                // be explicit
  if (&t1 == &t2) return ;
  assert(memcmp(&t1, &t2, sizeof(T)) == 0) ;
}


/*! \brief Element-wise comparison of two primitive arrays. */

template<class T>
  void ODSI::assert_same (const T* t1, const T* t2, int n, bool exact)

{ if (t1 == t2) return ;
  for (int i=0 ; i<n ; i++) assert_same(t1[i], t2[i], exact) ;
}


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
		inv_vec<basisel_struct>(b2.el),size) == 0) ;
}


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
  assert(c1.bnd == c2.bnd) ;
}


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

  //-- from consys.h::consys_struct

  // name can differ
  assert(!exact || c1.nme == c2.nme) ;
  assert(c1.parts == c2.parts) ;
  assert(c1.opts == c2.opts) ;
  assert(c1.varcnt == c2.varcnt) ;
  assert(c1.archvcnt == c2.archvcnt) ;
  assert(c1.logvcnt == c2.logvcnt) ;
  // FIXME loadproblem does not load these info
  assert(!exact || c1.intvcnt == c2.intvcnt) ;
  assert(!exact || c1.binvcnt == c2.binvcnt) ;
  assert(c1.maxcollen == c2.maxcollen) ;
  assert(c1.maxcolndx == c2.maxcolndx) ;
  assert(c1.concnt == c2.concnt) ;
  assert(c1.archccnt == c2.archccnt) ;
  assert(c1.cutccnt == c2.cutccnt) ;
  assert(c1.maxrowlen == c2.maxrowlen) ;
  assert(c1.maxrowndx == c2.maxrowndx) ;
  // capacity can differ
  assert(!exact || c1.colsze == c2.colsze) ;
  assert(!exact || c1.rowsze == c2.rowsze) ;
  // name can differ
  assert(!exact || c1.objnme == c2.objnme) ;
  assert(c1.objndx == c2.objndx) ;
  assert(c1.xzndx == c2.xzndx) ;
  
  int var_count = ODSI::idx(c1.varcnt) ;
  int con_count = ODSI::idx(c1.concnt) ;

  assert_same(c1.obj, c2.obj, var_count, exact) ;
  // FIXME loadproblem does not load these info
  if (exact) assert_same(c1.vtyp, c2.vtyp, var_count, true) ;
  assert_same(c1.vub, c2.vub, var_count, exact) ;
  assert_same(c1.vlb, c2.vlb, var_count, exact) ;
  assert_same(c1.rhs, c2.rhs, con_count, exact) ;
  assert_same(c1.rhslow, c2.rhslow, con_count, exact) ;
  assert_same(c1.ctyp, c2.ctyp, con_count, true) ;
  assert_same(c1.cub, c2.cub, con_count, exact) ;
  assert_same(c1.clb, c2.clb, con_count, exact) ;
  // check attvhdr_struct* attvecssrc
}


/*! \brief  Verify copy of a dylp lp problem (lpprob_struct)

  Essentially a matter of checking each individual component. Exact requires
  precise equivalence. Inexact allows for differences in the allocated capacity
  of the structures.
*/

void ODSI::assert_same (const lpprob_struct& l1, const lpprob_struct& l2,
			bool exact)

{ if (&l1 == &l2) return ;

  //-- from dylp.h::lpprob_struct

  assert(l1.ctlopts == l2.ctlopts) ;
  assert(l1.phase == l2.phase) ;
  assert(l1.lpret == l2.lpret) ;
  assert(l1.obj == l2.obj) ;
  assert(l1.iters == l2.iters) ;
  // capacity can differ
  assert(!exact || l1.colsze == l2.colsze) ;
  assert(!exact || l1.rowsze == l2.rowsze) ;

  int min_col = idx(std::min(l1.colsze, l2.colsze)) ;
  int min_row = idx(std::min(l1.rowsze, l2.rowsze)) ;

  assert_same(*l1.consys, *l2.consys, exact) ;
  assert_same(*l1.basis, *l2.basis, exact) ;
  assert_same(l1.status, l2.status, min_col, exact) ;
  assert_same(l1.x, l2.x, min_row, exact) ;
  assert_same(l1.y, l2.y, min_row, exact) ;
}


/*! \brief Verify copy of an ODSI object.

  Verify equivalence by checking each component in turn. Rely on the
  CoinPackedMatrix equivalence check to verify equivalence of the cached
  copies.
*/

void ODSI::assert_same (const OsiDylpSolverInterface& o1,
			const OsiDylpSolverInterface& o2, bool exact)

{ if (&o1 == &o2) return ;

  assert_same(*o1.lpprob, *o2.lpprob, exact) ;

  assert(!o1.lpprob || o1.consys == o1.lpprob->consys) ;
  assert(!o2.lpprob || o2.consys == o2.lpprob->consys) ;
  const CoinPackedMatrix* m1 = o1.getMatrixByCol() ;
  const CoinPackedMatrix* m2 = o2.getMatrixByCol() ;
  assert(!m1 || m1->isEquivalent(*m2)) ;

  assert_same(*o1.options, *o2.options, true) ;
  assert_same(*o1.statistics, *o2.statistics, true) ;
  assert_same(*o1.tolerances, *o2.tolerances, true) ;

  assert(o1.lp_retval == o2.lp_retval) ;
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

/*! \brief Given `file' or `file.ext1', construct `file.ext2'

  A utility routine for building related file names. Either of ext1 or ext2
  can be null or "" (the null string), in which case nothing is stripped or
  added, respectively.
*/

std::string ODSI::make_filename (const std::string filename, 
				 const char *ext1, const char *ext2)

{ string basename ;
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
  basename = filename ;
  if (ext1 && strlen(ext1) > 0)
  { string::size_type ext1pos = filename.rfind(ext1str) ;
    if (ext1pos < filename.npos)
    { basename = filename.substr(0,ext1pos) ; } }
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

   loadProblem and assignProblem will preserve existing options and
   tolerances across the problem change; if none exist, defaults will be
   installed. All variables are considered to be continuous. assignProblem
   destroys its arguments; loadProblem leaves them unchanged.

   When reading problem.mps, readMps intialises options and tolerances to
   default values and then looks for an options (.spc) file problem.spc and
   uses it to modify the default parameter settings. readMps will recognize
   continuous, binary, and general integer variables.
*/
//@{


/*!
  Read a problem definition in MPS format from the file basename.extension.

  This routine will also look in the same directory as the MPS file for a
  dylp options file `basename.spc'.  Options and tolerances are initialised to
  defaults and then modified according to problem.spc. `extension' can be
  null or "" (the null string).

  NOTE: dylp does not support direct reading of compressed files.
*/

int ODSI::readMps (const char* basename, const char* extension)

{ string mps,base,ext ;
  char *dotp ;

  base = basename ;
  ext = extension ;
  if (ext != "")
  { mps = make_filename(base.c_str(),"",ext.c_str()) ; }
  else
  { mps = base ; }

  destruct_problem() ;
  assert(!lpprob && !consys) ;

  construct_options() ;
  string spc = make_filename(base.c_str(),ext.c_str(),".spc") ;
  dylp_controlfile(spc.c_str(),true,false) ;
  bool r = dy_mpsin(mps.c_str(),&consys) ;
  assert(r) ;
  osi_probname = consys->nme ;

  construct_lpprob() ;

  return (0) ;
}


/*!
  Write a problem to the file basename.extension in MPS format.

  `extension' can be null or "" (the null string).

  NOTE: dylp does not support direct writing of compressed files.
*/

void ODSI::writeMps (const char *basename, const char *extension) const

{ CoinMpsIO mps ;
  int n = getNumCols(),
      m = getNumRows() ;
  char *vartyp = new char[n] ;
  typedef char *charp ;
  char **colnames = new charp[n],
       **rownames = new charp[n] ;
  
  string filename = basename ;
  if (extension && strlen(extension) > 0)
  { filename += '.' ;
    filename += extension ; }
  
  int i,j,errs ;

  for (j = 0 ; j < n ; j++) vartyp[j] = isInteger(j) ;

  for (i = 0 ; i < m ; i++)
    rownames[i] = consys_nme(consys,'c',idx(i),false,0) ;
  
  for (j = 0 ; j < n ; j++)
    colnames[j] = consys_nme(consys,'v',idx(j),false,0) ;

  mps.setMpsData(*getMatrixByRow(),
		 getInfinity(),
		 getColLower(),
		 getColUpper(),
		 getObjCoefficients(),
		 vartyp,
		 getRowLower(),
		 getRowUpper(),
		 colnames,
		 rownames) ;

/*
  We really need to work on symbolic names for these magic numbers.
*/
  errs = mps.writeMps(filename.c_str(),0,4,2) ;
/*
  Free up the arrays we allocated.
*/
  delete[] vartyp ;
  delete[] colnames ;
  delete[] rownames ;

  return ; }


/*!
  Write a problem to the file basename.extension in MPS format.

  `extension' can be null or "" (the null string).

  NOTE: dylp does not support direct writing of compressed files.
*/

int ODSI::writeMps (const char *filename,
		    const char **rownames, const char **colnames,
		    int formatType, int numberAcross) const

{ CoinMpsIO mps ;
  int n = getNumCols(),
      m = getNumRows() ;
  char *vartyp = new char[n] ;
  
  int j,errs ;

  for (j = 0 ; j < n ; j++) vartyp[j] = isInteger(j) ;

  mps.setMpsData(*getMatrixByRow(),
		 getInfinity(),
		 getColLower(),
		 getColUpper(),
		 getObjCoefficients(),
		 vartyp,
		 getRowLower(),
		 getRowUpper(),
		 colnames,
		 rownames) ;

/*
  We really need to work on symbolic names.
*/
  errs = mps.writeMps(filename,0,formatType,numberAcross) ;
/*
  Free up the arrays we allocated.
*/
  delete[] vartyp ;

  return (errs) ; }


/*!
  Constraints are described using an OSI packed matrix and upper and lower
  bounds for the left-hand-side.  The routine's parameters are destroyed once
  they have been copied to dylp structures. Existing options and tolerances
  are preserved if they exist; defaults are used otherwise.
*/

inline void
ODSI::assignProblem (CoinPackedMatrix*& matrix,
		     double*& col_lower, double*& col_upper, double*& obj, 
		     double*& row_lower, double*& row_upper)

{ loadProblem (*matrix,col_lower,col_upper,obj,row_lower,row_upper) ;
  delete matrix ; matrix = 0 ;
  delete [] col_lower ; col_lower = 0 ;
  delete [] col_upper ; col_upper = 0 ;
  delete [] obj ; obj = 0 ;
  delete [] row_lower ; row_lower = 0 ;
  delete [] row_upper ; row_upper = 0 ;
}


/*!
  Constraints are described using an OSI packed matrix, constraint sense,
  right-hand-side, and (optional) range. The routine's parameters are
  destroyed once they have been copied to dylp structures. Existing options
  and tolerances are preserved if they exist; defaults are used otherwise.
*/

inline void
ODSI::assignProblem (CoinPackedMatrix*& matrix, 
		     double*& lower, double*& upper, double*& obj, 
		     char*& sense, double*& rhs, double*& range)

{ loadProblem (*matrix, lower, upper, obj, sense, rhs, range) ;
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

inline void ODSI::loadProblem
(const CoinPackedMatrix& matrix,
 const double* col_lower, const double* col_upper, const double* obj,
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

inline void ODSI::loadProblem
(const CoinPackedMatrix& matrix, 
 const double* col_lower, const double* col_upper, const double* obj,
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

inline void ODSI::loadProblem
(const int colcnt, const int rowcnt,
 const int *start, const int *index, const double *value,
 const double* col_lower, const double* col_upper, const double* obj,
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

inline void ODSI::loadProblem
(const int colcnt, const int rowcnt,
 const int *start, const int *index, const double *value,
 const double* col_lower, const double* col_upper, const double* obj,
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
	monotonically.
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

  \todo Work out a general scheme for supporting dylp parameters.
*/

//@{

inline double ODSI::getInfinity () const { return DYLP_INFINITY ; }

bool ODSI::setDblParam (OsiDblParam key, double value)
/*
  All supported double parameters belong to the dylp tolerances structure.
  If there's no tolerance structure, create one, then change the appropriate
  entry.
*/
{ if (!tolerances)
  { lpopts_struct *throwaway = 0 ;
    dy_defaults(&throwaway,&tolerances) ;
    delete throwaway ; }
  
  switch (key)
  { case OsiDualTolerance:
    { osi_dual_tolerance = value ;
      tolerances->dfeas_scale = osi_dual_tolerance/tolerances->cost ;
      break ; }
    case OsiObjOffset:
    { osi_obj_offset = value ;
      break ; }
    case OsiPrimalTolerance:
    { osi_primal_tolerance = value ;
      tolerances->pfeas_scale = osi_primal_tolerance/tolerances->zero ;
      break ; }
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    default:
    { return false ; } }
  
  return (true) ; }


bool ODSI::getDblParam (OsiDblParam key, double& value) const
/*
  All supported double parameters are contained in the dylp tolerances
  structure.
*/
{ if (!tolerances) return (false) ;

  switch (key)
  { case OsiDualTolerance:
    { value = osi_dual_tolerance ;
      break ; }
    case OsiObjOffset:
    { value = osi_obj_offset ;
      break ; }
    case OsiPrimalTolerance:
    { value = osi_primal_tolerance ;
      break ; }
    case OsiDualObjectiveLimit:
    case OsiPrimalObjectiveLimit:
    default:
    { return false ; } }

  return (true) ; }


bool ODSI::setIntParam (OsiIntParam key, int value)
/*
  All supported integer parameters are contained in the options structure.
  If there's no options structure, create one, then change the appropriate
  entry.

  Dylp's iteration limit is normally enforced on a per-phase basis, with an
  overall limit of 3*(per-phase limit). Hence the factor of 3.
 */
{ if (!options)
  { lptols_struct *throwaway = 0 ;
    dy_defaults(&options,&throwaway) ;
    delete throwaway ; }

  switch (key)
  { case OsiMaxNumIteration:
    { osi_iterlim = value ;
      options->iterlim = value/3 ;
      break ; }
    case OsiMaxNumIterationHotStart:
    { osi_hotiterlim = value ;
      break ; }
    default:
    { return false ; } }
  
  return (true) ; }


bool ODSI::getIntParam (OsiIntParam key, int& value) const
/*
  All supported integer parameters are contained in the dylp options structure,
  so if it's missing, we're toast.
*/
{ if (!options) return (false) ;

  switch (key)
  { case OsiMaxNumIteration:
    { value = osi_iterlim ;
      break ; }
    case OsiMaxNumIterationHotStart:
    { value = osi_hotiterlim ;
      break ; }
    default:
    { return (false) ; } }
  
  return (true) ;  }


bool ODSI::setStrParam (OsiStrParam key, const std::string& value)
/*
  The only string parameter so far is OsiProbName. While it can be set and
  retrieved here, it will be overridden whenever an MPS file is loaded.
*/

{ switch (key)
  { case OsiProbName:
    { osi_probname = value ;
      break ; }
    default:
    { return (false) ; } }
  
  return (true) ; }

bool ODSI::getStrParam (OsiStrParam key, std::string& value) const

{ switch (key)
  { case OsiProbName:
    { value = osi_probname ;
      break ; }
    case OsiSolverName:
    { value="dylp" ;
      break ; }
    default:
    { return (false) ; } }
  
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
  The first action is to initialise the basis package, if it hasn't already
  been done. Next, check whether dylp is owned by another solver instance.
  If so, call detach_dylp to detach it.  The routine then sets various
  other dylp parameters as appropriate and calls dylp_dolp to solve the lp.

  Forcing fullsys here is a good choice in terms of the performance of the
  code, but it does reduce flexibility.
*/
void ODSI::initialSolve ()

{ bool fullsys_tmp,forcecold_tmp ;

  assert(lpprob && consys && options && tolerances && statistics) ;

  if (basis_ready == false)
  { int count = static_cast<int>((1.5*getNumRows()) + 2*getNumCols()) ;
    dy_initbasis(count,options->factor,0) ;
    basis_ready = true ; }
  if (dylp_owner != 0 && dylp_owner != this) dylp_owner->detach_dylp() ;

  if (isactive(local_logchn))
  { logchn = local_logchn ;
    gtxecho = local_gtxecho ; }
  lpprob->phase = dyINV ;
  forcecold_tmp = options->forcecold ;
  options->forcecold = true ;
  fullsys_tmp = options->fullsys ;
  options->fullsys = true ;
  lp_retval = dylp_dolp(lpprob,options,tolerances,statistics) ;
  if (flgon(lpprob->ctlopts,lpctlDYVALID)) dylp_owner = this ;
  options->fullsys = fullsys_tmp ;
  options->forcecold = forcecold_tmp ;
  destruct_row_cache() ;
  destruct_col_cache() ;
}

//@} // ColdStartMethods



/*! \defgroup SolverTerm Methods Returning Solver Termination Status */
//@{


inline int ODSI::getIterationCount () const

{ return (lpprob)?lpprob->iters:0 ; }


inline bool ODSI::isIterationLimitReached () const

{ return lp_retval == lpITERLIM ; }


inline bool ODSI::isProvenOptimal () const

{ return lp_retval == lpOPTIMAL ; }


inline bool ODSI::isProvenPrimalInfeasible () const

{ return lp_retval == lpINFEAS ; }


/*!
  Aka primal unbounded.
*/

inline bool ODSI::isProvenDualInfeasible () const

{ return lp_retval == lpUNBOUNDED ; }


/*!
  Returns true if dylp abandoned the problem. This could be due to numerical
  problems (accuracy check or singular basis), stalling, unexpected loss of
  feasibility, inability to allocate space, or other fatal internal error.
*/

inline bool ODSI::isAbandoned () const

{ if (lp_retval == lpACCCHK || lp_retval == lpSINGULAR ||
      lp_retval == lpSTALLED || lp_retval == lpNOSPACE ||
      lp_retval == lpLOSTFEAS || lp_retval == lpPUNT || lp_retval == lpFATAL)
    return true ;
  else
    return false ; }


//@} // SolverTerm



/*! \defgroup SolnInfo Methods to Get Solution Information */
//@{

/*!
  Returns +/- infinity (minimisation/maximisation) if there is no lp solution
  associated with this ODSI object.
*/

inline double ODSI::getObjValue () const
/*
  Of course, we need to correct for max/min.
*/
{ return (lpprob)?obj_sense*lpprob->obj:obj_sense*DYLP_INFINITY ; }


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
  Dylp returns a vector x of basic variables, in basis order, and a status
  vector. To assemble a complete primal solution, it's necessary to construct
  one vector that merges the two sources. The result is cached.

  NOTE the caution above about superbasic (vstatSB) variables. Nonbasic free
  variables (vstatNBFR) occur under the same termination conditions, but are
  legitimately zero. We really should return NaN for superbasics, using the
  function <limits>:std::numeric_limits::signaling_NaN(), but it's just not
  worth the hassle. Sun Workshop C++ (Rogue Wave Software) has it, but Gnu C++
  prior to v3 doesn't even have <limits>, and the local v3 installation's
  <limits> points to another, non-existent header.
*/
{ if (_col_x) return _col_x ;

  if (!(lpprob && lpprob->status && lpprob->x)) return (0) ;

  int n = getNumCols() ;
  flags statj ;
  double* solution = new double[n] ;

/*
  Walk the status vector, taking values for basic variables from lpprob.x and
  values for nonbasic variables from the appropriate bounds vector.
*/
  for (int j = 0 ; j < n ; j++)
  { statj = lpprob->status[idx(j)] ;
    if (((int) statj) < 0)
    { int k = -((int) statj) ;
      solution[j] = lpprob->x[k] ; }
    else
    { switch (statj)
      { case vstatNBLB:
	case vstatNBFX:
	{ solution[j] = consys->vlb[idx(j)] ;
	  break ; }
	case vstatNBUB:
	{ solution[j] = consys->vub[idx(j)] ;
	  break ; }
	case vstatNBFR:
	{ solution[j] = 0 ;
	  break ; }
	case vstatSB:
	{ solution[j] = -DYLP_INFINITY ;
	  break ; } } } }

  _col_x = solution ;
  return (_col_x) ; }


/*!
  Return a cached vector of dual variable values. If there is no cached copy,
  construct one and cache it.
*/

const double* ODSI::getRowPrice () const
/*
  Dylp reports dual variables in basis order, so we need to write a vector
  with the duals dropped into position by row index.
*/

{ if (_row_price) return _row_price ;

  if (!(lpprob && lpprob->basis && lpprob->y)) return (0) ;

  int n = getNumRows() ;
  double* price = new double[n] ;
  basis_struct* basis = lpprob->basis ;

  memset(price,0,n*sizeof(double)) ;

  for (int k = 1 ; k <= basis->len ; k++)
  { int i = basis->el[k].cndx ;
    price[inv(i)] = lpprob->y[k] ; }

  _row_price = price ;
  return price ;
}


/*!
  Return a cached vector of row activity. If there is no cached copy,
  calculate and cache lhs<i> = a<i>x.
*/

const double *ODSI::getRowActivity () const
/*
  This routine calculates the value of the left-hand-side of a constraint. The
  algorithm is straightforward: Retrieve the primal solution, then walk the
  variables, adding the contribution of each non-zero variable.
*/
{ assert(consys) ;
/*
  If we have a cached copy, we're done.
*/
  if (_row_lhs) return (_row_lhs) ;
/*
  Retrieve the primal solution.
*/
  const double *x = getColSolution() ;
  if (!x) return (0) ;
/*
  Create a vector to hold ax and clear it to 0.
*/
  double *lhs = new double[consys->concnt] ;
  memset(lhs,0,consys->concnt*sizeof(double)) ;
/*
  Walk the primal solution vector, adding in the contribution from non-zero
  variables.
*/
  pkvec_struct *aj = 0 ;
  bool r ;
  for (int j = 0 ; j < consys->varcnt ; j++)
  { if (x[j] != 0)
    { r = consys_getcol_pk(consys,idx(j),&aj) ;
      if (!r)
      { delete[] lhs ;
	if (aj) pkvec_free(aj) ;
	return (0) ; }
      for (int l = 0 ; l < aj->cnt ; l++)
      { int i = inv(aj->coeffs[l].ndx) ;
	lhs[i] += x[j]*aj->coeffs[l].val ; } } }
  if (aj) pkvec_free(aj) ;
/*
  Groom the vector to eliminate tiny values.
*/
  for (int i = 0 ; i < consys->concnt ; i++)
    setcleanzero(lhs[i],tolerances->zero) ;
/*
  Cache the result and return.
*/
  _row_lhs = lhs ;

  return (_row_lhs) ; }


/*!
  Return a cached vector of reduced costs. If there is no cached copy,
  calculate and cache cbar<j> = c<j> - ya<j>.
*/

const double *ODSI::getReducedCost () const
/*
  Calculate the reduced cost as cbar<j> = c<j> - ya<j>. It's tempting to use
  dy_pricenbvars, but that routine dives into the internal data structures
  of dylp. In an environment like COIN, where some other object might have
  usurped dylp's static structures, relying on object data is better.

  Dylp returns a vector of duals that matches the basis. The routine walks
  the basis, pulling out the constraint in each basis position and adding
  the contribution to the reduced costs.

  We need to compensate for maximisation once the calculation is complete.
*/

{ assert(lpprob && lpprob->basis && lpprob->y && consys && consys->obj) ;
/*
  Check for a cached value.
*/
  if (_col_cbar) return (_col_cbar) ;
/*
  Initialise the reduced cost vector by copying in the objective function.
*/
  double *cbar = 0 ;
  CLONE_VEC(double,INV_VEC(double,consys->obj),cbar,consys->varcnt);  
/*
  Now start to walk through the basis. For each entry, check that the
  corresponding dual (y[k] == y<i>) is nonzero. If so, pull the row i for this
  basis entry and add the contribution -y<i>a<ij> to c<j> for each non-zero
  coefficient.
*/
  int i,j ;
  basis_struct *basis = lpprob->basis ;
  double *y = lpprob->y ;
  pkvec_struct *ai = 0 ;
  bool r ;

  for (int k = 1 ; k <= basis->len ; k++)
  { i = basis->el[k].cndx ;
    if (y[k] != 0)
    { r = consys_getrow_pk(consys,i,&ai) ;
      if (!r)
      { delete[] cbar ;
	if (ai) pkvec_free(ai) ;
	return (0) ; }
      for (int l = 0 ; l < ai->cnt ; l++)
      { j = inv(ai->coeffs[l].ndx) ;
	cbar[j] -= y[k]*ai->coeffs[l].val ; } } }
  if (ai) pkvec_free(ai) ;
/*
  Groom the vector to eliminate tiny values.
*/
  for (j = 0 ; j < consys->varcnt ; j++)
    setcleanzero(cbar[j],tolerances->cost) ;
/*
  Compensate for maximisation, so that the caller sees the correct sign.
*/
  if (obj_sense < 0)
  { for (j = 0 ; j < consys->varcnt ; j++) cbar[j] = -cbar[j] ; }
/*
  Cache the result and return.
*/
  _col_cbar = cbar ;

  return (cbar) ; }



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
    return consys->vtyp[idx(i)] == vartypBIN ; }


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
    return consys->vtyp[idx(i)] == vartypINT ; }


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
    return consys->vtyp[idx(i)] == vartypCON ; }


/*! Returns a pointer to the dylp data structure, with indexing shift. */

inline const double* ODSI::getColLower () const

{ if (!consys || !consys->vlb)
    return (0) ;
  else
    return INV_VEC(double,consys->vlb) ; }


/*! Returns a pointer to the dylp data structure, with indexing shift. */

inline const double* ODSI::getColUpper () const

{ if (!consys || !consys->vub)
    return (0) ;
  else
    return INV_VEC(double,consys->vub) ; }


inline int ODSI::getNumCols () const

{ if (!consys)
    return 0 ;
  else
    return consys->varcnt ; }


inline int ODSI::getNumElements () const

{ if (!consys)
    return (0) ; 
  else
    return consys->mtx.coeffcnt ; }


inline int ODSI::getNumRows () const

{ if (!consys)
    return 0 ;
  else
    return consys->concnt ; }


/*! 
  Creates a cached copy, if it doesn't already exist, compensating for the
  objective sense. Returns a pointer to the cached copy.
*/

inline const double* ODSI::getObjCoefficients () const

{ if (!consys || !consys->obj) return (0) ;

  if (_col_obj) return _col_obj ;

  CLONE_VEC(double,INV_VEC(double,consys->obj),_col_obj,consys->varcnt) ;
  for (int i = 0 ; i < consys->varcnt ; i++) _col_obj[i] *= obj_sense ;

  return (_col_obj) ; }


inline double ODSI::getObjSense () const

{ return obj_sense ; }


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

{ if (_row_sense) return _row_sense ;
  
  int n = getNumRows() ;
  char* sense = new char[n] ;
  const contyp_enum* ctyp = INV_VEC(contyp_enum,consys->ctyp) ;

  for (int i = 0 ; i < n ; i++)
    sense[i] = type_to_sense(ctyp[i]) ;

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

{ if (_row_rhs) return _row_rhs ;
  
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
  algorithm. Alternating addition of constraints and variables means that the
  objective does not change monotonically.
*/

bool ODSI::isDualObjectiveLimitReached () const

{ UNSUPPORTED ;
  return false ; }


/*!
  Effectively unsupportable. See
  \link OsiDylpSolverInterface::isDualObjectiveLimitReached
	ODSI:isDualObjectiveLimitReached\endlink.
*/
bool ODSI::isPrimalObjectiveLimitReached () const

{ UNSUPPORTED ;
  return false ; }


/*!
  Supportable but dylp will pay no attention.
*/

void ODSI::setColSolution (const double* solution)

{ UNSUPPORTED ; }


/*!
  Supportable but dylp will pay no attention.
*/

void ODSI::setRowPrice (const double* price)

{ UNSUPPORTED ; }


/*!
  \todo Will require modification to dylp --- currently it does not report
	this information.
*/

vector<double*> ODSI::getDualRays (int) const

{ UNSUPPORTED ;
  return vector<double*>() ; }


/*!
  \todo Will require modification to dylp --- currently it does not report
	this information.
*/

vector<double*> ODSI::getPrimalRays (int) const

{ UNSUPPORTED ;
  return vector<double*>() ; }


/*!
  \todo Someday, maybe, but it's hardly worth trying to shoehorn bonsaiG into
	this framework. Better to start from scratch with OSI as the target,
	or use bonsaiG standalone.
*/

void ODSI::branchAndBound ()

{ UNSUPPORTED ; }


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
  constraints.  For ease of use, constraint status handled just like variable
  status. There's an array, one entry per constraint, coded as
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

  assert(lpprob && options && consys && consys->vlb && consys->vub) ;
/*
  Use a dynamic cast here to make sure we have an OsiDylpWarmStartBasis.
*/
  const OsiDylpWarmStartBasis *wsb =
			dynamic_cast<const OsiDylpWarmStartBasis *>(ws) ;
  if (!wsb) return (false) ;
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
  x, y, and actvars arrays have been resized if needed. If there's an options
  structure in existence, signal a warm start, and then we're done.
*/
  options->forcecold = false ;
  options->forcewarm = true ;
  
  return (true) ; }


void ODSI::resolve ()
/*
  This routine simply calls the solver, after making sure that the forcecold
  option is turned off. If we're reoptimising, then the basis should be ready
  and we should have warm start information. We do need to make sure the
  solver is ours to use.
*/
{ assert(lpprob && lpprob->basis && lpprob->status && basis_ready &&
	 consys && options && tolerances && statistics) ;

  if (dylp_owner != 0 && dylp_owner != this) dylp_owner->detach_dylp() ;

  if (isactive(local_logchn))
  { logchn = local_logchn ;
    gtxecho = local_gtxecho ; }

  lpprob->phase = dyINV ;
  options->forcecold = false ;

  lp_retval = dylp_dolp(lpprob,options,tolerances,statistics) ;  
  if (flgon(lpprob->ctlopts,lpctlDYVALID)) dylp_owner = this ;

  destruct_col_cache() ;
  destruct_row_cache() ; }

//@} // WarmStart



/*! \defgroup HotStartMethods Hot Start Methods
    \brief Methods for solving an LP from a hot start.

  In the absence of instructions to the contrary, dylp will always go for a
  hot start. The obligation of the caller is simply to insure that dylp's
  internal static data structures are intact from the previous call. In the
  context of OSI, this means no intervening calls to the solver by some
  other ODSI object.

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

{ assert(lpprob && options) ;

  assert(flgon(lpprob->ctlopts,lpctlDYVALID)) ;

  options->forcecold = false ;
  options->forcewarm = false ;

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

  assert(lpprob && lpprob->basis && lpprob->status && basis_ready &&
	 consys && options && tolerances && statistics) ;

  if (isactive(local_logchn))
  { logchn = local_logchn ;
    gtxecho = local_gtxecho ; }

  lpprob->phase = dyINV ;
  options->forcecold = false ;
  options->forcewarm = false ;
  if (osi_hotiterlim > 0)
  { tmp_iterlim = options->iterlim ;
    options->iterlim = osi_hotiterlim/3 ; }

  lp_retval = dylp_dolp(lpprob,options,tolerances,statistics) ;

  if (tmp_iterlim > 0) options->iterlim = tmp_iterlim ;

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
    main_lpopts = options ;
    main_lptols = tolerances ;
    bool r = (process_cmds(silent) != 0) ;
    (void) closefile(cmdchn) ;
    assert(r == cmdOK) ; }

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
  local_gtxecho = echo ; }

//@} // DylpMethods
