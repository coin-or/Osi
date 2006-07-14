
# AC_OSICBC_CONFIG(default_solver)

# Handles configuration of the underlying default solver in OsiCbc. The issue
# is that OsiCbc defines a default solver, used when the client does not
# specify a solver in the constructor. The default solver must, therefore,
# be present in the build. This macro checks that this is true, and sets the
# compile-time symbols OSICBC_DFLT_SOLVER, OSICBC_CLP_DFLT_SOLVER, and
# OSICBC_DFLT_SOLVER_HPP that control the build. The parameter default_solver
# should normally be clp, unless you're working on some other solver.

# Ideally, this macro would not require enumeration of solvers, but the
# effort required to avoid it is just not justified at present. One enumeration
# is hidden in AC_OSI_CANONICAL. The other is visible, the AM_CONDITIONAL
# list at the end.

# In an ideal world this macro would do absolutely nothing if Cbc is not
# present, but autotools is not an ideal world and we have to satisfy its
# requirements. In particular, the AM_CONDITIONAL macros need to execute or
# automake will complain. Really the only thing we need to suppress is the
# check that the default solver exists. All the rest is irrelevant when Cbc
# isn't present (hence OsiCbc will be configured but not actually compiled).

AC_DEFUN([AC_OSICBC_CONFIG],
[ 

# Process the with-osicbc-default-solver option.

  AC_ARG_WITH([osicbc-default-solver],
    AS_HELP_STRING([--with-osicbc-default-solver],
       [specify underlying solver for OsiCbc (default $1)]),
       [osicbc_with_solver=$withval],
       [osicbc_with_solver=$1])

# Get canonical forms of the solver name and an existence variable.

  AC_OSI_CANONICAL($osicbc_with_solver)

# Check that the requested solver is available. If we're not actually
# building OsiCbc, skip this check to avoid spurious failures in projects
# that don't have the default solver, Clp.

  if test $coin_has_cbc != unavailable &&
     test $coin_has_cbc != skipping; then
    if test $osi_exists_solver = no; then
      AC_MSG_ERROR([selected default solver $osicbc_with_solver is unavailable.
	Please select an available solver using the --with-osicbc-default-solver option.])
    fi
  fi

# State the result.

  AC_MSG_NOTICE([OsiCbc default solver is $osi_lc_solver])

# And set the configuration variables.

  AC_DEFINE_UNQUOTED([OSICBC_DFLT_SOLVER],
    [Osi${osi_mc_solver}SolverInterface],
    [define to the name of the default solver interface class, e.g.,
     OsiClpSolverInterface])
  AC_DEFINE_UNQUOTED([OSICBC_DFLT_SOLVER_HPP],
    ["Osi${osi_mc_solver}SolverInterface.hpp"],
    [define to the name of the .hpp file for the default solver interface
     class, e.g., "OsiClpSolverInterface.hpp" (include quotes)])
  if test $osi_mc_solver = "Clp"; then
    AC_DEFINE([OSICBC_DFLT_SOLVER_CLP],[1],
	      [define this symbol if clp is the default solver])
  fi

# Last but not least, we need automake conditionals.

  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_CLP],[test $osi_mc_solver = Clp])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_CPX],[test $osi_mc_solver = Cpx])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_DYLP],[test $osi_mc_solver = Dylp])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_FMP],[test $osi_mc_solver = Fmp])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_GLPK],[test $osi_mc_solver = Glpk])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_MSK],[test $osi_mc_solver = Msk])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_OSL],[test $osi_mc_solver = Osl])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SPX],[test $osi_mc_solver = Spx])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_SYM],[test $osi_mc_solver = Sym])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_VOL],[test $osi_mc_solver = Vol])
  AM_CONDITIONAL([OSICBC_DFLT_SOLVER_XPR],[test $osi_mc_solver = Xpr])
])
