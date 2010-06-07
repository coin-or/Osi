# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Osi],[trunk],[osi@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/OsiAuxInfo.cpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# A bit of initial setup (ADDLIBS, FADDLIBS)
AC_COIN_PROJECTDIR_INIT(Osi)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Osi)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_HAS_MODULE(CoinUtils, [coinutils = trunk])
if test $coin_has_coinutils != yes ; then
  #if we have to use the fallback routines, then obviously it will not find CoinUtils dependencies
  #thus, we check whether ThirdParty/Blas and ThirdParty/Lapack have been setup and
  #assume that if present then they are dependencies of CoinUtils
  AC_COIN_HAS_MODULE_FALLBACK(Blas,   [coinblas],   [../ThirdParty/Blas],   [required=0 dodefine=0 doconditional=0])
  AC_COIN_HAS_MODULE_FALLBACK(Lapack, [coinlapack], [../ThirdParty/Lapack], [required=0 dodefine=0 doconditional=0])
fi
AC_COIN_HAS_MODULE_FALLBACK(CoinUtils, [coinutils], [../CoinUtils])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
#so if we had to do find Blas and Lapack separately, then add their linker flags to the CoinUtils flags
#actually the libs may already be included in libCoinUtils, but we may still need the fortran library flags
COINUTILS_LIBS="$COINUTILS_LIBS $BLAS_LIBS $LAPACK_LIBS"

AC_COIN_HAS_MODULE(Sample, [coindatasample],   [required=0], [../Data/Sample])
AC_COIN_HAS_MODULE(Netlib, [coindatanetlib],   [required=0], [../Data/Netlib])
AC_COIN_HAS_MODULE(Glpk,   [coinglpk >= 4.39], [required=0], [../ThirdParty/Glpk])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available. Cplex detection has been
# tested on a semi-regular basis, and Mosek detection most likely works. For
# the rest, it's anyone's guess.  -- lh, 080529 --
# fixed detection of Mosek and Xpress -- stefan, 091003 (linux32, gcc) --

AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_makeenv])
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSinit])
AC_COIN_HAS_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBloadenv])

#############################################################################
#                            Examples solver                                #
#############################################################################

# choose a solver interface that can be used in the examples
if test $coin_has_glpk = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGlpkSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$GLPK_CFLAGS"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiGlpk $GLPK_LIBS"
elif test $coin_has_cpx = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiCpxSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$CPXINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiCpx $CPXLIB"
elif test $coin_has_grb = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGrbSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$GRBINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiGrb $GRBLIB"
elif test $coin_has_msk = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiMskSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$MSKINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiMsk $MSKLIB"
elif test $coin_has_xpr = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiXprSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$XPRINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiXpr $XPRLIB"
fi
AC_SUBST(OSI_EXAMPLES_SOLVER_NAME)
AC_SUBST(OSI_EXAMPLES_SOLVER_CFLAGS)
AC_SUBST(OSI_EXAMPLES_SOLVER_LIBS)

#############################################################################
#		  OSI version number                                        #
#############################################################################

AC_DEFINE_UNQUOTED(OSI_VERSION, "$PACKAGE_VERSION", "Version of OSI")

#############################################################################
#                 Configuration options for individual OSIs                 #
#############################################################################

# If any of the tests performed in this section actually require the presence
# of the solver (file presence, link checks, etc.) be sure to guard the call.

# we assume that GLPK is not too old
AC_DEFINE(GLPK_HAS_INTOPT, [1], [Define to 1 if GLPK has the advanced B&B solver lpx_intopt])

# For setting up osi-glpk.pc, it makes a difference whether Glpk is available because
# --with-glpk-lib has been specified or because a coinglpk.pc file was found.
# Normally, COIN_HAS_MODULE takes care of this, but since we set required=0 above,
# we have to find out here in which situation we are.
if test "$coin_has_glpk" = "yes" ; then
  AC_ARG_WITH([glpk-lib],[],
    [ OSIGLPKPC_LIBS=${GLPK_LIBS}
      OSIGLPKPC_CFLAGS=${GLPK_CFLAGS}
      OSIGLPKPC_REQUIRES=
    ],
    [ OSIGLPKPC_LIBS=
      OSIGLPKPC_CFLAGS=
      OSIGLPKPC_REQUIRES="coinglpk >= 4.39"
    ])
  AC_SUBST(OSIGLPKPC_LIBS)
  AC_SUBST(OSIGLPKPC_CFLAGS)
  AC_SUBST(OSIGLPKPC_REQUIRES)
fi

if test "$GXX" = "yes"; then
  case "$CXX" in
    icpc* | */icpc* | icc )
      ;;
    *)
      MSKCXXFLAGS="-Wno-long-long"
      ;;
  esac
fi
AC_SUBST(MSKCXXFLAGS)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/OsiCpx/Makefile
                 src/OsiGlpk/Makefile
                 src/OsiMsk/Makefile
                 src/OsiXpr/Makefile
                 src/OsiGrb/Makefile
                 src/unittests/Makefile
                 test/Makefile
                 osi.pc osi-uninstalled.pc
                 osi-unittests.pc osi-unittests-uninstalled.pc])

if test $coin_has_cpx = true ; then
  AC_CONFIG_FILES([src/OsiCpx/osi-cplex.pc
                   osi-cplex-uninstalled.pc:src/OsiCpx/osi-cplex-uninstalled.pc.in])
fi

if test $coin_has_glpk = yes ; then
  AC_CONFIG_FILES([src/OsiGlpk/osi-glpk.pc
                   osi-glpk-uninstalled.pc:src/OsiGlpk/osi-glpk-uninstalled.pc.in])
fi

if test $coin_has_grb = true ; then
  AC_CONFIG_FILES([src/OsiGrb/osi-gurobi.pc
                   osi-gurobi-uninstalled.pc:src/OsiGrb/osi-gurobi-uninstalled.pc.in])
fi

if test $coin_has_msk = true ; then
  AC_CONFIG_FILES([src/OsiMsk/osi-mosek.pc
                   osi-mosek-uninstalled.pc:src/OsiMsk/osi-mosek-uninstalled.pc.in])
fi

if test $coin_has_xpr = true ; then
  AC_CONFIG_FILES([src/OsiXpr/osi-xpress.pc
                   osi-xpress-uninstalled.pc:src/OsiXpr/osi-xpress-uninstalled.pc.in])
fi

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_osi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
