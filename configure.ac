# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Osi],[devel],[https://github.com/coin-or/Osi/issues/new],[coin-or-osi],[https://github.com/coin-or/Osi])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Osi which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test whether
# the package is actually there.
AC_CONFIG_SRCDIR(src/Osi/OsiSolverInterface.cpp)

# Do some initialization work (version numbers, change prefix default, ...)
AC_COIN_INITIALIZE

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.
AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.
AC_LANG_PUSH(C++)

# Initialize libtool.
AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h
AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

# Find components for OsiLib
AC_COIN_CHK_PKG(CoinUtils,[OsiLib OsiTest])
if test "$coin_has_coinutils" != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi

# Finalize libs and cflags for building OsiLib.
AC_COIN_FINALIZE_FLAGS(OsiLib)

# We just need the COIN_HAS_XXX and XXX_DATA variables for the data projects.
AC_COIN_CHK_PKG(Sample,,[coindatasample],[],dataonly)
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib],[],dataonly)

# Finalise lflags and cflags, then tack on information from OsiLib, which is
# not yet installed.
AC_COIN_FINALIZE_FLAGS(OsiCommonTestLib)
AC_COIN_CHK_HERE([OsiLib],[OsiCommonTestLib],[osi])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# These are third-party solvers with source distributions that can be
# downloaded and built. The general flow here is as above: Initialise the
# library, add external components, finalize the flags, and add any components
# being built here.
AC_COIN_CHK_PKG(Glpk,[OsiGlpkLib OsiTest],[coinglpk])

AC_LANG_PUSH(C++)
AC_COIN_CHK_LIBHDR(SoPlex,[OsiSpxLib OsiTest],[],[],[],
  [new soplex::DIdxSet()],
  [#include "soplex.h"
   #if SOPLEX_VERSION >= 400
   #error "OsiSpx requires SoPlex < 4.0"
   #endif
  ],
  [default_skip])
AC_LANG_POP(C++)

AC_COIN_CHK_LIBHDR(Cplex,[OsiCpxLib OsiTest],[-lcplex -lpthread -lm -ldl],[],[],
  [int err;
   #if CPX_VERSION >= 800
     CPXopenCPLEX(&err);
   #else
     CPXopenCPLEXdevelop(&err);
   #endif],
  [#include "cplex.h"])

AC_COIN_CHK_LIBHDR(Mosek,[OsiMskLib OsiTest],[-lmosek64 -lpthread],[],[],
  [MSKenv_t env_;
   #if MSK_VERSION_MAJOR >= 7
     MSK_makeenv(&env_,NULL);
   #else
     MSK_makeenv(&env_,NULL,NULL,NULL,NULL);
   #endif],
  [#include "mosek.h"
   #include <stdlib.h>])

AC_COIN_CHK_LIBHDR(Xpress,[OsiXprLib OsiTest],[-lxprs -lxprl],[],[],
  [XPRSinit(NULL)],
  [#include "xprs.h"
   #include <stdlib.h>])

AC_COIN_CHK_LIBHDR(Gurobi,[OsiGrbLib OsiTest],[-lgurobi -lpthread -lm],[],[],
  [GRBenv* env_;
   GRBloadenv(&env_, NULL)],
  [#include "gurobi_c.h"
   #include <stdlib.h>])

AC_COIN_FINALIZE_FLAGS([OsiGlpkLib OsiSpxLib OsiCpxLib OsiMskLib OsiXprLib OsiGrbLib OsiTest])

AC_COIN_CHK_HERE([OsiLib],[OsiGlpkLib OsiSpxLib OsiCpxLib OsiMskLib OsiXprLib OsiGrbLib],[osi])

#############################################################################
#                            Examples solver                                #
#############################################################################

# Choose a solver interface that can be used in the examples
if test "$coin_has_glpk" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGlpkSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSIGLPKLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSIGLPKLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-glpk
elif test "$coin_has_cplex" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiCpxSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSICPXLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSICPXLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-cplex
elif test "$coin_has_gurobi" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGrbSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSIGRBLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSIGRBLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-gurobi
elif test "$coin_has_mosek" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiMskSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSIMSKLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSIMSKLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-mosek
elif test "$coin_has_xpress" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiXprSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSIXPRLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSIXPRLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-xpress
elif test "$coin_has_soplex" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiSpxSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSISPXLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LFLAGS="$OSISPXLIB_LFLAGS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-soplex
else
  AC_MSG_WARN([No solver available, examples will not work.])
fi

AC_SUBST(OSI_EXAMPLES_SOLVER_NAME)
AC_SUBST(OSI_EXAMPLES_SOLVER_CFLAGS)
AC_SUBST(OSI_EXAMPLES_SOLVER_LFLAGS)
AC_SUBST(OSI_EXAMPLES_SOLVER_PCNAME)


#############################################################################
#                 Configuration options for individual OSIs                 #
#############################################################################

# If any of the tests performed in this section actually require the presence
# of the solver (file presence, link checks, etc.) be sure to guard the call.

# We assume that GLPK is not too old
AC_DEFINE(GLPK_HAS_INTOPT, [1],
  [Define to 1 if GLPK has the advanced B&B solver lpx_intopt])

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN([CoinUtils])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Osi/Makefile
                 src/OsiCpx/Makefile
                 src/OsiGlpk/Makefile
                 src/OsiMsk/Makefile
                 src/OsiXpr/Makefile
                 src/OsiGrb/Makefile
                 src/OsiSpx/Makefile
                 src/OsiCommonTest/Makefile
                 test/Makefile
                 osi.pc
                 osi-unittests.pc])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/Osi/config.h src/Osi/config_osi.h src/OsiCommonTest/config_osicommontest.h])

if test "$coin_has_cplex" = true ; then
  AC_CONFIG_FILES([osi-cplex.pc:src/OsiCpx/osi-cplex.pc.in])
  AC_CONFIG_HEADER([src/OsiCpx/config_osicpx.h])
fi

if test "$coin_has_glpk" = yes ; then
  AC_CONFIG_FILES([osi-glpk.pc:src/OsiGlpk/osi-glpk.pc.in])
  AC_CONFIG_HEADER([src/OsiGlpk/config_osiglpk.h])
fi

if test "$coin_has_gurobi" = yes ; then
  AC_CONFIG_FILES([osi-gurobi.pc:src/OsiGrb/osi-gurobi.pc.in])
  AC_CONFIG_HEADER([src/OsiGrb/config_osigrb.h])
fi

if test "$coin_has_mosek" = yes ; then
  AC_CONFIG_FILES([osi-mosek.pc:src/OsiMsk/osi-mosek.pc.in])
  AC_CONFIG_HEADER([src/OsiMsk/config_osimsk.h])
fi

if test "$coin_has_xpress" = yes ; then
  AC_CONFIG_FILES([osi-xpress.pc:src/OsiXpr/osi-xpress.pc.in])
  AC_CONFIG_HEADER([src/OsiXpr/config_osixpr.h])
fi

if test "$coin_has_soplex" = yes ; then
  AC_CONFIG_FILES([osi-soplex.pc:src/OsiSpx/osi-soplex.pc.in])
  AC_CONFIG_HEADER([src/OsiSpx/config_osispx.h])
fi

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
