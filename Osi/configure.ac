# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Osi],[split],[osi@list.coin-or.org])

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
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Osi)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_MODULE(CoinUtils,  [coinutils = 2.5BSPsplit])
AC_COIN_HAS_MODULE(Sample,     [coindatasample],   [required=0])
AC_COIN_HAS_MODULE(Netlib,     [coindatanetlib],   [required=0])
AC_COIN_HAS_MODULE(Glpk,       [coinglpk >= 4.39], [required=0])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available. Cplex detection has been
# tested on a semi-regular basis, and Mosek detection most likely works. For
# the rest, it's anyone's guess.  -- lh, 080529 --
# fixed detection of Mosek and Xpress -- stefan, 091003 (linux32, gcc) --

AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_makeenv])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSinit])
AC_COIN_HAS_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBloadenv])

#############################################################################
#		  Capability Tests					                                            #
#############################################################################

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
                 src/OsiFmp/Makefile
                 src/OsiGlpk/Makefile
                 src/OsiMsk/Makefile
                 src/OsiOsl/Makefile
                 src/OsiXpr/Makefile
                 src/OsiGrb/Makefile
                 src/unittests/Makefile
                 test/Makefile
                 osi_addlibs.txt
                 osi.pc osi-uninstalled.pc
                 osi-unittests.pc osi-unittests-uninstalled.pc])

if test $coin_has_cpx = true ; then
  AC_CONFIG_FILES([src/OsiCpx/osi-cplex.pc
                   osi-cplex-uninstalled.pc:src/OsiCpx/osi-cplex-uninstalled.pc.in])
fi

if test $coin_has_fmp = true ; then
  AC_CONFIG_FILES([src/OsiFmp/osi-fortmp.pc
                   osi-fortmp-uninstalled.pc:src/OsiFmp/osi-fortmp-uninstalled.pc.in])
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

if test $coin_has_osl = true ; then
  AC_CONFIG_FILES([src/OsiOsl/osi-osl.pc
                   osi-osl-uninstalled.pc:src/OsiOsl/osi-osl-uninstalled.pc.in])
fi

if test $coin_has_xpr = true ; then
  AC_CONFIG_FILES([src/OsiXpr/osi-xpress.pc
                   osi-xpress-uninstalled.pc:src/OsiXpr/osi-xpress-uninstalled.pc.in])
fi

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_osi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
