# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_INIT([Osi],[trunk],[osi@list.coin-or.org],[coin-or-osi],
	[https://projects.coin-or.org/Osi])
AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Osi which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test whether
# the package is actually there.

AC_CONFIG_SRCDIR(src/Osi/OsiSolverInterface.cpp)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE(Osi)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options.

AC_COIN_PROG_CXX

# This is a C++ package, set the language accordingly.

AC_LANG_PUSH(C++)

# Initialize libtool

AC_COIN_PROG_LIBTOOL

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h

AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

AC_COIN_CHK_PKG(CoinUtils,[OsiLib])
if test "$coin_has_coinutils" != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHK_PKG(Sample,,[coindatasample])
AC_COIN_CHK_PKG(Netlib,,[coindatanetlib])

# These are third-party solvers with source distributions that can be
# downloaded and built.

AC_COIN_CHK_PKG(Glpk,[OsiGlpkLib],,[build])
AC_COIN_CHK_PKG(SoPlex,[OsiSpxLib],,[build])

#############################################################################
#                            Third party solvers                            #
#############################################################################

# These are proprietary third-party solvers. Cplex detection has been
# tested on a semi-regular basis, and Mosek detection most likely works. For
# the rest, it's anyone's guess.  -- lh, 080529 --
# fixed detection of Mosek and Xpress -- stefan, 091003 (linux32, gcc) --

# AC_COIN_CHECK_USER_LIBRARY([Cplex], [CPX],[cplex.h],   [CPXgetstat])
# AC_COIN_CHECK_USER_LIBRARY([Mosek], [MSK],[mosek.h],   [MSK_makeenv])
# AC_COIN_CHECK_USER_LIBRARY([Xpress],[XPR],[xprs.h],    [XPRSinit])
# AC_COIN_CHECK_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBloadenv])

AM_CONDITIONAL(COIN_HAS_CPX,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_MSK,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_XPR,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_GRB,[test 0 = 1])

# (lh, 180119) DEPENDENCY_LINKING sucks all dependencies into the
# library. Seems like something to be avoided if possible.

AM_CONDITIONAL(DEPENDENCY_LINKING,[test 0 = 1])

#############################################################################
#                            Examples solver                                #
#############################################################################

# Finalize libs and cflags for building OsiLib.
AC_COIN_FINALIZE_FLAGS(OsiLib)

# choose a solver interface that can be used in the examples

if test "$coin_has_glpk" = yes ; then
  AC_COIN_FINALIZE_FLAGS(OsiGlpkLib)
  OSI_EXAMPLES_SOLVER_NAME=OsiGlpkSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$GLPK_CFLAGS"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiGlpk $GLPK_LIBS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-glpk
elif test "$coin_has_cpx" = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiCpxSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$CPXINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiCpx $CPXLIB"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-cplex
elif test "$coin_has_grb" = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGrbSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$GRBINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiGrb $GRBLIB"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-gurobi
elif test "$coin_has_msk" = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiMskSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$MSKINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiMsk $MSKLIB"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-mosek
elif test "$coin_has_xpr" = true ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiXprSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="-I$XPRINCDIR"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiXpr $XPRLIB"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-xpress
elif test "$coin_has_soplex" = yes ; then
  AC_COIN_FINALIZE_FLAGS(OsiSpxLib)
  OSI_EXAMPLES_SOLVER_NAME=OsiSpxSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$SOPLEX_CFLAGS"
  OSI_EXAMPLES_SOLVER_LIBS="-lOsiSpx $SOPLEX_LIBS"
  OSI_EXAMPLES_SOLVER_PCNAME=osi-soplex
else
  AC_MSG_WARN([No solver available, examples will not work.])
fi

# adjust linker flags for (i)cl compiler
if test x$coin_cxx_is_cl = xtrue || test x$coin_cc_is_cl = xtrue ;
then
  [OSI_EXAMPLES_SOLVER_LIBS=`echo " $OSI_EXAMPLES_SOLVER_LIBS" | sed -e 's/ -l\([^ ]*\)/ lib\1.lib/g' -e 's/ -L\([^ ]*\)/ -libpath:\`$(CYGPATH_W) \1\`/g'`]
fi

AC_SUBST(OSI_EXAMPLES_SOLVER_NAME)
AC_SUBST(OSI_EXAMPLES_SOLVER_CFLAGS)
AC_SUBST(OSI_EXAMPLES_SOLVER_LIBS)
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

AC_COIN_DOXYGEN(CoinUtils)

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

if test "$coin_has_cpx" = true ; then
  AC_CONFIG_FILES([osi-cplex.pc:src/OsiCpx/osi-cplex.pc.in])
fi

if test "$coin_has_glpk" = yes ; then
  AC_CONFIG_FILES([osi-glpk.pc:src/OsiGlpk/osi-glpk.pc.in])
fi

if test "$coin_has_grb" = yes ; then
  AC_CONFIG_FILES([osi-gurobi.pc:src/OsiGrb/osi-gurobi.pc.in])
fi

if test "$coin_has_msk" = yes ; then
  AC_CONFIG_FILES([osi-mosek.pc:src/OsiMsk/osi-mosek.pc.in])
fi

if test "$coin_has_xpr" = yes ; then
  AC_CONFIG_FILES([osi-xpress.pc:src/OsiXpr/osi-xpress.pc.in])
fi

if test "$coin_has_soplex" = yes ; then
  AC_CONFIG_FILES([osi-soplex.pc:src/OsiSpx/osi-soplex.pc.in])
fi

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/Osi/config.h src/Osi/config_osi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
