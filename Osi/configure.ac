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

# Initialize libtool.

AC_COIN_PROG_LIBTOOL

#############################################################################
#                   Math headers                                            #
#############################################################################

# Check for cmath/math.h, cfloat/float.h, cieeefp/ieeefp.h

AC_COIN_CHECK_MATH_HDRS

#############################################################################
#                              COIN-OR components                           #
#############################################################################

# Initialise variables for OsiLib

AC_COIN_INIT_LIBVARS([OsiLib],[-L\${libdir} -lOsi],[-I\${includedir}])

# Find components for OsiLib

AC_COIN_CHK_PKG(CoinUtils,[OsiLib])
if test "$coin_has_coinutils" != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi

# Finalize libs and cflags for building OsiLib.

AC_COIN_FINALIZE_FLAGS(OsiLib)

# Initialise variables for OsiCommonTestLib

AC_COIN_INIT_LIBVARS([OsiCommonTestLib],
		     [-L\${libdir} -lOsiCommonTest],[-I\${includedir}])

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

AC_COIN_INIT_LIBVARS([OsiGlpkLib],[-L\${libdir} -lOsiGlpk],[-I\${includedir}])
AC_COIN_CHK_PKG(Glpk,[OsiGlpkLib],,[build])
AC_COIN_INIT_LIBVARS([OsiSpxLib],[-L\${libdir} -lOsiSpx],[-I\${includedir}])
AC_COIN_CHK_PKG(SoPlex,[OsiSpxLib],[coinsoplex < 1.7],[build])
AC_COIN_FINALIZE_FLAGS([OsiGlpkLib OskSpxLib])
AC_COIN_CHK_HERE([OsiLib],[OsiGlpkLib OsiSpxLib],[osi])

# These are proprietary third-party solvers. Cplex detection has been
# tested on a semi-regular basis, and Mosek detection most likely works. For
# the rest, it's anyone's guess.  -- lh, 080529 --
# fixed detection of Mosek and Xpress -- stefan, 091003 (linux32, gcc) --

# AC_COIN_CHECK_USER_LIBRARY([Cplex], [CPX],[cplex.h],   [CPXgetstat])
# AC_COIN_CHECK_USER_LIBRARY([Mosek], [MSK],[mosek.h],   [MSK_makeenv])
# AC_COIN_CHECK_USER_LIBRARY([Xpress],[XPR],[xprs.h],    [XPRSinit])
# AC_COIN_CHECK_USER_LIBRARY([Gurobi],[GRB],[gurobi_c.h],[GRBloadenv])

AC_COIN_INIT_LIBVARS([OsiCpxLib],[-L\${libdir} -lOsiCpx],[-I\${includedir}])
AC_COIN_CHK_LIB(Cplex,[OsiCpxLib],[-lcplex],[-I/usr/include])
AC_COIN_INIT_LIBVARS([OsiMskLib],[-L\${libdir} -lOsiMsk],[-I\${includedir}])
AC_COIN_CHK_LIB(Mosek,[OsiMskLib],[-lmosek],[-I/usr/include])
AC_COIN_INIT_LIBVARS([OsiXprLib],[-L\${libdir} -lOsiXpr],[-I\${includedir}])
AC_COIN_CHK_LIB(Xpress,[OsiXprLib],[-lxpress],[-I/usr/include])
AC_COIN_INIT_LIBVARS([OsiGrbLib],[-L\${libdir} -lOsiGrb],[-I\${includedir}])
AC_COIN_CHK_LIB(Gurobi,[OsiGrbLib],[-lgurobi],[-I/usr/include])
AC_COIN_FINALIZE_FLAGS([OsiCpxLib OsiMskLib OsiXprLib OsiGrbLib])
AC_COIN_CHK_HERE([OsiLib],[OsiCpxLib OsiMskLib OsiXprLib OsiGrbLib],[osi])

# Force the automake conditionals to false until third-party link check is
# reimplemented.

AM_CONDITIONAL(COIN_HAS_CPX,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_MSK,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_XPR,[test 0 = 1])
AM_CONDITIONAL(COIN_HAS_GRB,[test 0 = 1])


#############################################################################
#                            Examples solver                                #
#############################################################################

# Choose a solver interface that can be used in the examples

if test "$coin_has_glpk" = yes ; then
  OSI_EXAMPLES_SOLVER_NAME=OsiGlpkSolverInterface
  OSI_EXAMPLES_SOLVER_CFLAGS="$OSIGLPKLIB_CFLAGS"
  OSI_EXAMPLES_SOLVER_LIBS="$OSIGLPKLIB_LFLAGS"
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

