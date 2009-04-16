# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Osi],[0.100.0],[coin-osi@list.coin-or.org])

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

AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(DyLP)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(SYMPHONY)
AC_COIN_HAS_PROJECT(Vol)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available. Cplex detection has been
# tested on a semi-regular basis, and Mosek detection most likely works. For
# the rest, it's anyone's guess.  -- lh, 080529 --

AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

# Glpk can be supplied as a ThirdParty code, or as an external library in the
# same was as the previous solvers.

AC_COIN_HAS_GLPK(Glpk)

# We need Ampl in case we're using Cbc and it's build with Ampl.

AC_COIN_HAS_ASL


#############################################################################
#		  Capability Tests					    #
#############################################################################

#############################################################################
#                 Configuration options for individual OSIs                 #
#############################################################################

# If any of the tests performed in this section actually require the presence
# of the solver (file presence, link checks, etc.) be sure to guard the call.

AC_OSICBC_CONFIG(clp)
AC_ODSI_INFO(yes)
AC_ODSI_STATISTICS(no)
AC_ODSI_PARANOIA(1)

# If glpk is available, choose the MIP solver. lpx_intopt is preferred;
# lpx_integer will be used as a fallback. coin_has_glpk and glpk_libcheck
# are defined and set by COIN_HAS_GLPK. If we're building ThirdParty/Glpk,
# we already know the answer. Make it look like we're checking.

if test x"$use_thirdpartyglpk" = xbuild; then
  AC_MSG_CHECKING([for glpk mip solver ])
  AC_DEFINE(GLPK_HAS_INTOPT,[1],
	    [Define to 1 if GLPK has the advanced B&B solver lpx_intopt])
  AC_MSG_RESULT([lpx_intopt])
else
  if test "x$coin_has_glpk" = "xtrue" && test "x$glpk_libcheck" = "xyes" ; then
    AC_OGSI_MIPSOLVER()
  fi
fi

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/OsiCbc/Makefile
                 src/OsiClp/Makefile
                 src/OsiCpx/Makefile
                 src/OsiDylp/Makefile
                 src/OsiFmp/Makefile
                 src/OsiGlpk/Makefile
                 src/OsiMsk/Makefile
                 src/OsiOsl/Makefile
                 src/OsiSpx/Makefile
                 src/OsiSym/Makefile
                 src/OsiVol/Makefile
                 src/OsiXpr/Makefile
                 test/Makefile
		 osi_addlibs.txt])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([inc/config_osi.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
