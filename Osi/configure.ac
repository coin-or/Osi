# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: configure.ac,v 1.1.2.1 2006/04/19 23:25:12 andreasw Exp $

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Osi],[0.1],[coin-osi@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

AC_REVISION(0.9)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/OsiAuxInfo.cpp)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_SRCDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE

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
AC_COIN_HAS_PROJECT(Sym)
AC_COIN_HAS_PROJECT(Vol)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosekdl.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

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
AC_OUTPUT

# ...and tell the user, that the configuration finished
AC_MSG_NOTICE([Osi configuration successful])
