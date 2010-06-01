# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                             libOsiCbc                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiCbc.la

# List all source files for this library, including headers
libOsiCbc_la_SOURCES = \
	OsiCbcSolverInterface.cpp OsiCbcSolverInterface.hpp

# This is for libtool (on Windows)
libOsiCbc_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.

# Needed by any OsiCbc configuration. Clp is inextricably wired into Cbc,
# hence OsiCbc always needs the Clp include directories. -- lh, 100601 --

AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(srcdir)/..`

# Plus the appropriate include directories for the default solver. The
# underlying configuration macro (Osi/m4/ac_osicbc_config.m4) is set up to
# handle any Osi solver. Individual OsiXxx maintainers can add their solver
# here when they've done the testing. Typically, you want to duplicate the
# include directories specified for AM_CPPFLAGS in OsiXXX/Makefile.am, keeping
# in mind the common directories above do not need to be specified again.
# The current set matches existing Osi interfaces (excluding OSL and Vol) as of
# January 2008. Most of these are untested!
# -- lh, January 2008 --

if OSICBC_DFLT_SOLVER_CLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
endif

if OSICBC_DFLT_SOLVER_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

if OSICBC_DFLT_SOLVER_DYLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiDylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
endif

if OSICBC_DFLT_SOLVER_FMP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiFmp` \
	-I`$(CYGPATH_W) $(FMPINCDIR)`
endif

if OSICBC_DFLT_SOLVER_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiGlpk` \
	-I`$(CYGPATH_W) $(GLPKINCDIR)`
endif

if OSICBC_DFLT_SOLVER_GLPK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiGlpk` \
	-I`$(CYGPATH_W) $(GLPKINCDIR)`
endif

if OSICBC_DFLT_SOLVER_MSK
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiMsk` \
	-I`$(CYGPATH_W) $(MSKINCDIR)`
endif

if OSICBC_DFLT_SOLVER_SPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiSpx` \
	-I`$(CYGPATH_W) $(SPXINCDIR)`
endif

if OSICBC_DFLT_SOLVER_SYM
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiSym` \
	-I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/include`
endif

if OSICBC_DFLT_SOLVER_XPR
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiXpr` \
	-I`$(CYGPATH_W) $(XPRINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	OsiCbcSolverInterface.hpp
