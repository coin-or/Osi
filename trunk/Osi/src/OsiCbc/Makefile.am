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

# Needed by any OsiCbc configuration.

AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(srcdir)/..`

# Plus the appropriate include directories for the default solver. The
# underlying configuration macro (Osi/m4/ac_osicbc_config.m4) is set up to
# handle any Osi solver, but I'll leave it to individual OsiXxx maintainers to
# add their solver here when they've done the testing.  -- lh, July 2006 --

if OSICBC_DFLT_SOLVER_CLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiClp` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
endif

if OSICBC_DFLT_SOLVER_DYLP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../OsiDylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = OsiCbcSolverInterface.hpp
