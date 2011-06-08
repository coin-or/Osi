# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

########################################################################
#                             libOsiSpx                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiSpx.la

# List all source files for this library, including headers
libOsiSpx_la_SOURCES = OsiSpxSolverInterface.cpp OsiSpxSolverInterface.hpp

# This is for libtool
libOsiSpx_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = $(COINUTILS_CFLAGS) $(SOPLEX_CFLAGS) \
	-I`$(CYGPATH_W) $(srcdir)/../Osi`

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Osi

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = OsiSpxSolverInterface.hpp
