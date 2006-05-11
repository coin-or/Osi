# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                             libOsiFmp                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiFmp.la

# List all source files for this library, including headers
libOsiFmp_la_SOURCES = \
	OsiFmpSolverInterface.cpp OsiFmpSolverInterface.hpp

# This is for libtool (on Windows)
libOsiFmp_la_LDFLAGS = -no-undefined

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(COINBASESRCDIR)/src` \
	-I`$(CYGPATH_W) $(FMPINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/..`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = OsiFmpSolverInterface.hpp
