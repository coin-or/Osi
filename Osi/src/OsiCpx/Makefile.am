# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiCpx                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiCpx.la

# List all source files for this library, including headers
libOsiCpx_la_SOURCES = \
	OsiCpxSolverInterface.cpp OsiCpxSolverInterface.hpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libOsiCpx_la_LIBADD = $(OSILIB_LFLAGS) ../Osi/libOsi.la $(OSICPXLIB_LFLAGS)
endif

# This is for libtool (on Windows)
libOsiCpx_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.

AM_CPPFLAGS = \
	`$(CYGPATH_W) $(OSICPXLIB_CFLAGS)` \
	`$(CYGPATH_W) -I$(srcdir)/../Osi`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiCpxSolverInterface.hpp
