# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Lou Hafer	SFU	2010-11-20

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                    Common Test library for Osi                       #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in $libdir
lib_LTLIBRARIES = libOsiCommonTest.la

# List all source files for this library, including headers

libOsiCommonTest_la_SOURCES = \
    OsiCommonTest.hpp \
	OsiColCutTest.cpp \
	OsiCutsTest.cpp \
	OsiRowCutDebuggerTest.cpp \
	OsiRowCutTest.cpp \
	OsiSimplexAPITest.cpp \
	OsiNetlibTest.cpp \
	OsiUnitTestUtils.cpp \
	OsiSolverInterfaceTest.cpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libOsiCommonTest_la_LIBADD = $(OSICOMMONTESTLIB_LFLAGS)
endif

# Libtool flags
libOsiCommonTest_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory. This
# "cygpath" stuff is necessary to compile with native compilers on Cygwin.

AM_CPPFLAGS = `$(CYGPATH_W) -I$(srcdir)/../Osi $(OSICOMMONTESTLIB_CFLAGS)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,

includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiUnitTests.hpp
