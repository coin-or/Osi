# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 1465 2010-02-07 16:36:26Z stefan $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                    Common Test library for Osi                       #
########################################################################

# Name of the library compiled in this directory.
# We want it to be installed in the 'lib/coin' directory

libcoindir = $(libdir)/coin
libcoin_LTLIBRARIES = libOsiCommonTests.la

libOsiCommonTests_la_SOURCES = OsiCommonTests.hpp \
	OsiColCutTest.cpp \
	OsiCutsTest.cpp \
	OsiRowCutDebuggerTest.cpp \
	OsiRowCutTest.cpp \
	OsiSimplexAPITest.cpp \
	OsiNetlibTest.cpp \
	OsiUnitTestUtils.cpp \
	OsiSolverInterfaceTest.cpp

includecoindir = $(includedir)/coin
includecoin_HEADERS = OsiUnitTests.hpp

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows

AM_CPPFLAGS = $(COINUTILS_CFLAGS) \
	-I`$(CYGPATH_W) $(srcdir)/../Osi`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

