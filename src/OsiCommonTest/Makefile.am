# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Lou Hafer  SFU   2010-11-20

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
libOsiCommonTest_la_LIBADD = $(OSICOMMONTESTLIB_LFLAGS) ../Osi/libOsi.la

# Libtool flags
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSICOMMONTESTLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiUnitTests.hpp

install-exec-local:
	$(install_sh_DATA) config_osicommontest.h $(DESTDIR)$(includecoindir)/OsiUnitTestsConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiUnitTestsConfig.h
