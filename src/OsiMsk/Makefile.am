# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiMsk                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiMsk.la

# List all source files for this library, including headers
libOsiMsk_la_SOURCES = OsiMskSolverInterface.cpp OsiMskSolverInterface.hpp

# List all additionally required libraries
libOsiMsk_la_LIBADD = ../Osi/libOsi.la $(OSIMSKLIB_LFLAGS)

# This is for libtool (on Windows)
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSIMSKLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiMskSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osimsk.h $(DESTDIR)$(includecoindir)/OsiMskConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiMskConfig.h
