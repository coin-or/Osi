# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiCpx                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiCpx.la

# List all source files for this library, including headers
libOsiCpx_la_SOURCES = OsiCpxSolverInterface.cpp OsiCpxSolverInterface.hpp

# List all additionally required libraries
libOsiCpx_la_LIBADD = ../Osi/libOsi.la $(OSICPXLIB_LFLAGS)

# This is for libtool (on Windows)
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSICPXLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiCpxSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osicpx.h $(DESTDIR)$(includecoindir)/OsiCpxConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiCpxConfig.h
