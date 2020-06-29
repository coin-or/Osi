# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.


include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiSpx                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiSpx.la

# List all source files for this library, including headers
libOsiSpx_la_SOURCES = OsiSpxSolverInterface.cpp OsiSpxSolverInterface.hpp

# List all additionally required libraries
libOsiSpx_la_LIBADD = ../Osi/libOsi.la $(OSISPXLIB_LFLAGS)

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSISPXLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiSpxSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osispx.h $(DESTDIR)$(includecoindir)/OsiSpxConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiSpxConfig.h
