# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                            libOsiGlpk                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiGlpk.la

# List all source files for this library, including headers
libOsiGlpk_la_SOURCES = OsiGlpkSolverInterface.cpp OsiGlpkSolverInterface.hpp

# List all additionally required libraries
libOsiGlpk_la_LIBADD = ../Osi/libOsi.la $(OSIGLPKLIB_LFLAGS)

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSIGLPKLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiGlpkSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osiglpk.h $(DESTDIR)$(includecoindir)/OsiGlpkConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiGlpkConfig.h
