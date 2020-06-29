# Copyright (C) 2009 Stefan Vigerske and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.


include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiGrb                                 
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiGrb.la

# List all source files for this library, including headers
libOsiGrb_la_SOURCES = OsiGrbSolverInterface.cpp OsiGrbSolverInterface.hpp

# List all additionally required libraries
libOsiGrb_la_LIBADD = ../Osi/libOsi.la $(OSIGRBLIB_LFLAGS)

# This is for libtool (on Windows)
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi $(OSIGRBLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiGrbSolverInterface.hpp

install-exec-local:
	$(install_sh_DATA) config_osigrb.h $(DESTDIR)$(includecoindir)/OsiGrbConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiGrbConfig.h
