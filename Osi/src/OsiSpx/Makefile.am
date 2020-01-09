# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

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
libOsiSpx_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi -I$(srcdir)/../OsiCommonTest $(OSISPXLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiSpxSolverInterface.hpp
