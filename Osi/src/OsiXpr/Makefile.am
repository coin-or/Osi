# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                             libOsiXpr                                #
########################################################################

# Name of the library compiled in this directory.
lib_LTLIBRARIES = libOsiXpr.la

# List all source files for this library, including headers
libOsiXpr_la_SOURCES = \
	OsiXprSolverInterface.cpp OsiXprSolverInterface.hpp

# List all additionally required libraries
libOsiXpr_la_LIBADD = ../Osi/libOsi.la $(OSIXPRLIB_LFLAGS)

# This is for libtool (on Windows)
libOsiXpr_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/../Osi -I$(srcdir)/../OsiCommonTest $(OSIXPRLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(pkgincludedir)
includecoin_HEADERS = OsiXprSolverInterface.hpp
