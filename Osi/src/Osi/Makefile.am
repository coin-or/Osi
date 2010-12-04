# Copyright (C) 2010 Lou Hafer
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Lou Hafer        SFU        2010-07-29

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libOsi                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib/coin' directory

libcoindir = $(libdir)/coin
libcoin_LTLIBRARIES = libOsi.la

# List all source files for this library, including headers
libOsi_la_SOURCES = \
	OsiConfig.h \
	OsiAuxInfo.cpp OsiAuxInfo.hpp \
	OsiBranchingObject.cpp OsiBranchingObject.hpp \
	OsiChooseVariable.cpp OsiChooseVariable.hpp \
	OsiColCut.cpp OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.cpp OsiCut.hpp \
	OsiCuts.cpp OsiCuts.hpp \
	OsiNames.cpp \
	OsiPresolve.cpp OsiPresolve.hpp \
	OsiRowCut.cpp OsiRowCut.hpp \
	OsiRowCutDebugger.cpp OsiRowCutDebugger.hpp \
	OsiSolverBranch.cpp OsiSolverBranch.hpp \
	OsiSolverInterface.cpp OsiSolverInterface.hpp \
	OsiSolverParameters.hpp

# This is for libtool
libOsi_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags.
AM_CPPFLAGS = $(COINUTILS_CFLAGS)

# This line is necessary to allow VPATH compilation with MS compilers on
# Cygwin.  This "cygpath" stuff is necessary to compile with native compilers
# on Windows.

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'includedir/coin'

includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	../../inc/config_osi.h \
	OsiAuxInfo.hpp \
	OsiBranchingObject.hpp \
	OsiChooseVariable.hpp \
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.hpp \
	OsiCuts.hpp \
	OsiPresolve.hpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.hpp \
	OsiSolverBranch.hpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = OsiConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)

