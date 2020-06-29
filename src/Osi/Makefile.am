# Copyright (C) 2010 Lou Hafer
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Lou Hafer        SFU        2010-07-29

include $(top_srcdir)/BuildTools/Makemain.inc

########################################################################
#                               libOsi                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in $libdir
lib_LTLIBRARIES = libOsi.la

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
	OsiSolverParameters.hpp \
	OsiFeatures.cpp OsiFeatures.hpp

# List all additionally required libraries
libOsi_la_LIBADD = $(OSILIB_LFLAGS)

# This is for libtool
AM_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags.
AM_CPPFLAGS = $(OSILIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in pkgincludedir
includecoindir = $(pkgincludedir)
includecoin_HEADERS = \
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
	OsiSolverParameters.hpp \
	OsiFeatures.hpp

install-exec-local:
	$(install_sh_DATA) config_osi.h $(DESTDIR)$(includecoindir)/OsiConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/OsiConfig.h
