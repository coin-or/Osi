# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libOsi                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
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
	OsiOpbdpSolve.cpp OsiOpbdpSolve.hpp \
	OsiPresolve.cpp OsiPresolve.hpp \
	OsiRowCut.cpp OsiRowCut.hpp \
	OsiRowCutDebugger.cpp OsiRowCutDebugger.hpp \
	OsiSolverBranch.cpp OsiSolverBranch.hpp \
	OsiSolverInterface.cpp OsiSolverInterface.hpp \
	OsiSolverParameters.hpp

# This is for libtool (on Windows)
libOsi_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Now we check into which directories to recurse, and which subdirectory
# libraries to include
SUBDIRS = .
if COIN_HAS_CBC
  SUBDIRS += OsiCbc
endif
if COIN_HAS_CLP
  SUBDIRS += OsiClp
endif
if COIN_HAS_CPX
  SUBDIRS += OsiCpx
endif
if COIN_HAS_DYLP
  SUBDIRS += OsiDylp
endif
if COIN_HAS_FMP
  SUBDIRS += OsiFmp
endif
if COIN_HAS_GLPK
  SUBDIRS += OsiGlpk
endif
if COIN_HAS_MSK
  SUBDIRS += OsiMsk
endif
if COIN_HAS_OSL
  SUBDIRS += OsiOsl
endif
if COIN_HAS_SPX
  SUBDIRS += OsiSpx
endif
if COIN_HAS_SYMPHONY
  SUBDIRS += OsiSym
endif
if COIN_HAS_VOL
  SUBDIRS += OsiVol
endif
if COIN_HAS_XPR
  SUBDIRS += OsiXpr
endif

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	../inc/config_osi.h \
	OsiAuxInfo.hpp \
	OsiBranchingObject.hpp \
	OsiChooseVariable.hpp \
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.hpp \
	OsiCuts.hpp \
	OsiOpbdpSolve.hpp \
	OsiPresolve.hpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.hpp

#	OsiSolverBranch.hpp \
#OsiAuxInfo.hpp
#OsiOpbdpSolve.hpp

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
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
