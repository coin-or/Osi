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
	OsiAuxInfo.cpp OsiAuxInfo.hpp \
	OsiColCut.cpp OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.cpp OsiCut.hpp \
	OsiCuts.cpp OsiCuts.hpp \
	OsiOpbdpSolve.cpp OsiOpbdpSolve.hpp \
	OsiPresolve.cpp OsiPresolve.hpp \
	OsiRowCut.cpp OsiRowCut.hpp \
	OsiRowCutDebugger.cpp OsiRowCutDebugger.hpp \
	OsiSolverBranch.cpp OsiSolverBranch.hpp \
	OsiSolverInterface.cpp OsiSolverInterface.hpp \
	OsiSolverParameters.hpp

# This is for libtool (on Windows)
libOsi_la_LDFLAGS = # -no-undefined

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
# "top_srcdir" refers to the basic directory for the main package that is
# being compiled.
CXXFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src`

# Now we check into which directories to recurse, and which subdirectory
# libraries to include
SUBDIRS =
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
if COIN_HAS_SYM
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
	OsiColCut.hpp \
	OsiCollections.hpp \
	OsiCut.hpp \
	OsiCuts.hpp \
	OsiPresolve.hpp \
	OsiSolverInterface.hpp \
	OsiSolverParameters.hpp \
	OsiRowCut.hpp \
	OsiRowCutDebugger.hpp

#	OsiSolverBranch.hpp \
#OsiAuxInfo.hpp
#OsiOpbdpSolve.hpp

########################################################################
#     Generating the configure header that is to be distributed        #
########################################################################

# You should only need to adapt the following two lines
CONFIG_DIST = OsiConfig.hpp
CONFIG_H = config_osi.h
$(libOsi_la_OBJECTS): ../inc/$(CONFIG_DIST)

nodist_include_HEADERS = ../inc/$(CONFIG_DIST)


../inc/$(CONFIG_DIST): ../inc/$(CONFIG_H)
	echo "/** Required defines from $(CONFIG_H) */" >bla.tmp
	echo "#ifndef __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	echo "#define __$(CONFIG_DIST)__" | sed -e 's/\./_/g' >>bla.tmp
	$(EGREP) 'COIN_' $< >>bla.tmp
	echo "#endif" >>bla.tmp
	if test -r ../inc/$(CONFIG_DIST); then \
	  if diff ../inc/$(CONFIG_DIST) bla.tmp >/dev/null 2>&1 ; then \
	    rm -f bla.tmp; \
	  else \
	    rm -f ../inc/$(CONFIG_DIST); \
	    mv bla.tmp ../inc/$(CONFIG_DIST); \
	  fi; \
	else \
	  mv bla.tmp ../inc/$(CONFIG_DIST); \
	fi

CLEANFILES = bla.tmp

DISTCLEANFILES = ../inc/$(CONFIG_DIST)
