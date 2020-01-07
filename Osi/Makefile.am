# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

include BuildTools/Makemain.inc

########################################################################
#            Subdirectories and installation of .pc files              #
########################################################################

SUBDIRS = src/Osi src/OsiCommonTest

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = osi.pc osi-unittests.pc

if COIN_HAS_CPLEX
  SUBDIRS += src/OsiCpx
  pkgconfiglib_DATA += osi-cplex.pc
endif
if COIN_HAS_GLPK
  SUBDIRS += src/OsiGlpk
  pkgconfiglib_DATA += osi-glpk.pc
endif
if COIN_HAS_MOSEK
  SUBDIRS += src/OsiMsk
  pkgconfiglib_DATA += osi-mosek.pc
endif
if COIN_HAS_XPRESS
  SUBDIRS += src/OsiXpr
  pkgconfiglib_DATA += osi-xpress.pc
endif
if COIN_HAS_GUROBI
  SUBDIRS += src/OsiGrb
  pkgconfiglib_DATA += osi-gurobi.pc
endif
if COIN_HAS_SOPLEX
  SUBDIRS += src/OsiSpx
  pkgconfiglib_DATA += osi-soplex.pc
endif

# We don't want to compile the test subdirectory, unless the test target
# is specified.  But we need to list it as subdirectory to make sure that
# it's included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	examples/basic2.cpp \
	examples/basic.cpp \
	examples/build.cpp \
	examples/Makefile.in \
	examples/parameters.cpp \
	examples/query.cpp \
	examples/README \
	examples/specific.cpp

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean

DISTCLEANFILES =

