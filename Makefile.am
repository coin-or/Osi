# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include BuildTools/Makemain.inc

########################################################################
#            Subdirectories and installation of .pc files              #
########################################################################

SUBDIRS = src/Osi src/OsiCommonTest

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

########################################################################
#                           Extra Targets                              #
########################################################################

test: all
	cd test; $(MAKE) test

unitTest: test

clean-local: clean-doxygen-docs
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc uninstall-doxygen-docs

.PHONY: test unitTest
