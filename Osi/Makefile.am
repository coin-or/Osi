# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

SUBDIRS = src/Osi src/OsiCommonTest
if COIN_HAS_CPX
  SUBDIRS += src/OsiCpx
endif
if COIN_HAS_GLPK
  SUBDIRS += src/OsiGlpk
endif
if COIN_HAS_MSK
  SUBDIRS += src/OsiMsk
endif
if COIN_HAS_XPR
  SUBDIRS += src/OsiXpr
endif
if COIN_HAS_GRB
  SUBDIRS += src/OsiGrb
endif
if COIN_HAS_SPX
  SUBDIRS += src/OsiSpx
endif

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

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

# Doxygen documentation

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )


clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

distclean-local:
	if test -r test/Makefile; then cd test; $(MAKE) distclean; fi

install-exec-local: install-doc

uninstall-local: uninstall-doc

.PHONY: test unitTest doxydoc

########################################################################
#                  Installation of the addlibs file                    #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = osi.pc osi-unittests.pc

if COIN_HAS_PKGCONFIG
addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Osi
	
install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs osi > $(addlibsdir)/osi_addlibs.txt

uninstall-hook:
	rm -f $(addlibsdir)/osi_addlibs.txt
endif

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =

include BuildTools/Makemain.inc
