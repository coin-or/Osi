# Copyright (C) 2006, 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                         unitTest for Coin                            #
########################################################################

noinst_PROGRAMS = unitTest

unitTest_SOURCES = \
	OsiColCutTest.cpp \
	OsiCutsTest.cpp \
	OsiRowCutDebuggerTest.cpp \
	OsiRowCutTest.cpp \
	OsiSolverInterfaceTest.cpp \
	unitTest.cpp

# List libraries of COIN projects
unitTest_LDADD =

# List additional libraries
#LIBS +=

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../src`

if COIN_HAS_COINUTILS_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(COINUTILSINSTDIR)/include/coin`
else
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Depending of what solvers are available, we add the corresponding files,
# libraries and include dirs

if COIN_HAS_CLP
  unitTest_SOURCES += OsiClpSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiClp`
  unitTest_LDADD += ../src/OsiClp/libOsiClp.la
if COIN_HAS_CLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CLPINSTDIR)/lib` -lClp \
	`cat $(CLPDOCDIR)/clp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc`
  unitTest_LDADD += $(CLPOBJDIR)/src/libClp.la
  LIBS  += `cat $(CLPOBJDIR)/clp_addlibs.txt`
endif
endif

if COIN_HAS_CBC
  unitTest_SOURCES += OsiCbcSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiCbc`
  unitTest_LDADD += ../src/OsiCbc/libOsiCbc.la
if COIN_HAS_CBC_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(CBCINSTDIR)/lib` -lCbc \
	`cat $(CBCDOCDIR)/cbc_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc`
  unitTest_LDADD += $(CBCOBJDIR)/src/libCbc.la
  LIBS += `cat $(CBCOBJDIR)/cbc_addlibs.txt`
endif
if COIN_HAS_ASL
  unitTest_LDADD += $(ASLLIB)
  AM_CPPFLAGS += $(ASL_CPPFLAGS) -DCOIN_HAS_ASL
endif
endif

if COIN_HAS_SYMPHONY
  unitTest_SOURCES += OsiSymSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiSym`
  unitTest_LDADD += ../src/OsiSym/libOsiSym.la
if COIN_HAS_SYMPHONY_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(SYMPHONYINSTDIR)/lib` -lSym \
	`cat $(SYMPHONYDOCDIR)/sym_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SYMPHONYSRCDIR)/include`
  unitTest_LDADD += $(SYMPHONYOBJDIR)/src/libSym.la
  LIBS += `cat $(SYMPHONYOBJDIR)/sym_addlibs.txt`
endif
endif

if COIN_HAS_CGL
if COIN_HAS_CGL_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(CGLINSTDIR)/lib` -lCgl \
	`cat $(CGLDOCDIR)/cgl_addlibs.txt`
else
  unitTest_LDADD += $(CGLOBJDIR)/src/libCgl.la
  LIBS += `cat $(CGLOBJDIR)/cgl_addlibs.txt`
endif
endif

if COIN_HAS_CPX
  unitTest_SOURCES += OsiCpxSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiCpx`
  unitTest_LDADD += ../src/OsiCpx/libOsiCpx.la
endif

if COIN_HAS_DYLP
  unitTest_SOURCES += OsiDylpSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiDylp`
  unitTest_LDADD += ../src/OsiDylp/libOsiDylp.la
if COIN_HAS_DYLP_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(DYLPINSTDIR)/lib` -lDylp \
	`cat $(DYLPDOCDIR)/dylp_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(DYLPSRCDIR)/src/DylpStdLib` \
	-I`$(CYGPATH_W) $(DYLPSRCDIR)/src/Dylp` \
	-I`$(CYGPATH_W) $(DYLPOBJDIR)/inc`
  unitTest_LDADD += $(DYLPOBJDIR)/src/Dylp/libDylp.la
  LIBS += `cat $(DYLPOBJDIR)/dylp_addlibs.txt`
endif
endif

if COIN_HAS_FMP
  unitTest_SOURCES += OsiFmpSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiFmp`
  unitTest_LDADD += ../src/OsiFmp/libOsiFmp.la
endif

if COIN_HAS_GLPK
  unitTest_SOURCES += OsiGlpkSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(GLPKINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiGlpk`
  unitTest_LDADD += ../src/OsiGlpk/libOsiGlpk.la
endif

if COIN_HAS_MSK
  unitTest_SOURCES += OsiMskSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(MSKINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiMsk`
  unitTest_LDADD += ../src/OsiMsk/libOsiMsk.la
endif

if COIN_HAS_OSL
  unitTest_SOURCES += OsiOslSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(OSLINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiOsl`
  unitTest_LDADD += ../src/OsiOsl/libOsiOsl.la
endif

if COIN_HAS_SPX
  unitTest_SOURCES += OsiSpxSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(SPXINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiSpx`
  unitTest_LDADD += ../src/OsiSpx/libOsiSpx.la
endif

if COIN_HAS_VOL
  unitTest_SOURCES += OsiVolSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../src/OsiVol`
  unitTest_LDADD += ../src/OsiVol/libOsiVol.la
if COIN_HAS_VOL_PREINSTALLED
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLINSTDIR)/include/coin`
  LIBS += -L`$(CYGPATH_W) $(VOLINSTDIR)/lib` -lVol \
	`cat $(VOLDOCDIR)/vol_addlibs.txt`
else
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(VOLSRCDIR)/src`
  unitTest_LDADD += $(VOLOBJDIR)/src/libVol.la
  LIBS += `cat $(VOLOBJDIR)/vol_addlibs.txt`
endif
endif

if COIN_HAS_XPR
  unitTest_SOURCES += OsiXprSolverInterfaceTest.cpp
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(XPRINCDIR)` \
	-I`$(CYGPATH_W) $(srcdir)/../src/OsiXpr`
  unitTest_LDADD += ../src/OsiXpr/libOsiXpr.la
endif

# Finally the Coin library (many of the above need it)
unitTest_LDADD += ../src/libOsi.la
LIBS += `cat $(OSIOBJDIR)/osi_addlibs.txt` $(ADDLIBS)
if COIN_HAS_COINUTILS_PREINSTALLED
  LIBS += -L`$(CYGPATH_W) $(COINUTILSINSTDIR)/lib` -lCoinUtils \
	`cat $(COINUTILSDOCDIR)/coinutils_addlibs.txt`
else
  unitTest_LDADD += $(COINUTILSOBJDIR)/src/libCoinUtils.la
  LIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`
endif

test: unitTest$(EXEEXT)
	@if test x"$(DIRNETLIB)" != x; then \
		NETLIBARG="-netlibDir $(DIRNETLIB) -testOsiSolverInterface"; \
	fi; \
	./unitTest$(EXEEXT) $$NETLIBARG

.PHONY: test

########################################################################
#                          Cleaning stuff                              #
########################################################################

# Here we list everything that is not generated by the compiler, e.g.,
# output files of a program

DISTCLEANFILES = *.mps *.mps.gz *.lp test2out *.out.gz *.out
