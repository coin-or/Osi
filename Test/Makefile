# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export MakefileDir := $(CoinDir)/Common/make
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

IncDir  := ${OsiIncDir}
LibDir  := ${OsiLibDir}
LibName := ${OsiLibName}
Define  :=
TESTSRC :=

ifeq ($(findstring libOsiClp,${SOLVERLIBS}),libOsiClp)
    IncDir += ${ClpIncDir}
    LibDir += ${ClpLibDir}
    LibName += ${ClpLibName} libOsiClp.so
    Define += COIN_USE_CLP
    TESTSRC += OsiClpSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiCpx,${SOLVERLIBS}),libOsiCpx)
    IncDir += ${CpxIncDir}
    LibDir += ${CpxLibDir}
    LibName += ${CpxLibName} libOsiCpx.so
    Define += COIN_USE_CPX
    TESTSRC += OsiCpxSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiDylp,${SOLVERLIBS}),libOsiDylp)
    IncDir += ${DylpIncDir}
    LibDir += ${DylpLibDir}
    LibName += ${DylpLibName} libOsiDylp.so
    Define += COIN_USE_DYLP
    TESTSRC += OsiDylpSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiOsl,${SOLVERLIBS}),libOsiOsl)
    IncDir += ${OslIncDir}
    LibDir += ${OslLibDir}
    LibName += ${OslLibName} libOsiOsl.so
    Define += COIN_USE_OSL
    TESTSRC += OsiOslSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiSpx,${SOLVERLIBS}),libOsiSpx)
    IncDir += ${SpxIncDir}
    LibDir += ${SpxLibDir}
    LibName += ${SpxLibName} libOsiSpx.so
    Define += COIN_USE_SPX
    TESTSRC += OsiSpxSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiVol,${SOLVERLIBS}),libOsiVol)
    IncDir += ${VolIncDir}
    LibDir += ${VolLibDir}
    LibName += ${VolLibName} libOsiVol.so
    Define += COIN_USE_VOL
    TESTSRC += OsiVolSolverInterfaceTest.cpp
endif

ifeq ($(findstring libOsiXpr,${SOLVERLIBS}),libOsiXpr)
    IncDir += ${XprIncDir}
    LibDir += ${XprLibDir}
    LibName += ${XprLibName} libOsiXpr.so
    Define += COIN_USE_XPR
    TESTSRC += OsiXprSolverInterfaceTest.cpp
endif

##############################################################################

TESTSRC += OsiColCutTest.cpp
TESTSRC += OsiCutsTest.cpp
TESTSRC += OsiIndexedVectorTest.cpp
TESTSRC += OsiMpsReaderTest.cpp
TESTSRC += OsiPackedMatrixTest.cpp
TESTSRC += OsiPackedVectorTest.cpp
TESTSRC += OsiRowCutDebuggerTest.cpp
TESTSRC += OsiRowCutTest.cpp
TESTSRC += OsiShallowPackedVectorTest.cpp
TESTSRC += OsiTripleCompareTest.cpp
TESTSRC += OsiSolverInterfaceTest.cpp
TESTSRC += unitTest.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################

CXXFLAGS += $(OPTFLAG)
ifeq ($(OptLevel),-g)
    CXXFLAGS += -DOSI_DEBUG
endif

###############################################################################

space:= $(empty) $(empty)
OptVersion := $(subst $(space),_,$(OptLevel))

TARGETDIR := $(UNAME)$(OptVersion)
DEPDIR := dep

VPATH := . : include : Junk : ${TARGETDIR} : ${DEPDIR}

TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

#########################################################################

CXXFLAGS += $(addprefix -I,${IncDir})
CXXFLAGS += $(addprefix -D,${Define})

LIBDIRS := ${LibDir}
LIBS    := ${LibName}

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(call ifprefix,$(SHLINKPREFIX),$(LIBDIRS))
LDFLAGS += $(patsubst lib%,-l%,$(basename $(LIBS)))

###############################################################################
# This one must be right before the rules

include ${MakefileDir}/Makefile.rules

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean doc unitTest libdepend

default: unitTest

###############################################################################

unitTest: $(TARGETDIR)/unitTest

libdepend :
	(cd .. && ${MAKE} install)

$(TARGETDIR)/unitTest : libdepend $(TESTOBJ)
	@rm -rf Junk
	@echo ""
	@echo Creating unitTest
	@mkdir -p $(TARGETDIR)
	@rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(TESTOBJ) $(LDFLAGS) $(SYSLD) -lm
	${CP} $@ ..

###############################################################################

doc:
	doxygen $(MakefileDir)/doxygen.conf

clean :
	rm -rf Junk
	@rm -rf $(DEPDIR)
	@rm -rf $(TARGETDIR)

###############################################################################

%::
	@mkdir -p Junk
	touch Junk/$(notdir $@)

###############################################################################

-include $(TESTDEP)
