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

IncDir :=
IncDir += ${zlibIncDir}
IncDir += ${bzlibIncDir}
IncDir += ${OsiSolverIncDir}
IncDir += ${OsiIncDir}
IncDir += ${OslIncDir}
IncDir += ${VolIncDir}
IncDir += ${CpxIncDir}
IncDir += ${SpxIncDir}
IncDir += ${XprIncDir}
IncDir += ${ClpIncDir}

LibDir :=
LibDir += ${zlibLibDir}
LibDir += ${bzlibLibDir}
LibDir += ${OsiSolverLibDir}
LibDir += ${OsiLibDir}
LibDir += ${OslLibDir}
LibDir += ${VolLibDir}
LibDir += ${CpxLibDir}
LibDir += ${SpxLibDir}
LibDir += ${XprLibDir}
LibDir += ${ClpLibDir}

LibName :=
LibName += ${zlibLibName}
LibName += ${bzlibLibName}
LibName += ${OsiSolverLibName}
LibName += ${OsiLibName}
LibName += ${OslLibName}
LibName += ${VolLibName}
LibName += ${CpxLibName}
LibName += ${SpxLibName}
LibName += ${XprLibName}
LibName += ${ClpLibName}

Define :=
Define += ${zlibDefine}
Define += ${bzlibDefine}
Define += ${OslDefine}
Define += ${VolDefine}
Define += ${CpxDefine}
Define += ${SpxDefine}
Define += ${XprDefine}
Define += ${ClpDefine}
Define += ${OsiDefine}

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

#########################################################################

CXXFLAGS += $(addprefix -I,${IncDir})
CXXFLAGS += $(addprefix -D,${Define})

LIBDIRS := ${LibDir}
LIBS    := ${LibName}

LDFLAGS := $(addprefix -L,$(LIBDIRS))
LDFLAGS += $(call ifprefix,$(SHLINKPREFIX),$(LIBDIRS))
LDFLAGS += $(patsubst lib%,-l%,$(basename $(LIBS)))

###############################################################################

TESTSRC := $(wildcard *Test.cpp)
TESTOBJ := $(addprefix $(TARGETDIR)/, $(TESTSRC:.cpp=.o))
TESTDEP := $(addprefix $(DEPDIR)/, $(TESTSRC:.cpp=.d))

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
