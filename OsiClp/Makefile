# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -g

LIBNAME := OsiClp

LIBSRC :=
LIBSRC += OsiClpSolverInterface.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export MakefileDir := $(CoinDir)/Common/make
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

OsiIncDir := $(CoinDir)/include
OsiLibDir := $(CoinDir)/lib

export ExtraIncDir := ${ClpIncDir} ${OsiIncDir}
export ExtraLibDir := ${ClpLibDir} ${OsiLibDir}
export ExtraLibName := ${ClpLibName} ${OsiLibName}
export ExtraDefine := ${ClpDefine}

export LibType OptLevel LIBNAME LIBSRC

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libosi library clean doc

default: install

install library : % :
	(cd $(CoinDir)/Clp && $(MAKE) install)
	(cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

libOsiClp :
	(cd $(CoinDir)/Osi && $(MAKE) libOsi)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

clean doc : % :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*
