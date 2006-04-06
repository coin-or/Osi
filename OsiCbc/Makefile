# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
OptLevel := -O
# OptLevel := -g

LIBNAME := OsiCbc

LIBSRC :=
LIBSRC += OsiCbcSolverInterface.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

export ExtraIncDir  := \
    ${CbcIncDir}  ${OsiIncDir}  ${CoinIncDir}  ${zlibIncDir}  ${bzlibIncDir} 
export ExtraLibDir  := \
    ${CbcLibDir}  ${OsiLibDir}  ${CoinLibDir}  ${zlibLibDir}  ${bzlibLibDir} 
export ExtraLibName := \
    ${CbcLibName} ${OsiLibName} ${CoinLibName} ${zlibLibName} ${bzlibLibName}
export ExtraDefine  := \
    ${CbcDefine}  ${OsiDefine}  ${CoinDefine}  ${zlibDefine}  ${bzlibDefine} 

# By default, OsiCbc uses OsiClp as the underlying OSI solver. You can change
# that here. If you want OsiClp, LEAVE THESE LINES COMMENTED OUT AND DO
# NOTHING ELSE. If you want a solver other than dylp, change as appropriate.
# Consult Makefiles/Makefile.location for your solver.

# ExtraDefine += OSICBC_DFLT_SOLVER=OsiDylpSolverInterface
# ExtraDefine += OSICBC_DFLT_SOLVER_HPP=\"OsiDylpSolverInterface.hpp\"
# ExtraIncDir += $(DylpIncDir) $(OsiDylpIncDir)
# ExtraLibDir += $(DylpLibDir) $(OsiDylpLibDir)
# ExtraLibName += $(DylpLibName) $(OsiDylpLibName)
# ExtraDefine += $(DylpDefine) $(OsiDylpDefine)

export LibType OptLevel LIBNAME LIBSRC

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libosi library clean doc

default: install

install library : % :
	(cd $(CoinDir)/Cbc && $(MAKE) -f Makefile.lightweight install)
	(cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

libOsiCbc :
	(cd $(CoinDir)/Osi && $(MAKE) libOsi)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

# Make a decent clean that works for just OsiCbc. Running clean up in Osi will
# work, but it also cleans out all other solvers and the main OSI.  We must
# remove the libraries and include files from from the central COIN
# directories, otherwise changes will not be seen.

clean : 
	@ rm -rf $(UNAME)*
	@ rm -rf dep
	@ rm -rf Doc
	@ (cd $(CoinLibDir) && rm -f libOsiCbc*)
	@ (cd $(CoinIncDir) && rm -f OsiCbc*.hpp)

# Make a doc target that works.

doc :
	@ doxygen $(MakefileDir)/doxygen.conf
