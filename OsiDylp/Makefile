# Static or shared libraries should be built (STATIC or SHARED)?
LibType := SHARED

# Select optimization (-O or -g). -O will be automatically bumped up to the 
# highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1 or -O2
#OptLevel := -g
OptLevel := -O

LIBNAME := OsiDylp

LIBSRC :=
LIBSRC += OsiDylpMessages.cpp \
	  OsiDylpPresolve.cpp \
	  OsiDylpSolverInterface.cpp \
	  OsiDylpWarmStartBasis.cpp

##############################################################################
# You should not need to edit below this line.
##############################################################################
# The location of the customized Makefiles
export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

export ExtraIncDir  := \
    ${DylpIncDir}  ${OsiIncDir}  ${CoinIncDir}  ${zlibIncDir}  ${bzlibIncDir} 
export ExtraLibDir  := \
    ${DylpLibDir}  ${OsiLibDir}  ${CoinLibDir}  ${zlibLibDir}  ${bzlibLibDir} 
export ExtraLibName := \
    ${DylpLibName} ${OsiLibName} ${CoinLibName} ${zlibLibName} ${bzlibLibName}
export ExtraDefine  := \
    ${DylpDefine}  ${OsiDefine}  ${CoinDefine}  ${zlibDefine}  ${bzlibDefine} 

export LibType OptLevel LIBNAME LIBSRC

# This is normally defined in Makefile.lib. To make it easier to add extra
# actions for install and clean, I've redefined it here. The two definitions
# should match, obviously. After a bit of experience, perhaps some of the extra
# actions can be moved to Makefile.lib.

InstallDir := $(CoinDir)

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install libosi library clean doc

# Define this to get an mps dump of problem when dylp fails

# CXXFLAGS += -DDYLP_POSTMORTEM

default: install

library :
	(cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)
	$(MAKE) -f ${MakefileDir}/Makefile.lib library 

# The general sequence here is to make sure the main OSI library is installed,
# then install OsiDylp. Sccs likes to make files non-writeable unless they're
# actually being edited. I need to make the copies in $(InstallDir)/include
# writeable before Makefile.lib attempts to delete them.

install :
	(cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)
	@(cd $(InstallDir)/include && \
	  ([ -f OsiDylpSolverInterface.hpp ] && chmod a+w OsiDylp*.hpp) || \
	true)
	$(MAKE) -f ${MakefileDir}/Makefile.lib install

libOsiDylp :
	(cd $(CoinDir)/Osi && $(MAKE) libOsi)
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*

# Make a decent clean that works for just OsiDylp. Running clean up in Osi
# will work, but it also cleans out all other solvers and the main OSI.

clean : 
	$(MAKE) -f ${MakefileDir}/Makefile.lib clean
	rm -rf $(UNAME)*
	rm -rf dep
	(cd $(InstallDir)/lib && rm libOsiDylp*)

doc :
	$(MAKE) -f ${MakefileDir}/Makefile.lib $*
