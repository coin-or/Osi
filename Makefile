# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export CoinDir := $(shell cd ..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

# Comment in the libraries you want to build besides libOsi.

SOLVERLIBS :=
ifneq ($(filter COIN_libVol,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiVol
endif
ifneq ($(filter COIN_libOsl,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiOsl
endif
ifneq ($(filter COIN_libClp,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiClp
endif
ifneq ($(filter COIN_libCpx,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiCpx
endif
ifneq ($(filter COIN_libSpx,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiSpx
endif
ifneq ($(filter COIN_libXpr,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiXpr
endif
ifneq ($(filter COIN_libDylp,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiDylp
endif
ifneq ($(filter COIN_libGlpk,$(CoinLibsDefined)),)
    SOLVERLIBS += libOsiGlpk
endif

export SOLVERLIBS

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install clean doc unitTest
.PHONY: inst-libOsi $(addprefix inst-,$(SOLVERLIBS))
.PHONY: clean-libOsi $(addprefix clean-,$(SOLVERLIBS))
.PHONY: doc-libOsi $(addprefix doc-,$(SOLVERLIBS))

default: install

install: inst-libOsi $(addprefix inst-,$(SOLVERLIBS))

clean: clean-libOsi $(addprefix clean-,$(SOLVERLIBS))

doc:
	doxygen $(MakefileDir)/doxygen.conf

###############################################################################

unitTest : install
	(cd Test && ${MAKE} unitTest)

libOsi : 
	(cd $(CoinDir)/Coin && $(MAKE))
	${MAKE} -f Makefile.Osi library

inst-libOsi : libOsi
	${MAKE} -f Makefile.Osi install

$(addprefix inst-,$(SOLVERLIBS)) : inst-lib% :
	(cd $* && ${MAKE} install)

clean-libOsi :
	@rm -rf Junk
	@rm -rf $(UNAME)
	@rm -rf dep

$(addprefix clean-,$(SOLVERLIBS)) : clean-lib% :
	(cd $* && ${MAKE} clean)
