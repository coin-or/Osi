# Comment in the libraries you want to build besides libOsi.

SOLVERLIBS :=
SOLVERLIBS += libOsiOsl
#SOLVERLIBS += libOsiVol
SOLVERLIBS += libOsiClp
# SOLVERLIBS += libOsiCpx
# SOLVERLIBS += libOsiSpx
# SOLVERLIBS += libOsiXpr
# SOLVERLIBS += libOsiDylp

# Look at and if necessary edit the following files:
# - ../Common/make/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export SOLVERLIBS

.DELETE_ON_ERROR:

export CoinDir = $(shell cd ..; pwd)

.PHONY: default install clean unitTest
.PHONY: inst-libOsi $(addprefix inst-,$(SOLVERLIBS))
.PHONY: clean-libOsi $(addprefix clean-,$(SOLVERLIBS))

default: install

install: inst-libOsi $(addprefix inst-,$(SOLVERLIBS))

clean: clean-libOsi $(addprefix clean-,$(SOLVERLIBS))

###############################################################################

unitTest : install
	(cd Test && ${MAKE} unitTest)

libOsi : 
	${MAKE} -f Makefile.Osi library

$(SOLVERLIBS) : lib% :
	(cd $* && ${MAKE} library)

inst-libOsi :
	${MAKE} -f Makefile.Osi install

$(addprefix inst-,$(SOLVERLIBS)) : inst-lib% :
	(cd $* && ${MAKE} install)

clean-libOsi :
	${MAKE} -f Makefile.Osi clean

$(addprefix clean-,$(SOLVERLIBS)) : clean-lib% :
	(cd $* && ${MAKE} clean)
