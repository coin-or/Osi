# Look at and if necessary edit the following files:
# - ../Makefiles/Makefile.location
# - Makefile.Osi 
# - Osi*/Makefile for the libs you have specified above

###############################################################################

export CoinDir := $(shell cd ../..; pwd)
export MakefileDir := $(CoinDir)/Makefiles
include ${MakefileDir}/Makefile.coin
include ${MakefileDir}/Makefile.location

###############################################################################

.DELETE_ON_ERROR:

.PHONY: default clean unitTest

default: unitTest

unitTest clean :
	(cd .. && ${MAKE} install)
	${MAKE} -f Makefile.test $*
