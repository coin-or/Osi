SOLVERLIBS += libOsiOsl
SOLVERLIBS += libOsiVol
#LIBS += osicpx
#LIBS += osispx
#LIBS += osixpr
#LIBS += osidylp

.DELETE_ON_ERROR:

export CoinDir = $(shell cd ..; pwd)

.PHONY: 

default: install

install: inst-libOsi $(addprefix inst-,$(SOLVERLIBS))

clean: clean-libOsi $(addprefix clean-,$(SOLVERLIBS))

###############################################################################

libOsi : 
	${MAKE} -f Makefile.Osi library

$(SOLVERLIBS) : lib% :
	(cd $* && ${MAKE} -f Makefile.$* library)

inst-libOsi :
	${MAKE} -f Makefile.Osi install

$(addprefix inst-,$(SOLVERLIBS)) : inst-lib% :
	(cd $* && ${MAKE} -f Makefile.$* install)

clean-libOsi :
	${MAKE} -f Makefile.Osi clean

$(addprefix clean-,$(SOLVERLIBS)) : clean-lib% :
	(cd $* && ${MAKE} -f Makefile.$* clean)
