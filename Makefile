LIBS =
LIBS += osi
LIBS += osiosl
LIBS += osivol
#LIBS += osicpx
#LIBS += osispx
#LIBS += osixpr
#LIBS += osidylp

.DELETE_ON_ERROR:

.PHONY: 

default: all-lib

install: all-instlib

clean: all-cleanlib

all-lib: $(addprefix lib,$(LIBS))
all-instlib: $(addprefix instlib,$(LIBS))
all-cleanlib: $(addprefix cleanlib,$(LIBS))

lib% :
	${MAKE} -f Makefile.$* library

instlib% :
	${MAKE} -f Makefile.$* install

cleanlib% :
	${MAKE} -f Makefile.$* clean

