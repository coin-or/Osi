# This is the makefile for OsiDylp, the OSI layer for the dylp LP solver.

# cvs/svn: $Id$

# Specify a shared (SHARED) or static (STATIC) library. Shared is the preferred
# choice unless you have some reason to want a static library.

export LibType ?= SHARED

# Ask the compiler to include debugging information (-g) and/or optimisation
# (-O). These are not mutually exclusive, but be sure you understand what
# you're asking for if you ask for both. -O will be automatically bumped up to
# the highest level of optimization the compiler supports. If want something in
# between then specify the exact level you want, e.g., -O1, -O2, etc.

export OptLevel ?= -O

# If OsiDylp is being built as a dependency for the Dylp sample main programs,
# this symbol will be defined. If it's not defined, the two options that make
# a difference for OsiDylp are stats and info. Define info if you want
# information prints, and stats if you want statistics. Note that both OsiDylp
# and dylp must be built with statistics enabled. dylp will ignore an incoming
# statistics structure if built without statistics, and OsiDylp will not supply
# one if built without statistics. The default is neither. To get both, use
# DYLP_OPTIONS ?= stats info. (Use `:=' to override an existing definition.)

export DYLP_OPTIONS ?= stats info


##############################################################################
# You should not need to edit below this line.
##############################################################################

export LIBNAME := OsiDylp

export LIBSRC :=
LIBSRC += OsiDylpMessages.cpp \
	  OsiDylpPresolve.cpp \
	  OsiDylpSolverInterface.cpp \
	  OsiDylpWarmStartBasis.cpp

# The location of the customized Makefiles

export CoinDir ?= $(shell cd ../..; pwd)
export MakefileDir ?= $(CoinDir)/Makefiles

include ${MakefileDir}/Makefile.coin

# This is normally defined in Makefile.lib. To make it easier to add extra
# actions for install and clean, I've redefined it here. The two definitions
# should match, obviously.

InstallDir := $(CoinDir)

# In the event that we somehow arrive here without the proper Dylp defines,
# add them.

ifeq (,$(filter COIN_USE_DYLP, $(CXXFLAGS)))
  CXXFLAGS += $(patsubst %, -D%, $(DylpDefine) $(OsiDylpDefine))
  CXXFLAGS += $(patsubst %, -I%, $(DylpIncDir) $(OsiDylpIncDir))
endif

# This bit of makefile is adapted from Utils/Makefile.Config in the dylp
# distribution. It's duplicated here because we're not guaranteed that dylp
# is actually present.

# Adjust CXXFLAGS correspond to possible OsiDylp compile time
# options as requested by DYLP_OPTIONS. Known keywords are:
#  * stats	enable statistics collection
#  * info	enable informational printing

# info is a bit tricky --- enabling info in dylp means *not* defining NDEBUG.
# Since COIN boilerplate will insert -DNDEBUG in CXXFLAGS when the user asks
# for optimisation, we need to strip it. In the case where we don't want
# info prints, no need to add it twice.

ifdef DYLP_OPTIONS
  ifeq ($(filter stats,$(DYLP_OPTIONS)),stats)
    CXXFLAGS += -DDYLP_STATISTICS
  endif
  ifeq ($(filter info,$(DYLP_OPTIONS)),info)
    CXXFLAGS := $(filter-out -DNDEBUG, $(CXXFLAGS))
  else
    ifeq ($(filter -DNDEBUG, $(CXXFLAGS)),)
      CXXFLAGS += -DNDEBUG
    endif
  endif
endif


###############################################################################

.DELETE_ON_ERROR:

.PHONY: default install library libDylp libOsi clean doc

# Define this to get an mps dump of problem if dylp fails

# CXXFLAGS += -DDYLP_POSTMORTEM

default: install

# Make libDylp if required, then libOsi and libOsiDylp.

libDylp: LIBNAME := Dylp
libDylp:
	@ (cd $(DylpHome)/Dylp && $(MAKE) -f Makefile install)

libOsi:
	@ (cd $(CoinDir)/Osi && $(MAKE) -f Makefile.Osi install)

library : libDylp libOsi
	@ $(MAKE) -f ${MakefileDir}/Makefile.lib library 

# The general sequence here is to make sure the main OSI library is installed,
# then install OsiDylp. Make sure the copies of .hpp files in
# $(InstallDir)/include are writeable before Makefile.lib attempts to delete
# them. (A hedge against version control systems like sccs that keep files as
# read-only unless they're checked out for editing.) It's necessary to remove
# the existing libOsiDylp so that Makefile.lib will refresh it. (Otherwise
# you can be left with a binary for the wrong target architecture.)

install : library
	@ ( cd $(CoinIncDir) && \
	  ( [ -f OsiDylpSolverInterface.hpp ] && \
	   chmod a+w OsiDylp*.hpp ) || \
	   true )
	@ ( cd $(CoinLibDir) && rm -f libOsiDylp* )
	@ $(MAKE) -f $(MakefileDir)/Makefile.lib install

# Make only libOsiDylp.

libOsiDylp : libOsi
	@ $(MAKE) -f ${MakefileDir}/Makefile.lib library

# Make a decent clean that works for just OsiDylp. Running clean up in Osi will
# work, but it also cleans out all other solvers and the main OSI.  We must
# remove the libraries and include files from from the central COIN
# directories, otherwise changes will not be seen.

clean : 
	@ rm -rf $(UNAME)*
	@ rm -rf dep
	@ rm -rf Doc
	@ (cd $(CoinLibDir) && rm -f libOsiDylp*)
	@ (cd $(CoinIncDir) && rm -f OsiDylp*.hpp)

doc :
	@ doxygen $(MakefileDir)/doxygen.conf
