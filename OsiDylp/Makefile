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

# export OptLevel ?= -g
export OptLevel ?= -O

# If OsiDylp is being built as a dependency for the Dylp sample programs, the
# symbol DYLP_OPTIONS will be defined on entry to this makefile (but may well
# be defined to be empty, hence the use of the origin function). If it's not
# defined, the two options that make a difference for OsiDylp are stats and
# info. Request info if you want information prints, and stats if you want
# statistics. Note that both OsiDylp and dylp must be built with statistics
# enabled before you'll actually see any stats. dylp will ignore an incoming
# statistics structure if built without statistics, and OsiDylp will not supply
# one if built without statistics. The default is neither info or stats.

# If you are triggering the build of Dylp from here (in other words, if you
# haven't already issued a `make samples' command in the Dylp distribution
# directory), the entire Dylp build will be controlled by the options you
# request here. If you want separate build options for the two, the best way
# is to build dylp first (but *not* the sample programs, which require
# OsiDylp). Then build OsiDylp, etc.

ifeq ("$(origin DYLP_OPTIONS)","undefined")
# export DYLP_OPTIONS ?= paranoia stats info
  export DYLP_OPTIONS ?=
endif


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

# Having trouble with the build? Changes to .hpp files don't seem to be
# taking effect? That's because by default you're building with the cached
# copies of the include files stashed in ${CoinDir}/include. Generally, this
# is a good thing --- this build sees the same .hpp files the rest of COIN is
# using.  If you'd like to play with changes before exposing them to the world,
# uncomment the following line.

# CXXFLAGS += -I./include

ifeq (,$(filter COIN_libDylp, $(CoinLibsDefined)))
  NO_DYLP_MSG := COIN_HAS_DYLP is not defined in Makefiles/Makefile.location.
  NO_DYLP_MSG += Probably the line 'CoinLibsDefined += COIN_libDylp' is
  NO_DYLP_MSG += commented out.
  $(error $(NO_DYLP_MSG))
endif

# Add the defines and includes for Dylp and OsiDylp.

CXXFLAGS += $(patsubst %, -D%, $(DylpDefine) $(OsiDylpDefine))
CXXFLAGS += $(patsubst %, -I%, $(DylpIncDir) $(OsiDylpIncDir))

# This bit of makefile is adapted from Utils/Makefile.Config in the dylp
# distribution.
# Adjust CXXFLAGS to correspond to the possible OsiDylp compile time
# options as requested by DYLP_OPTIONS. Known keywords are:
#  * stats	enable statistics collection
#  * info	enable informational printing
#  * paranoia	enable paranoid checks

# info is a bit tricky --- enabling info in dylp means *not* defining NDEBUG.
# Since COIN boilerplate will insert -DNDEBUG in CXXFLAGS when the user asks
# for optimisation, we need to strip it. Avoid duplicating existing defines.

ifdef DYLP_OPTIONS
  ifeq ($(filter paranoia,$(DYLP_OPTIONS)),paranoia)
    ifeq ($(filter -DPARANOIA, $(CXXFLAGS)),)
      CXXFLAGS += -DPARANOIA
    endif
  endif
  ifeq ($(filter stats,$(DYLP_OPTIONS)),stats)
    ifeq ($(filter -DDYLP_STATISTICS, $(CXXFLAGS)),)
      CXXFLAGS += -DDYLP_STATISTICS
    endif
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
