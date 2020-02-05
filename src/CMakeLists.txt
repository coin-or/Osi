set(NAME "OSI")

# PTHREAD
if(PTHREADS_FOUND)
  set(${NAME}_PTHREADS "1" CACHE INTERNAL "Use pthread")
endif()

set(COIN_${NAME}_CHECKLEVEL "0" CACHE INTERNAL
  "${NAME} check level")
set(COIN_${NAME}_VERBOSITY "0" CACHE INTERNAL
  "${NAME} verbosity level")
configure_file(config.h.cmake.in config.h)
configure_file(config_osi.h.cmake.in config_osi.h)

set(_SRCS
  Osi/OsiAuxInfo.cpp
  Osi/OsiBranchingObject.cpp
  Osi/OsiChooseVariable.cpp
  Osi/OsiColCut.cpp
  Osi/OsiCut.cpp
  Osi/OsiCuts.cpp
  Osi/OsiNames.cpp
  Osi/OsiPresolve.cpp
  Osi/OsiRowCut.cpp
  Osi/OsiRowCutDebugger.cpp
  Osi/OsiSolverBranch.cpp
  Osi/OsiSolverInterface.cpp)

set(_HDRS
  Osi/OsiConfig.h
  Osi/OsiAuxInfo.hpp
  Osi/OsiBranchingObject.hpp
  Osi/OsiChooseVariable.hpp
  Osi/OsiColCut.hpp
  Osi/OsiCollections.hpp
  Osi/OsiCut.hpp
  Osi/OsiCuts.hpp
  Osi/OsiPresolve.hpp
  Osi/OsiRowCut.hpp
  Osi/OsiRowCutDebugger.hpp
  Osi/OsiSolverBranch.hpp
  Osi/OsiSolverInterface.hpp
  Osi/OsiSolverParameters.hpp)

add_library(Osi ${_SRCS} ${_HDRS})
target_include_directories(Osi PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(Osi
  PUBLIC HAVE_CONFIG_H
  PRIVATE OSI_BUILD)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET Osi PROPERTY CXX_STANDARD 11)
  set_property(TARGET Osi PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(Osi PUBLIC cxx_std_11)
endif()
target_link_libraries(Osi PUBLIC Coin::CoinUtils)
set_target_properties(Osi PROPERTIES
  PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_osi.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::Osi ALIAS Osi)

# Install
include(GNUInstallDirs)
install(TARGETS Osi
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
