add_executable(OsiTest)
target_sources(OsiTest PRIVATE
  OsiTestSolver.cpp
  OsiTestSolverInterface.cpp
  OsiTestSolverInterfaceIO.cpp
  OsiTestSolverInterfaceTest.cpp
	unitTest.cpp)
target_include_directories(OsiTest PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(OsiTest
  PUBLIC HAVE_CONFIG_H
  PRIVATE OSI_BUILD)
target_compile_features(OsiTest PUBLIC cxx_std_11)
if(APPLE)
  set_target_properties(OsiTest PROPERTIES
    INSTALL_RPATH "@loader_path;@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(OsiTest PROPERTIES
    INSTALL_RPATH "$ORIGIN:$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
target_link_libraries(OsiTest PRIVATE
  Coin::Osi
  Coin::OsiCommonTest
)
set_target_properties(OsiTest PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_executable(Coin::OsiTest ALIAS OsiTest)
add_test(NAME OsiTest COMMAND OsiTest)
