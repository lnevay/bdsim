if (DEFINED ENV{ROOTSYS})
  message(STATUS "Use ROOTSYS from environment: $ENV{ROOTSYS}")
  set (ROOTSYS $ENV{ROOTSYS})
elseif(USE_AFS)
  # add ROOT afs library to PATH (only in cmake)
  if (APPLE)
    set(ROOTSYS "/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-mac1011-clang70-opt/root")
  elseif (RHL6)
    set(ROOTSYS "/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root")
  endif()
  if($ENV{VERBOSE})
    message(STATUS "ROOTSYS: ${ROOTSYS}")
  endif()
endif()

find_package(ROOT REQUIRED)
# major ROOT version
string(SUBSTRING ${ROOT_VERSION} 0 1 ROOT_MAJOR_VERSION)
# add ROOT include directory
include_directories(${ROOT_INCLUDE_DIR})

option(ROOT_DOUBLE_OUTPUT "Double precision root output" OFF)
if(ROOT_DOUBLE_OUTPUT)
  add_definitions("-D__ROOTDOUBLE__")
  set(PREPROCESSOR_DEFS "-D__ROOTBUILD__;-D__ROOTDOUBLE__")
else()
  set(PREPROCESSOR_DEFS "-D__ROOTBUILD__")
endif()

# Make Dictionaries
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/root)
file(GLOB linkHeaders ${CMAKE_CURRENT_SOURCE_DIR}/include/*LinkDef.hh)
# for loop over link definitions
foreach(header ${linkHeaders})
  # remove LinkDef.hh
  string(FIND ${header} "LinkDef.hh" pos REVERSE)
  string(FIND ${header} "/" dir REVERSE)
  MATH(EXPR beginpos "${dir}+1")
  MATH(EXPR length "${pos}-${beginpos}")
  string(SUBSTRING ${header} ${beginpos} ${length} className)
  #message(STATUS "${ROOTCINT_EXECUTABLE} -f ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc -c ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}")
  if (${ROOT_MAJOR_VERSION} STREQUAL "6")
    # add preprocessor flag __ROOTBUILD__ as used in BDSIM.
    # root6 .pcm file needs to be at library level
    # for discussion see e.g. this forum topic: https://root.cern.ch/phpBB3/viewtopic.php?t=17654
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc
      ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict_rdict.pcm
      COMMAND ${ROOTCINT_EXECUTABLE}
      ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc -noIncludePaths -inlineInputHeader
      ${PREPROCESSOR_DEFS}  -I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_SOURCE_DIR}/parser/
      ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      COMMENT "Generate ROOT Dictionary for ${className}"
      )
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${className}Dict_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc
      )
    # add arbitrary dependence to make it run at compile time (since DEPENDS at custom_command seems not to work)
    set_property(
      SOURCE ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc
      APPEND PROPERTY OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${className}Dict_rdict.pcm
      )
    # keep list of ROOT pcms
    set(root_pcms ${root_pcms} ${CMAKE_CURRENT_BINARY_DIR}/${className}Dict_rdict.pcm)
  else()
    # ROOT 5
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.h
      COMMAND ${ROOTCINT_EXECUTABLE}
      ARGS -f ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc
      -c -p -I${CMAKE_CURRENT_SOURCE_DIR}
      ${PREPROCESSOR_DEFS}
      ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/include/${className}.hh ${header}
      COMMENT "Generate ROOT Dictionary for ${className}"
      )
    set(root_dicts_headers ${root_dicts_headers} ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.h)
  endif()
  # keep list of ROOT dictionaries
  set(root_files ${root_files} ${CMAKE_CURRENT_SOURCE_DIR}/src/${className}.cc)
  set(root_dicts ${root_dicts} ${CMAKE_CURRENT_BINARY_DIR}/root/${className}Dict.cc)
endforeach()
