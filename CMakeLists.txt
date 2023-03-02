cmake_minimum_required(VERSION 3.13)

project(MassCCS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()

  set(doxyfile ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

  add_custom_target(docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# RelWithDebInfo use O3 by default to allow vectorization profiling
option(RELWITHDEBINFO_FORCE_O3 "Force O3 with RelWithDebInfo" ON)
if(RELWITHDEBINFO_FORCE_O3)
  string(REPLACE "-O2" "-O3" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
endif()

option(VERBOSE_VECTORIZE "Show vectorization report" OFF)
if(VERBOSE_VECTORIZE)
  if("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")
  elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopt-info-vec-missed")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-missed")
  elseif("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -opt-report5")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -opt-report5")
  else()
    message("Disable vectorization report: Unknown compiler ${CMAKE_C_COMPILER_ID}")
  endif()
endif()

option(USE_FASTMATH "Enable fast-math optimizations" ON)
if(NOT USE_FASTMATH)
  message("Manually disabled fast math compilation optimizations")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")
else()
  message("Not building with fast-math optimizations: compiler not supported.")
endif()

option(USE_OPENMP "Enable OpenMP" ON)
if(USE_OPENMP)
  find_package(OpenMP QUIET)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else()
    message("Not building with multithreading support: OpenMP not found in system.")
  endif()
else()
  message("Not building with multithreading support: manually disabled.")
endif()

option(USE_NATIVE "Build with -march=native" ON)
if(USE_NATIVE)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    message("Building with -march=native")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  else()
    message("Not building with -march=native: not supported by compiler")
  endif()
else()
  message("Not building with -march=native: not enabled.")
endif()

file(GLOB
  CHECK_SOURCE_FILES
  src/*.cpp
  src/headers/*.h
)
include("cmake/ClangFormat.cmake")

include_directories(src/headers)

add_executable(massccs
  src/main.cpp
  src/System.cpp
  src/Input.cpp
  src/RandomNumber.cpp
  src/Math.cpp
  src/Time.cpp
  src/MoleculeTarget.cpp
  src/GasBuffer.cpp
  src/Equipotential.cpp
  src/LinkedCell.cpp
  src/Force.cpp
)

