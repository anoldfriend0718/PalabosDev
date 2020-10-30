cmake_minimum_required(VERSION 3.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    message("ccache used.")
endif() 

project(mypalabos)
enable_language(CXX)

# This compiles Palabos library and all the examples
# It is used mainly for Continuous Integration (CI)
set(CI ON CACHE BOOL "Compile Everything (Continuous Integration)")
set(PALABOS_SOURCE_REPO "/home/anoldfriend/Workspace/MyRepo/palabos")
message("PALABOS_SOURCE_REPO:" ${PALABOS_SOURCE_REPO})
set(CUSTOMIZED_PARABOS_SOURCE_REPO ${CMAKE_SOURCE_DIR})
message("CUSTOMIZED_PARABOS_SOURCE_REPO: " ${CUSTOMIZED_PARABOS_SOURCE_REPO})
set(PALABOS_LIB_DIRECTORY_BASE  "${CMAKE_BINARY_DIR}/lib")
message("PALABOS_LIB_DIRECTORY_BASE: " ${PALABOS_LIB_DIRECTORY_BASE})


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
    set(PALABOS_LIB_DIRECTORY  ${PALABOS_LIB_DIRECTORY_BASE}/debug)
    set(PALABOS_LIB_NAME palabosd)    
else()
    set(PALABOS_LIB_DIRECTORY  ${PALABOS_LIB_DIRECTORY_BASE}/release)
    set(PALABOS_LIB_NAME palabosr)
endif()
LINK_DIRECTORIES(${PALABOS_LIB_DIRECTORY})

message("PALABOS_LIB_DIRECTORY:" ${PALABOS_LIB_DIRECTORY})

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif()


###############################################################################
# Palabos Library
###############################################################################
include_directories("${PALABOS_SOURCE_REPO}/src")
include_directories("${PALABOS_SOURCE_REPO}/externalLibraries")
include_directories("${PALABOS_SOURCE_REPO}/externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC ${PALABOS_SOURCE_REPO}/src/*.cpp)
file(GLOB_RECURSE EXT_SRC ${PALABOS_SOURCE_REPO}/externalLibraries/tinyxml/*.cpp)

if(${CMAKE_BUILD_TYPE} STREQUAL Debug)
    set(PALABOS_LIB_NAME palabosd) 
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "lib/debug/")
    add_library(${PALABOS_LIB_NAME} STATIC ${PALABOS_SRC} ${EXT_SRC})
else()
    set(PALABOS_LIB_NAME palabosr)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "lib/release/")
    add_library(${PALABOS_LIB_NAME} STATIC ${PALABOS_SRC} ${EXT_SRC})
endif()

find_path(PLOG_INCLUDE_DIRS "plog/Appenders/AndroidAppender.h")
include_directories(${PLOG_INCLUDE_DIRS})
###############################################################################
# CUSTOMIZED PALABOS DIRECTORIES
###############################################################################
include_directories("${CUSTOMIZED_PARABOS_SOURCE_REPO}/src")

# To link palabos in directories other than this
cmake_policy(SET CMP0079 NEW)

add_subdirectory(apps/cylinder2d)
add_subdirectory(apps/pmsinglephaseflow)

enable_testing()
add_subdirectory(tests)