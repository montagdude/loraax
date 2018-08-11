cmake_minimum_required(VERSION 3.0)
set(LORAAX_VERSION 0.1)
project(loraax)

# Set default flags
add_definitions(-DLORAAX_VERSION=\"${LORAAX_VERSION}\")

set(ENABLE_OPENMP TRUE
    CACHE BOOL "Whether to build with OpenMP support.")
set(OPENMP_FLAG "-fopenmp"
    CACHE STRING "Compiler flag for OpenMP.")

if (NOT ENABLE_OPENMP)
	set(OPENMP_FLAG "")
endif (NOT ENABLE_OPENMP)

if (CMAKE_BUILD_TYPE MATCHES "Release")
	add_definitions(-UDEBUG)
	if (NOT ENABLE_OPENMP)
		message(WARNING "Disabling OpenMP support since ENABLE_OPENMP=FALSE.")
	endif (NOT ENABLE_OPENMP)
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 ${OPENMP_FLAG}")
elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
	add_definitions(-DDEBUG)
	if (ENABLE_OPENMP)
		message(WARNING "Enabling OpenMP for debug build since ENABLE_OPENMP=TRUE.")
	endif (ENABLE_OPENMP)
	set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall ${OPENMP_FLAG}")
else (CMAKE_BUILD_TYPE MATCHES "Release")
	set(CMAKE_CXX_FLAGS "-O2 ${OPENMP_FLAG}")
endif (CMAKE_BUILD_TYPE MATCHES "Release")

# Executable
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_executable(loraax ${SOURCES})

# Find required packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(tinyxml2 REQUIRED)
find_package(LibXfoil REQUIRED)
if (LIBXFOIL_FOUND)
	include_directories(${LIBXFOIL_INCLUDE_DIR})
	target_link_libraries(loraax ${LIBXFOIL_LIBRARY})
endif (LIBXFOIL_FOUND)
target_link_libraries(loraax Eigen3::Eigen tinyxml2) 

# Install
install(TARGETS loraax DESTINATION bin)
