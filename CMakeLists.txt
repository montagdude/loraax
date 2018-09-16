cmake_minimum_required(VERSION 3.0)
set(LORAAX_VERSION 0.1)
project(loraax CXX)

# Executable
include_directories(include)
file(GLOB SOURCES "src/*.cpp")
add_executable(loraax ${SOURCES})

# Set default flags
add_definitions(-DLORAAX_VERSION=\"${LORAAX_VERSION}\")

set(ENABLE_OPENMP TRUE
    CACHE BOOL "Whether to build with OpenMP support.")
set(OPENMP_FLAG "-fopenmp"
    CACHE STRING "Compiler flag for OpenMP.")
set(BUILD_DOCS TRUE
    CACHE BOOL "Whether to build documentation.")

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

# MinGW build definitions
if (MINGW)
	add_definitions(-DISMINGW=1)
else (MINGW)
	add_definitions(-UISMINGW)
endif (MINGW)

# Paths
set(DATADIR ${CMAKE_INSTALL_PREFIX}/share/loraax
    CACHE PATH "Directory to install documentation and examples.")

# Find packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(tinyxml2 REQUIRED)
if (TINYXML2_FOUND)
	include_directories(${TINYXML2_INCLUDE_DIR})
	target_link_libraries(loraax ${TINYXML2_LIBRARY})
endif (TINYXML2_FOUND)
find_package(LibXfoil REQUIRED)
if (LIBXFOIL_FOUND)
	include_directories(${LIBXFOIL_INCLUDE_DIR})
	target_link_libraries(loraax ${LIBXFOIL_LIBRARY} gfortran quadmath)
endif (LIBXFOIL_FOUND)
target_link_libraries(loraax Eigen3::Eigen) 

# Optionally build documentation (needs pdflatex)
if (BUILD_DOCS)
	find_package(LATEX COMPONENTS PDFLATEX)
	if (NOT LATEX_FOUND)
		message(FATAL_ERROR "BUILD_DOCS is enabled but could not find pdflatex.")
	endif (NOT LATEX_FOUND)
	# Build 3 times to ensure all references are resolved
	add_custom_command(OUTPUT user_guide.pdf
	                   COMMAND ${PDFLATEX_COMPILER} "${CMAKE_SOURCE_DIR}/doc/user_guide.tex"
	                   COMMAND ${PDFLATEX_COMPILER} "${CMAKE_SOURCE_DIR}/doc/user_guide.tex"
	                   COMMAND ${PDFLATEX_COMPILER} "${CMAKE_SOURCE_DIR}/doc/user_guide.tex")
	add_custom_target(user_guide ALL
	                  DEPENDS user_guide.pdf)
endif (BUILD_DOCS)

# Install
install(TARGETS loraax DESTINATION bin)
install(FILES
        "sample_cases/StreakWing.xml"
        "sample_cases/StreakWing_inputs.xml"
        "sample_cases/clarky.dat"
        "sample_cases/fx63137sm.dat"
        "sample_cases/mh32_smoothed.dat"
        "sample_cases/naca0012.xml"
        "sample_cases/naca0012_inputs.xml"
        "sample_cases/wonky_aircraft.xml"
        "sample_cases/wonky_inputs.xml"
        "sample_cases/wortmann_ar4.xml"
        "sample_cases/wortmann_ar4_inputs.xml"
        DESTINATION ${DATADIR}/sample_cases)
if (BUILD_DOCS)
	install(FILES ${CMAKE_BINARY_DIR}/user_guide.pdf DESTINATION ${DATADIR})
endif (BUILD_DOCS)
