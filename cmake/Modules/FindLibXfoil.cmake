find_path(LIBXFOIL_INCLUDE_DIR libxfoil.h)
find_library(LIBXFOIL_LIBRARY NAMES xfoil)

if (LIBXFOIL_INCLUDE_DIR AND LIBXFOIL_LIBRARY)
	set(LIBXFOIL_FOUND TRUE)
endif (LIBXFOIL_INCLUDE_DIR AND LIBXFOIL_LIBRARY)

if (LIBXFOIL_FOUND)
	message(STATUS "Found libxfoil header in ${LIBXFOIL_INCLUDE_DIR}")
	message(STATUS "Found libxfoil library: ${LIBXFOIL_LIBRARY}")
else (LIBXFOIL_FOUND)
	if (NOT LIBXFOIL_INCLUDE_DIR)
		message(FATAL_ERROR "Could NOT find libxfoil.h!")
	endif (NOT LIBXFOIL_INCLUDE_DIR)

	if (NOT LIBXFOIL_LIBRARY)
		message(FATAL_ERROR "Could NOT find libxfoil library!")
	endif (NOT LIBXFOIL_LIBRARY)
endif (LIBXFOIL_FOUND)
