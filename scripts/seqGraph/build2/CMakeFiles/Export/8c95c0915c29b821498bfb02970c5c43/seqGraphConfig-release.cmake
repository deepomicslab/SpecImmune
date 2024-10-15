#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "seqGraph" for configuration "RELEASE"
set_property(TARGET seqGraph APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(seqGraph PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libseqGraph.so.1.0.1"
  IMPORTED_SONAME_RELEASE "libseqGraph.so.1"
  )

list(APPEND _cmake_import_check_targets seqGraph )
list(APPEND _cmake_import_check_files_for_seqGraph "${_IMPORT_PREFIX}/lib64/libseqGraph.so.1.0.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
