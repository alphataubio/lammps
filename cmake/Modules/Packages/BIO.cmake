
find_package(libsbml QUIET)
if(NOT libsbml_FOUND)
  set(libsbml_URL "https://github.com/sbmlteam/libsbml/archive/refs/tags/v5.20.4.tar.gz" CACHE STRING "URL for libsbml tarball")
  set(libsbml_MD5 "98a04a2806d500bbc18c68ae2969b88c" CACHE STRING "MD5 checksum of libsbml tarball")
  mark_as_advanced(libsbml_URL)
  mark_as_advanced(libsbml_MD5)

  set(WITH_CPP_NAMESPACE ON)
  set(WITH_SWIG OFF)
  mark_as_advanced(WITH_CPP_NAMESPACE)
  mark_as_advanced(WITH_SWIG)

  # download and build a local copy of libyaml
  include(ExternalCMakeProject)
  ExternalCMakeProject(libsbml ${libsbml_URL} ${libsbml_MD5} libsbml . "")
  include_directories(${libsbml_SOURCE_DIR}/src)
  include_directories(${libsbml_BINARY_DIR}/src)

  add_library(libsbml ALIAS sbml-static)

endif()

target_link_libraries(lammps PRIVATE libsbml)
