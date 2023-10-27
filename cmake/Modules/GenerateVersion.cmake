#######################################################################
# generate git-version signature
#######################################################################
include_directories(${oofem_BINARY_DIR})
find_package(Git)
add_custom_target(version
  ${CMAKE_COMMAND} -D SRC=${oofem_SOURCE_DIR}/src/oofem_version.h.in
                   -D DST=${oofem_BINARY_DIR}/oofem_version.h
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -P ${CMAKE_SOURCE_DIR}/src/GenerateVersionHeader.cmake
  )