# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/caponaroni/hapticmouse/firmware/m1_logger")
  file(MAKE_DIRECTORY "/Users/caponaroni/hapticmouse/firmware/m1_logger")
endif()
file(MAKE_DIRECTORY
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/m1_logger"
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix"
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/tmp"
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/src/m1_logger-stamp"
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/src"
  "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/src/m1_logger-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/src/m1_logger-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/caponaroni/hapticmouse/firmware/m1_logger/build/_sysbuild/sysbuild/images/m1_logger-prefix/src/m1_logger-stamp${cfgdir}") # cfgdir has leading slash
endif()
