#!/bin/sh
#  This is automatically generated file.
#  Edit "cmake/setup.sh.cmake" file.

# set environment variables LD_LIBRARY_PATH and DYLD_LIBRARY_PATH

@LD_LIBRARY_PATH@_OLD=${@LD_LIBRARY_PATH@}
@DYLD_LIBRARY_PATH@_OLD=${@DYLD_LIBRARY_PATH@}

@LD_LIBRARY_PATH@="@ROOT_LIBRARY_DIR@:@CLHEP_LIBRARY_DIR@:@CMAKE_INSTALL_PREFIX@/@install_lib_dir@"
if [ -n "${@LD_LIBRARY_PATH@_OLD}" ]; then
  @LD_LIBRARY_PATH@="${@LD_LIBRARY_PATH@}:${@LD_LIBRARY_PATH@_OLD}"
fi
export @LD_LIBRARY_PATH@

@DYLD_LIBRARY_PATH@="@ROOT_LIBRARY_DIR@:@CLHEP_LIBRARY_DIR@:@CMAKE_INSTALL_PREFIX@/@install_lib_dir@"
if [ -n "${@DYLD_LIBRARY_PATH@_OLD}" ]; then
  @DYLD_LIBRARY_PATH@="${@DYLD_LIBRARY_PATH@}:${@DYLD_LIBRARY_PATH@_OLD}"
fi
export @DYLD_LIBRARY_PATH@

