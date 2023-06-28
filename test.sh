#!/bin/bash

set -xeuo pipefail

export FC=${FC:-gfortran}

BUILD_DIR=build
cmake -B $BUILD_DIR
cmake --build $BUILD_DIR
ctest --test-dir $BUILD_DIR
