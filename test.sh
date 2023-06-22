#!/bin/bash

set -xeuo pipefail

export FC=${FC:-gfortran}

cmake -B build
cmake --build build
