#!/bin/bash

set -xeuo pipefail

export FC=${FC:-gfortran-10}

cmake -B build
cmake --build build
