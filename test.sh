#!/bin/bash

set -xeuo pipefail

export FC=${FC:-gfortran-10}

# FENGSHA
make -C src clean
make -C src

# Wrapper
$FC -c driver/ccpp/catchem_dust_wrapper.F90 -I./src
