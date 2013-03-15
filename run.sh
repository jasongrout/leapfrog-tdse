#!/bin/bash

set -o verbose

gfortran leapfrog.f90 -o leapfrog.o
./leapfrog.o
gnuplot animate.gnuplot