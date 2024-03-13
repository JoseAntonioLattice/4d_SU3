#! /bin/bash

FC=gfortran

SRC=analysis

$FC -o bin/analysis.exe $SRC/analysis.f90 

bin/analysis.exe <<< input_parameters.par

