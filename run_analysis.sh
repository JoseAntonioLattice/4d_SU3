#! /bin/bash

FC=gfortran

SRC=analysis

$FC -J bin -c src/statistics.f90 -o bin/statistics.o
$FC -J bin -c src/number2string_mod.f90 -o bin/number2string_mod.o
$FC -I bin -c $SRC/analysis.f90 -o bin/analysis.o
$FC  bin/statistics.o bin/number2string_mod.o bin/analysis.o -o bin/analysis.exe  

bin/analysis.exe <<< input_parameters.par

