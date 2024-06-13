#! /bin/bash

FC=gfortran

SRC=analysis

#$FC -J bin -c src/statistics.f90 -o bin/statistics.o
$FC -J bin -c src/number2string_mod.f90 -o bin/number2string_mod.o
$FC -I bin -I include -c $SRC/analysis.f90 -o bin/analysis.o
$FC bin/number2string_mod.o bin/analysis.o -L lib -lstats -o bin/analysis.exe   

bin/analysis.exe <<< "analysis_parameters.dat"
