#! /bin/bash

FC=gfortran

SRC=analysis

FFLAGS=' '

#$FC -J bin -c src/statistics.f90 -o bin/statistics.o
$FC $FFLAGS -J bin -c src/number2string_mod.f90 -o bin/number2string_mod.o
$FC $FFLAGS -I bin -I include -c $SRC/analysis.f90 -o bin/analysis.o
$FC bin/number2string_mod.o bin/analysis.o -L lib -lstats -o bin/analysis.exe

bin/analysis.exe <<< "analysis_parameters.dat"
