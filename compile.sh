#!/bin/bash
rm -f *.o
rm -f *.mod

OPT="-funroll-loops -O3"
FC=ifort

# OPT="-g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow"
echo "Compile constants.f90..."
$FC $OPT -c constants.f90
echo "Compile params.f90..."
$FC $OPT -c params.f90
echo "Compile mathfunc.f90..."
$FC $OPT -c mathfunc.f90
echo "Compile fourier_interpolation.f90..."
$FC $OPT -c fourier_interpolation.f90
echo "Compile matsfunc_ops.f90..."
$FC $OPT -c matsfunc_ops.f90
echo "Compile gw.f90..."
$FC $OPT -c gw.f90
echo "Compile anacont.f90..."
$FC $OPT -c anacont.f90
echo "Compile io_ops.f90..."
$FC $OPT -c io_ops.f90
echo "Compile solver.f90..."
$FC $OPT -c solver.f90
echo "Compile observables.f90..."
$FC $OPT -c observables.f90
echo "Compile main.f90..."
$FC $OPT -c main.f90

echo "Linking main executable..."
$FC constants.o params.o mathfunc.o fourier_interpolation.o matsfunc_ops.o gw.o anacont.o io_ops.o solver.o observables.o main.o -o gw_dmft -llapack
#$FC constants.o params.o mathfunc.o fourier_interpolation.o matsfunc_ops.o gw.o anacont.o io_ops.o solver.o observables.o main.o -o gw_dmft_unfold2 -llapack

rm -f *.o
rm -f *.mod

# ifort -llapack constants.o io_ops.o matsfunc_ops.o solver.o main.o -o main
