#!/bin/bash
echo "Filling:"
h5dump -d /densities input.out.h5 | tail -n 5 | head -n 2

echo "Sign:"
h5dump -d /Sign input.out.h5 | tail -n 4 | head -n 1

h5dump -d /SigmaG_omega  input.out.h5 > sgf.dat
h5dump -d /SigmaG_legendre input.out.h5 > gf_legendre.dat
cp ../g0R_DCA.dat ./g0_solver.dat
python getoutput_ctint.py
