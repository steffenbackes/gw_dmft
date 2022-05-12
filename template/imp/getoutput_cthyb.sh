#!/bin/bash
echo "Filling:"
h5dump -d /simulation/results/n/mean/value input.out.h5 | tail -n 7 | head -n 3

echo "Sign:"
h5dump -d /Sign input.out.h5 | tail -n 4 | head -n 1

h5dump -d /gf/data  input.out.h5 > gf.dat
h5dump -d /G1_LEGENDRE input.out.h5 > gf_legendre.dat
python getoutput.py
