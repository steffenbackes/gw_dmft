#!/bin/bash
rm Swl.dat
rm Gwl.dat
rm nnt.dat
rm nnw.dat
python splitInput.py

for A in 0 1
do
   cp Uijkl_5orb.dat  imp${A}/Uijkl.dat
   cp delta0_matrix_solver_A${A}.dat  imp${A}/delta0_matrix_solver.dat
   cp mu_matrix_solver_A${A}.dat  imp${A}/mu_matrix_solver.dat

   cd imp${A}
   sed -i "s/seed.*/seed=`shuf -i 1-1000 -n 1`/" input.ini
   mv Swl.dat Swl.dat_old
   mv Gwl.dat Gwl.dat_old
   mv input.out.h5 input.out.h5_old

   echo 'Running runsolver.sh, sleep 5'
   sleep 5

   export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/steffen/lib"
   ncpus=`wc -l < $PBS_NODEFILE`
   mpirun -np $ncpus hybmat input.ini > out.solver 2>&1

   echo 'Solver done, sleep 5'
   sleep 5
   ./getoutput_cthyb.sh
   ############################################################
   cd ..

python combineOutput.py

#cp imp/Gwl.dat ./
#cp imp/nnt.dat ./
#cp imp/nnw.dat ./
sleep 5
