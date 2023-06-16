#!/bin/bash
noAtoms=1
noCores=48

for a in {1..${noAtoms}}
do
   cd imp${a}
   sed -i "s/seed.*/seed=`shuf -i 1-1000 -n 1`/" input.ini
   sed -i "s/SEED.*/SEED = `shuf -i 1-1000 -n 1`/" alps_parm.dat
   mv Swl.dat Swl.dat_old
   mv Gwl.dat Gwl.dat_old
   mv nnt.dat nnt.dat_old
   mv nnw.dat nnw.dat_old
   mv input.out.h5 input.out.h5_old
   
   echo 'Running solver for atom '${a}
   sleep 5

   ############################################################
   #mpirun -np ${noCores} hybmat input.ini > out.solver 2>&1
   #mpirun -np ${noCores} ctint_real input.ini > out.solver 2>&1
   mpirun -np ${noCores} alps_cthyb alps_parm.dat > out.solver 2>&1

   echo 'Solver done, sleep 5'
   sleep 5
   #./getoutput_cthyb.sh
   #./getoutput_ctint.sh
   python getoutput_seg.py
   ############################################################
   
   cd ..
done
