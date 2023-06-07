import numpy as np

A0 = [0,1,2,3,4]
A1 = [5,6,7,8,9]
Atoms = [A0,A1]

for A in Atoms:
   # hybridization function
   inf = open('delta0_matrix_solver.dat','r')
   outA = open('delta0_matrix_solver_A'+str(Atoms.index(A))+'.dat','w')
   for line in inf:
      data = line.split()
      if int(data[1]) in A and int(data[2]) in A:
         outA.write(line)
   inf.close()
   outA.close()

   # local potential 
   inf = open('mu_matrix_solver.dat','r')
   outA = open('mu_matrix_solver_A'+str(Atoms.index(A))+'.dat','w')
   for line in inf:
      data = line.split()
      if int(data[0]) in A and int(data[1]) in A:
         outA.write(line)
   inf.close()
   outA.close()

