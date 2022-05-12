#PBS -S /bin/bash
#PBS -q workq

#PBS -l nodes=1:ppn=48
#PBS -l walltime=48:00:00
# #PBS -l mem=64gb
#PBS -N SrVO3_dmft

#PBS -o job.out
#PBS -e job.err
#PBS -k oe 

cd $PBS_O_WORKDIR/

gw_dmft > out.gwdmft 2>&1
