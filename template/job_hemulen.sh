#PBS -S /bin/bash
#PBS -q GroupA

#PBS -l nodes=2:ppn=8
#PBS -l walltime=48:00:00
# #PBS -l mem=64gb
#PBS -N SrVO3_dmft

cd $PBS_O_WORKDIR/
. /opt/intel/composer_xe_2015/bin/compilervars.sh intel64 &>> /dev/null
. /opt/intel/impi/5.1.3.210/intel64/bin/mpivars.sh
. /opt/CTQMC/env.sh

gw_dmft > out.gwdmft 2>&1
