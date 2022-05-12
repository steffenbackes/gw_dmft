#!/bin/bash
#MSUB -r gen_h-vxc_dc2shiftw0                # jobname
#MSUB -n  300                      # No. tasks
###MSUB -N 1                       # No. nodes 
###MSUB -c 1                           # Cores per task
###MSUB -M 4000                    # Memory per core in Mbyte
#MSUB -Q long                     # this allows 3 days max, otherwise 24h
#MSUB -T 259000                # Max. runtime in sec.
#MSUB -o job.out
#MSUB -e job.err
#MSUB -q skylake
####MSUB -q hybrid
####MSUB -q xlarge
#MSUB -A gen1393
#MSUB -m scratch

set -x
cd ${BRIDGE_MSUB_PWD}
export OMP_NUM_THREADS=1
#
#
module load hdf5

#h5dump --help

gw_dmft
