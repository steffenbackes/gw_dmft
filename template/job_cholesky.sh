#! /bin/bash
#
# #SBATCH --partition=cpu_seq
#SBATCH --partition=cpu_shared  
# #SBATCH --partition=cpu_dist

#SBATCH --ntasks=20
# #SBATCH --nodes=1
# #SBATCH --ntasks-per-node=1
# #SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=99:00:00

#SBATCH --job-name=DMFT_cluster
#SBATCH --error=job.err
#SBATCH --output=job.out
#SBATCH --account=scem
#
#
export OMP_NUM_THREADS=1

module load Eigen
module load gcc
module load fftw
module load hdf5
module load nfft
module load openmpi
module load lapack
module load openblas
module load python

h5dump

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/mnt/beegfs/softs/opt/gcc_10.2.0/lapack/3.10.0/lib64/"
gw_dmft
