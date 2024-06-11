#PBS -N bios_1000pts
#PBS -l wd
#PBS -P tm70
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=48
#PBS -l mem=48GB
#PBS -l storage=gdata/rp23
#PBS -j oe

./setup.sh
module load intel-mpi/2019.5.281
mpiexec -n $PBS_NCPUS ./cable-mpi > logs/log_out_cable.txt
