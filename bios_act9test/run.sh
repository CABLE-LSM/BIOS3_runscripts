#PBS -N bios_act9test
#PBS -l wd
#PBS -P tm70
#PBS -q normal
#PBS -l walltime=00:30:00
#PBS -l ncpus=4
#PBS -l mem=5GB
#PBS -l storage=gdata/rp23
#PBS -j oe

./setup.sh
module load intel-mpi/2019.5.281
mpiexec -n $PBS_NCPUS ./cable-mpi > logs/log_out_cable.txt
