source /usr/share/Modules/init/csh

#module purge
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load netcdf/oneapi-2024.2/hdf5-1.14.4/intel-mpi/4.9.2
module load hdf5/oneapi-2024.2/intel-mpi/1.14.4
module load intel-mkl/2024.2

which mpiifort
which mpiicc

#module load anaconda3/2024.10 #to use python/xarray
