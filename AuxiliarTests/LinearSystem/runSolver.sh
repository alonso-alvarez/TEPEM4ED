##################################################
#
# python -m pip install [--user] numpy mpi4py
# python -m pip install [--user] petsc petsc4py

#######################
# Executable
#
# 1. install: python -m pip install pyinstaller
# 2. python -m PyInstaller test_parallel.py -> create a dist folder
# 3. dist/test_parallel: mpirun -np 8 --use-hwthread-cpus test_parallel

export HWLOC_HIDE_ERRORS=2

## max number of processors -> nproc --all
npmax=$(nproc)
echo '** Max number of processors: '$npmax

np=4
if [ "$np" -le "$npmax" ]; then
	mpirun -np $np --use-hwthread-cpus python3 PySolverTEPEM.py -n 10
else
	mpirun -np $np --oversubscribe --use-hwthread-cpus python3 PySolverTEPEM.py
fi;