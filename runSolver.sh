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

idCase=${HOME}'/Projects/TEPEM4ED/NumericalExamples/StenosedArtery'

export HWLOC_HIDE_ERRORS=2

cd PySolver
	## Compile
cd ..
cp ./PySolver/PySolverTEPEM.py .

cp $idCase'/Mesh.txt' .
cp $idCase'/Basparam.txt' .

## max number of processors -> nproc
np=7; mpirun -np $np --use-hwthread-cpus python3 PySolverTEPEM.py



rm Mesh.txt