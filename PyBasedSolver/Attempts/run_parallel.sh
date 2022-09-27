# python -m pip install [--user] numpy mpi4py
# python -m pip install [--user] petsc petsc4py

np=8

mpirun -np $np --use-hwthread-cpus python3 test_parallel.py