flags='-Og -g -fopenmp -fcheck=all -fbounds-check -o'
flags0='-O3 -g -fopenmp -o'
stat=''

## Post-processing files
gfortran Post3D_VolumeData.f90 Common_GeometricalBasis.f90 Common_FieldBasis.f90 $flags0 Post3D_VolumeData
