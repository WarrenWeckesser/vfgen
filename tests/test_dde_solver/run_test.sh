set -e

echo "Fetching dde_solver_m.f90 from github."
curl -s -S -O https://raw.githubusercontent.com/WarrenWeckesser/dde_solver/master/dde_solver_m.f90

echo "Running vfgen to generate Fortran routines for DDE_SOLVER."
vfgen dde_solver ../vf/sdd.vf

echo "Showing gfortran version."
gfortran --version

echo "Building test_dde_solver."
gfortran -g -c dde_solver_m.f90
gfortran -g -c sdd.f90
gfortran -g test_dde_solver.f90 sdd.o dde_solver_m.o -o test_dde_solver

echo "Running test_dde_solver."
./test_dde_solver 2>out

echo "Showing out."
cat out

if [ -s out ]; then
    exit -1
fi
