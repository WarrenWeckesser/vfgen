set -e

echo "Running vfgen to generate Fortran routines for LSODA."
vfgen lsoda ../vf/pidecay.vf

echo "Building test_lsoda."
gfortran --version
gfortran -c -std=legacy -w opkdmain.f opkda1.f opkda2.f
gfortran -c test_lsoda_pidecay.f
gfortran -c pidecay_rhs.f
gfortran test_lsoda_pidecay.o pidecay_rhs.o opkdmain.o opkda1.o opkda2.o -o test_lsoda_pidecay

echo "Running test_lsoda_pidecay."
./test_lsoda_pidecay > out 2>&1

if [ -s out ]; then
    exit -1
fi
