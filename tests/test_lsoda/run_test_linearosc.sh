set -e

echo "Running vfgen to generate Fortran routines for LSODA."
vfgen lsoda ../vf/linearosc.vf

echo "Building test_lsoda."
gfortran --version
gfortran -c -std=legacy -w opkdmain.f opkda1.f opkda2.f
gfortran -c test_lsoda_linearosc.f
gfortran -c linearosc_rhs.f
gfortran test_lsoda_linearosc.o linearosc_rhs.o opkdmain.o opkda1.o opkda2.o -o test_lsoda_linearosc

echo "Running test_lsoda."
./test_lsoda_linearosc > out 2>&1

if [ -s out ]; then
    exit -1
fi
