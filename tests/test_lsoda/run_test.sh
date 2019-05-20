echo "Fetching LSODA from netlib."
curl -s -S -O "http://www.netlib.org/odepack/{opkdmain.f,opkda1.f,opkda2.f}"

echo "Running vfgen to generate Fortran routines for LSODA."
vfgen lsoda ../vf/linearosc.vf

echo "Building test_lsoda."
gfortran -c -w opkdmain.f opkda1.f opkda2.f
gfortran -c test_lsoda.f
gfortran -c linearosc_rhs.f
gfortran test_lsoda.o linearosc_rhs.o opkdmain.o opkda1.o opkda2.o -o test_lsoda

echo "Running test_lsoda."
./test_lsoda > out 2>&1

if [ -s out ]; then
    exit -1
fi
