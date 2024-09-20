set -e

echo "Running vfgen evf ../vf/rossler.vf > rossler_evf.vf"
vfgen evf ../vf/rossler.vf > rossler_evf.vf

echo "Running vfgen gsl rossler_evf.vf"
vfgen gsl rossler_evf.vf

echo "Building test_evf"
gcc -c rossler_evf_gvf.c
gcc -c test_evf.c 
gcc -o test_evf test_evf.o rossler_evf_gvf.o -lgsl -lgslcblas -lm

echo "Running test_evf"
./test_evf
