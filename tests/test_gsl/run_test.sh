set -e

echo "Installing libgsl"
apt-get -yq install libgsl-dev

echo "Running vfgen to generated code for gsl"
vfgen gsl ../vf/linearosc.vf

echo "Building test_gsl"
gcc -c linearosc_gvf.c
gcc -c test_gsl.c 
gcc -o test_gsl test_gsl.o linearosc_gvf.o -lgsl -lgslcblas -lm

echo "Running test_gsl"
./test_gsl 
