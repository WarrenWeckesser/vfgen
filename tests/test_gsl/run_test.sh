set -e

echo "Generating for linearosc.vf"
vfgen gsl ../vf/linearosc.vf

echo "Building test_linearosc"
gcc -c linearosc_gvf.c
gcc -c test_linearosc.c 
gcc -o test_linearosc test_linearosc.o linearosc_gvf.o -lgsl -lgslcblas -lm

echo "Running test_linearosc"
./test_linearosc 

echo "Generating for linearoscp.vf"
vfgen gsl:func=yes ../vf/linearoscp.vf

echo "Building test_linearoscp"
gcc -c linearoscp_gvf.c
gcc -c test_linearoscp.c 
gcc -o test_linearoscp test_linearoscp.o linearoscp_gvf.o -lgsl -lgslcblas -lm

echo "Running test_linearoscp"
./test_linearoscp
