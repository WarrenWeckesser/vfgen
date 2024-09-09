set -e

echo "Running vfgen cvode7"
vfgen cvode7 ../vf/linearosc.vf

echo "Building the test program test_linearosc_cv7"
mkdir build
cd build
cmake ..
make

echo "Running test_linearosc_cv7"
./test_linearosc_cv7
