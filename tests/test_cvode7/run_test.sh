set -e

cd test1

echo "Running vfgen with linearosc.vf"
vfgen cvode7 ../../vf/linearosc.vf

echo "Building the test program test_linearosc_cv7"
mkdir build
cd build
cmake ..
make

echo "Running test_linearosc_cv7"
./test_linearosc_cv7

echo "Cleaning up"
cd ..
rm -rf build
rm linearosc_cv7.c linearosc_cv7.h
cd ..

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cd test2

echo "Running vfgen with linearoscp.vf"
vfgen cvode7:func=yes ../../vf/linearoscp.vf

echo "Building the test program test_linearoscp_cv7"
mkdir build
cd build
cmake ..
make

echo "Running test_linearoscp_cv7"
./test_linearoscp_cv7

echo "Cleaning up"
cd ..
rm -rf build
rm linearoscp_cv7.c linearoscp_cv7.h
cd ..
