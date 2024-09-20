set -e

echo "Generating for linearosc.vf"
vfgen adolc ../vf/linearoscp.vf

echo "Building test_linearosc"
g++ -c linearoscp_adolc.cpp
g++ -c test_linearoscp.cpp
g++ -o test_linearoscp test_linearoscp.o linearoscp_adolc.o -ladolc -lm

echo "Running test_linearosc"
./test_linearoscp
