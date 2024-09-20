set -e

echo "Generating for linearoscp.vf"
vfgen adolc ../vf/linearoscp.vf

echo "Building test_linearoscp"
g++ -c linearoscp_adolc.cpp
g++ -c test_linearoscp.cpp
g++ -o test_linearoscp test_linearoscp.o linearoscp_adolc.o -ladolc -lm

echo "Running test_linearoscp"
./test_linearoscp
