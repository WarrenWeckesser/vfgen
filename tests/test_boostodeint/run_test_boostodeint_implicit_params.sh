set -e

echo "Running vfgen."
vfgen boostodeint:system=implicit ../vf/linearoscp.vf

echo "Building test_boostodeint_implicit_params."
g++ -c linearoscp_vf.cpp
g++ -c test_boostodeint_implicit_params.cpp
g++ test_boostodeint_implicit_params.o linearoscp_vf.o -o test_boostodeint_implicit_params

echo "Running test_boostodeint_implicit_params."
./test_boostodeint_implicit_params
