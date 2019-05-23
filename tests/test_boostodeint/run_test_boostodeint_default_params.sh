
echo "Running vfgen."
vfgen boostodeint ../vf/linearoscp.vf

echo "Building test_boostodeint_default_params."
g++ -c linearoscp_vf.cpp
g++ -c test_boostodeint_default_params.cpp
g++ test_boostodeint_default_params.o linearoscp_vf.o -o test_boostodeint_default_params

echo "Running test_boostodeint_default_params."
./test_boostodeint_default_params
