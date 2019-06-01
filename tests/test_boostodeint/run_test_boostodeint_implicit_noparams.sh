
echo "Running vfgen."
vfgen boostodeint:system=implicit ../vf/linearosc.vf

echo "Building test_boostodeint_implicit_noparams."
g++ -c linearosc_vf.cpp
g++ -c test_boostodeint_implicit_noparams.cpp
g++ test_boostodeint_implicit_noparams.o linearosc_vf.o -o test_boostodeint_implicit_noparams

echo "Running test_boostodeint_default_noparams."
./test_boostodeint_implicit_noparams
