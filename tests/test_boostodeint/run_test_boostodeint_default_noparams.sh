set -e

echo "Running vfgen."
vfgen boostodeint ../vf/linearosc.vf

echo "Building test_boostodeint_default_noparams."
g++ -c linearosc_vf.cpp
g++ -c test_boostodeint_default_noparams.cpp
g++ test_boostodeint_default_noparams.o linearosc_vf.o -o test_boostodeint_default_noparams

echo "Running test_boostodeint_default_noparams."
./test_boostodeint_default_noparams
