set -e

echo "Generating matlab code to be tested with octave"
vfgen matlab ../vf/linearosc.vf

echo "Running test_octave script with functions generated for matlab"
octave --no-gui --no-window-system test_matlab_with_octave.m
