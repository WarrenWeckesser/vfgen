
echo "Installing octave"
apt-get -yqq install octave

echo "Testing octave"
cd test_octave
vfgen octave ../vf/linearosc.vf
octave --no-gui --no-window-system test_octave.m
