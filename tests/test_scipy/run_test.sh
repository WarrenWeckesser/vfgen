set -e

echo "Installing scipy"
apt-get -yqq install python3-numpy python3-scipy

echo "Testing scipy"
vfgen scipy ../vf/linearosc.vf
python3 test_scipy.py
