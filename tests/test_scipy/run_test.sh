
echo "Installing scipy"
apt-get -yqq install python-scipy

echo "Testing scipy"
vfgen scipy ../vf/linearosc.vf
python test_scipy.py
