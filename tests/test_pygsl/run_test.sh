set -e

echo "Installing swig, setuptools, wheel, numpy and pygsl"
apt-get -yqq install swig
python3 -m pip install setuptools wheel
python3 -m pip install numpy
python3 -m pip install pygsl

echo "Testing pygsl"
vfgen pygsl ../vf/linearosc.vf
python3 test_pygsl.py
