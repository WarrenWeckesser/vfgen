set -e

echo "Installing swig, setuptools, wheel, numpy and pygsl"
apt-get -yqq install swig
python3 -m pip install setuptools wheel
python3 -m venv venv
source venv/bin/activate
python3 -m pip install numpy==1.26.4
git clone https://github.com/pygsl/pygsl.git
cd pygsl
pip3 install .
cd ..

echo "Testing pygsl"
vfgen pygsl ../vf/linearosc.vf
python3 test_pygsl.py

deactivate
