set -e

echo "Testing scipy"
vfgen scipy ../vf/linearosc.vf
python3 test_scipy.py
