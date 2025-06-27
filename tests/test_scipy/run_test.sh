set -e

echo "Testing scipy"
vfgen scipy ../vf/linearosc.vf
python3 test_scipy.py

echo "Testing scip with tfirst"
vfgen scipy:tfirst=yes ../vf/linearosc.vf
python3 test_scipy_tfirst.py