set -e

echo "Running vfgen"
vfgen r ../vf/linearosc.vf

echo "Running R --version"
R --version

echo "Running R test code"
R -f test_r.R --vanilla --slave
