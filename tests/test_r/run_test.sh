set -e

echo "Running vfgen"
vfgen r ../vf/linearosc.vf
vfgen r ../vf/delayed_logistic.vf

echo "Running R --version"
R --version

echo "Running R test ode code"
R -f test_r_ode.R --vanilla --slave

echo "Running R test dede code"
R -f test_r_dede.R --vanilla --slave
