set -e

echo "------------------------------"
echo "Testing help"
echo "------------------------------"
vfgen help help

echo "------------------------------"
echo "Testing octave"
echo "------------------------------"
cd test_octave
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing dde_solver"
echo "------------------------------"
cd test_dde_solver
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing lsoda"
echo "------------------------------"
cd test_lsoda
bash download_lsoda.sh
bash run_test_linearosc.sh
bash run_test_pidecay.sh
cd ..

echo "------------------------------"
echo "Testing gsl"
echo "------------------------------"
cd test_gsl
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing evf"
echo "------------------------------"
cd test_evf
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing boostodeint"
echo "------------------------------"
cd test_boostodeint
apt-get -y install libboost-math-dev
bash run_test_boostodeint_default_noparams.sh
bash run_test_boostodeint_default_params.sh
bash run_test_boostodeint_implicit_noparams.sh
bash run_test_boostodeint_implicit_params.sh
cd ..

echo "------------------------------"
echo "Testing scilab"
echo "------------------------------"
cd test_scilab
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing scipy"
echo "------------------------------"
cd test_scipy
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing pygsl"
echo "------------------------------"
cd test_pygsl
bash run_test.sh
cd ..

echo "------------------------------"
echo "Testing R"
echo "------------------------------"
cd test_r
bash run_test.sh
cd ..
