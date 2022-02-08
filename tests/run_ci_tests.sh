echo "Testing help"
vfgen help help

echo "Testing octave"
cd test_octave
bash run_test.sh
cd ..

echo "Testing dde_solver"
cd test_dde_solver
bash run_test.sh
cd ..

echo "Testing lsoda"
cd test_lsoda
bash run_test.sh
cd ..

echo "Testing gsl"
cd test_gsl
bash run_test.sh
cd ..

echo "Testing evf"
cd test_evf
bash run_test.sh
cd ..

echo "Testing boostodeint"
cd test_boostodeint
apt-get -y install libboost-math-dev
bash run_test_boostodeint_default_noparams.sh
bash run_test_boostodeint_default_params.sh
bash run_test_boostodeint_implicit_noparams.sh
bash run_test_boostodeint_implicit_params.sh
cd ..

echo "Testing scilab"
cd test_scilab
bash run_test.sh
cd ..

echo "Testing scipy"
cd test_scipy
bash run_test.sh
cd ..

echo "Testing R"
cd test_r
bash run_test.sh
cd ..
