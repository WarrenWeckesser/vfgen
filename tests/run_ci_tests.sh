echo "Testing help"
vfgen help help
if [ $? -ne 0 ]
then
    exit -1
fi

echo "Testing octave"
cd test_octave
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing dde_solver"
cd test_dde_solver
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing lsoda"
cd test_lsoda
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing gsl"
cd test_gsl
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing evf"
cd test_evf
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing boostodeint"
cd test_boostodeint
apt-get -y install libboost-math-dev
bash run_test_boostodeint_default_noparams.sh
if [ $? -ne 0 ]
then
    exit -1
fi
bash run_test_boostodeint_default_params.sh
if [ $? -ne 0 ]
then
    exit -1
fi
bash run_test_boostodeint_implicit_noparams.sh
if [ $? -ne 0 ]
then
    exit -1
fi
bash run_test_boostodeint_implicit_params.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing scilab"
cd test_scilab
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "Testing scipy"
cd test_scipy
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

# echo "Testing R"
# cd test_r
# bash run_test.sh
# if [ $? -ne 0 ]
# then
#     exit -1
# fi
# cd ..
