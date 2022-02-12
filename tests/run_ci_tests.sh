echo "------------------------------"
echo "Testing help"
echo "------------------------------"
vfgen help help
if [ $? -ne 0 ]
then
    exit -1
fi

echo "------------------------------"
echo "Testing octave"
echo "------------------------------"
cd test_octave
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing dde_solver"
echo "------------------------------"
cd test_dde_solver
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing lsoda"
echo "------------------------------"
cd test_lsoda
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing gsl"
echo "------------------------------"
cd test_gsl
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing evf"
echo "------------------------------"
cd test_evf
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing boostodeint"
echo "------------------------------"
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

echo "------------------------------"
echo "Testing scilab"
echo "------------------------------"
cd test_scilab
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing scipy"
echo "------------------------------"
cd test_scipy
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

echo "------------------------------"
echo "Testing pygsl"
echo "------------------------------"
cd test_pygsl
bash run_test.sh
if [ $? -ne 0 ]
then
    exit -1
fi
cd ..

# echo "------------------------------"
# echo "Testing R"
# echo "------------------------------"
# cd test_r
# bash run_test.sh
# if [ $? -ne 0 ]
# then
#     exit -1
# fi
# cd ..
