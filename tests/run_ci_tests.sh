set -e

echo "------------------------------"
echo "Testing help"
echo "------------------------------"
vfgen help help


echo "------------------------------"
echo "Testing R"
echo "------------------------------"
cd test_r
bash run_test.sh
cd ..
