
echo "Installing R and the desolve package"
sudo add-apt-repository ppa:marutter/c2d4u
sudo apt-get update
apt-get -y install r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
