
echo "Installing R and the desolve package"
apt-get -yqq install lsb_release
deb https://ppa.launchpadcontent.net/marutter/c2d4u/ubuntu $(lsb_release -c -s) main
deb-src https://ppa.launchpadcontent.net/marutter/c2d4u/ubuntu $(lsb_release -c -s) main
sudo add-apt-repository ppa:marutter/c2d4u
sudo apt-get update
apt-get -y install r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
