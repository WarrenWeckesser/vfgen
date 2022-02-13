set -e

echo "Installing R and the desolve package"
apt-get update
apt-get -y install r-base r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
