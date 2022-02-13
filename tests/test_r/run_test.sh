set -e

echo "Installing R and the desolve package"
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
apt-get update
apt-get -y install r-base r-cran-desolve

echo "Testing R"
vfgen r ../vf/linearosc.vf

R -f test_r.R --vanilla --slave
