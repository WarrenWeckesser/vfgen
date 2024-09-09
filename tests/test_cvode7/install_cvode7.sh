set -e

mkdir cvode_install
cd cvode_install

echo "Fetching CVODE 7.1.1"
curl -L -s -S -O "https://github.com/LLNL/sundials/releases/download/v7.1.1/cvode-7.1.1.tar.gz"
tar xf cvode-7.1.1.tar.gz

echo "Building CVODE 7.1.1"
cd cvode-7.1.1
mkdir build
cd build
cmake ..
make

echo "Installing CVODE 7.1.1"
sudo make install
