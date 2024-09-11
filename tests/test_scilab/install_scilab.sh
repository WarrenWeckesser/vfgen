set -e

echo "Installing scilab"
apt-get -yqq install scilab-cli
scilab-cli -version

