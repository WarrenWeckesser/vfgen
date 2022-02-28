set -e

echo "Fetching LSODA from netlib."
curl -s -S -O "http://www.netlib.org/odepack/{opkdmain.f,opkda1.f,opkda2.f}"
