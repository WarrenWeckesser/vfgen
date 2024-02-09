set -e

echo "Fetching LSODA from netlib."
curl -s -S -O "https://www.netlib.org/odepack/opkdmain.f"
curl -s -S -O "https://www.netlib.org/odepack/opkda1.f"
curl -s -S -O "https://www.netlib.org/odepack/opkda2.f"
