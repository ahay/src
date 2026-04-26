#!/bin/bash
# Example: build a synthetic trace, filter it, window it, transpose, and summarize.
# Run with:
#   bash example-pipeline.sh
set -euo pipefail
mkdir -p /tmp/using-sf-programs-demo
cd /tmp/using-sf-programs-demo

echo "=== stage 1: synthesize a 2D dataset (20 traces x 1000 samples) ==="
sfspike n1=1000 n2=20 nsp=2 k1=300,700 mag=1,0.5 > spikes.rsf

echo "=== stage 2: bandpass each trace (low-pass at 4 Hz) ==="
< spikes.rsf sfbandpass fhi=4 phase=y > filtered.rsf

echo "=== stage 3: window out the middle 500 samples ==="
< filtered.rsf sfwindow n1=500 f1=250 > windowed.rsf

echo "=== stage 4: transpose (now axis 1 is trace number) ==="
< windowed.rsf sftransp > transposed.rsf

echo "=== stage 5: summarize ==="
< transposed.rsf sfattr

echo
echo "=== same thing as a single pipe ==="
sfspike n1=1000 n2=20 nsp=2 k1=300,700 mag=1,0.5 \
  | sfbandpass fhi=4 phase=y \
  | sfwindow n1=500 f1=250 \
  | sftransp \
  | sfattr
