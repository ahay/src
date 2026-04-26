#!/bin/bash
# Example: inspect, modify, and recover RSF data.
set -euo pipefail
mkdir -p /tmp/working-with-rsf-data-demo
cd /tmp/working-with-rsf-data-demo

echo "=== 1. create a 2D RSF ==="
sfspike n1=100 n2=5 k1=50 mag=1 > data.rsf

echo "=== 2. header is plain text ==="
cat data.rsf

echo "=== 3. sfin tells you the shape ==="
sfin data.rsf

echo "=== 4. sfattr summarizes data values ==="
sfattr < data.rsf

echo "=== 5. modify header metadata (add units and labels) ==="
< data.rsf sfput label1="Time" unit1="s" label2="Trace" unit2="" > labeled.rsf
sfin labeled.rsf

echo "=== 6. where is the binary? ==="
grep "^	in=" data.rsf
# Binary lives under $DATAPATH. If you did 'rm *.rsf' you'd orphan it.

echo "=== 7. clean up the RIGHT way ==="
sfrm data.rsf labeled.rsf
ls *.rsf 2>/dev/null && echo "leaked headers" || echo "clean"
