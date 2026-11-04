#!/bin/sh
# Run the model on THIS board and package everything needed for comparison.
#
#   sh collect.sh            -> writes ./bms_result_<hostname>/
#   sh collect.sh DIR        -> writes ./DIR/
#
# Send that folder back. Nothing is installed and nothing outside the folder
# is touched (apart from the compiled binary in the current directory).

set -e

HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:-"bms_result_$(hostname 2>/dev/null || uname -n)"}

if command -v cc >/dev/null 2>&1; then
    CC=cc
elif command -v gcc >/dev/null 2>&1; then
    CC=gcc
elif command -v clang >/dev/null 2>&1; then
    CC=clang
else
    echo "No C compiler found."
    echo "  Debian/Ubuntu : sudo apt install build-essential"
    echo "  Fedora        : sudo dnf install gcc"
    exit 1
fi

echo "compiler : $CC"
"$CC" -std=c99 -O2 -o bms_gaussian "$HERE/bms_gaussian.c" -lm
echo "built    : $(pwd)/bms_gaussian"

mkdir -p "$OUT"

# Three plain runs. The plain run is the timing measurement: no CSV writing,
# no frame stream, identical workload to every other board.
echo
echo "plain runs (timing)"
for i in 1 2 3; do
    ./bms_gaussian > "$OUT/run_$i.txt" 2>&1
    grep -h "wall    :" "$OUT/run_$i.txt" | sed "s/^/  run $i  /"
done

# One export run, which is what the heatmaps are drawn from. Its timing is a
# different workload and is not used for comparisons.
echo
echo "export run (heatmaps)"
./bms_gaussian --export "$OUT/data" > "$OUT/run_export.txt" 2>&1
grep -h "wall    :" "$OUT/run_export.txt" | sed 's/^/  export /'

# Everything needed to explain the numbers: what the CPU is, how it is being
# clocked, and how the compiler was invoked.
{
    echo "== board =="
    uname -a
    echo
    echo "== cpu =="
    nproc 2>/dev/null || echo "nproc unavailable"
    grep -m4 -E "model name|Hardware|Processor|CPU part" /proc/cpuinfo 2>/dev/null || true
    echo
    echo "== cpu frequency =="
    cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo "no cpufreq"
    cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq 2>/dev/null || true
    cat /sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq 2>/dev/null || true
    echo
    echo "== memory =="
    grep -E "^MemTotal" /proc/meminfo 2>/dev/null || true
    echo
    echo "== compiler =="
    "$CC" --version | head -1
    echo
    echo "== binary =="
    ls -l bms_gaussian
    echo
    echo "== load at start =="
    cat /proc/loadavg 2>/dev/null || true
} > "$OUT/environment.txt" 2>&1

cp "$HERE/bms_gaussian.c" "$OUT/bms_gaussian.c"
cp "$HERE/collect.sh" "$OUT/collect.sh" 2>/dev/null || true

echo
echo "done -> $OUT"
echo "send back the whole folder."
