#! /bin/bash

if [ $# -ne 2 ]; then
  echo "Incorrect arguments"
  echo "Usage: "
  echo "$0 \"command\" <ntrys>"
  exit 1
fi

CMDLINE=$1
NTRYS=$2

ncpus=`grep "physical id" /proc/cpuinfo | sort -u | wc -l`
corespercpu=`grep "core id" /proc/cpuinfo | sort -u | wc -l`
let ncores=ncpus*corespercpu
nthreads=`grep "processor" /proc/cpuinfo | wc -l`
let threadspercore=nthreads/ncores

echo "========================================================================"
echo "Performing strong scaling tests on [$CMDLINE]"
echo "Running each case $NTRYS times"
echo "System info:"
echo "    hostname: `hostname`"
echo "    ncpus: $ncpus"
echo "    cores per cpu: $corespercpu"
echo "    ncores: $ncores"
echo "    nthreads: $nthreads"
echo "    threads per core: $threadspercore"
echo "========================================================================"

nthreads_list="1 2 4 8 12 16 24 32 40 48 56 64 80 96 112 128 160 192 224 256"

function run_test {
  export OMP_NUM_THREADS=$1
  echo "Running with $1 threads"

  echo "----------------------------------------------------------------------"
  for i in `seq 1 $2` ;
  do
    $CMDLINE
    echo ""
  done
  echo "----------------------------------------------------------------------"
}

# run tests on number of threads from nthreads_list
for i in $nthreads_list ;
do
  if [ "$i" -gt $ncores ]; then
    break
  fi
  last_val=$i

  run_test $i $NTRYS
done

# run tests on full machine if ncores is not in nthreads_list
if [ "$last_val" -ne $ncores ]; then
  run_test $ncores $NTRYS
fi

echo ""
echo "=== SMP performance Tests ==="

# run tests on HT threads
let next=ncores*2
for i in `seq $next $ncores $nthreads` ;
do
  run_test $i $NTRYS
done

