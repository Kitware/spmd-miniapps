#! /bin/bash

if [ $# -ne 2 ]; then
  echo "Incorrect arguments"
  echo "Usage: "
  echo "$0 <bin_path> <log_postfix>"
  exit 1
fi

SPATH=`dirname $0`
BPATH=$1
DPATH=$HOME/Devs/sujin/data
POSTFIX=$2

echo "1] Marching Cubes, scalar_2.1"
$SPATH/strongscaling.sh \
  "$BPATH/MarchingCubes $DPATH/large.image.vtk /dev/null 1.3 scalar_2.1" 4 \
  | tee "mc-scalar_2.1-$POSTFIX" | $SPATH/gentable.py
echo ""

echo "2] Marching Cubes, simd_2.1"
$SPATH/strongscaling.sh \
  "$BPATH/MarchingCubes $DPATH/large.image.vtk /dev/null 1.3 simd_2.1" 4 \
  | tee "mc-simd_2.1-$POSTFIX" | $SPATH/gentable.py
echo ""

echo "3] Marching Tetrahedra, scalar"
$SPATH/strongscaling.sh \
  "$BPATH/TetmeshContour $DPATH/large.tet.vtk /dev/null 1.3 scalar" 4 \
  | tee "mt-scalar-$POSTFIX" | $SPATH/gentable.py
echo ""

echo "4] Marching Tetrahedra, simd"
$SPATH/strongscaling.sh \
  "$BPATH/TetmeshContour $DPATH/large.tet.vtk /dev/null 1.3 simd" 4 \
  | tee "mt-simd-$POSTFIX" | $SPATH/gentable.py
echo ""

echo "5] Marching Tetrahedra, simd_2"
$SPATH/strongscaling.sh \
  "$BPATH/TetmeshContour $DPATH/large.tet.vtk /dev/null 1.3 simd_2" 4 \
  | tee "mt-simd_2-$POSTFIX" | $SPATH/gentable.py
echo ""

