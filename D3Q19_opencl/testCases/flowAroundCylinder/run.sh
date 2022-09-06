#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cp cylinder.stl ../../build/
cd ../../build/
rm -rf walls.stl
rm -rf movingWalls.stl
ln -s cylinder.stl walls.stl
./ns.opencl | tee log
