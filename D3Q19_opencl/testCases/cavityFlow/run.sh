#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cd ../../build/
rm -rf walls.stl
rm -rf movingWalls.stl
./ns.opencl | tee log
