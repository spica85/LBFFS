#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cp backStepBoxes.stl ../../build/
cd ../../build/
rm -rf walls.stl movingWalls.stl inlets.stl outlets.stl
ln -s backStepBoxes.stl walls.stl
./ns.opencl | tee log
