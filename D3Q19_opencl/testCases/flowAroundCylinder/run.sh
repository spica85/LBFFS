#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cp cylinder.stl outlets.stl ../../build/
cd ../../build/
rm -rf walls.stl movingWalls.stl inlets.stl
ln -s cylinder.stl walls.stl
./ns.opencl | tee log
