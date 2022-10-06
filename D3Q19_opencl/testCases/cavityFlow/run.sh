#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cd ../../build/
rm -rf walls.stl movingWalls.stl inlets.stl outlets.stl
./ns.opencl | tee log
