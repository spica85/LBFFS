#!/bin/bash
cp input*.txt ../../build/input.txt
cp boundaryConditions*.txt ../../build/boundaryConditions.txt
cp simpleCar.stl ../../build/
cd ../../build/
rm -rf walls.stl movingWalls.stl inlets.stl outlets.stl
ln -s simpleCar.stl walls.stl
./ns.opencl | tee log
