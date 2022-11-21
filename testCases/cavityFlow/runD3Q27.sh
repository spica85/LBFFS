#!/bin/bash
runDir=../../D3Q27_opencl/build
cp input*.txt $runDir/input.txt
cp boundaryConditions*.txt $runDir/boundaryConditions.txt
cd $runDir
rm -rf walls.stl movingWalls.stl inlets.stl outlets.stl
./ns.opencl | tee log
