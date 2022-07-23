# LBFFS
LBFFS is a flow simulator based on Lattice Boltzmann Method, written in C++ and OpenCL languages to run on a single GPU.

## Features
* D3Q19 discrete velocity model
* Recursive regularized collision model
* Large eddy simulation based on Smagorinsky sub-grid scale model with damping near walls
* On-lattice boundary condition on wall by half-way Bounce-Back method
* Off-lattice boundary condition on wall by Filippova & Hanel’s Interpolated Bounce-Back method
* Off-lattice boundary setting by importing a STL file
* Outlet boundary condition which suppresses wave reflections [Geier et al., Comput. Math. Appl. (2015), Appendix F]
* Spongezone which suppress wave reflections at outlet boundary condition

## Test cases
* Lid driven cavity flow
* Poiseuille flow
* Flow around a cylinder
* Backward facing step flow
* Flow around a car
