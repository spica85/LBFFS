# LBFFS
LBFFS is a flow simulator based on Lattice Boltzmann Method, written in C++ and OpenCL languages to run on a single GPU

## Features
* D3Q19 discrete velocity model
* Recursive regularized collision model
* Large eddy simulation based on Smagorinsky sub-grid scale model with damping near walls
* On-lattice boundary condition on wall by half-way Bounce-Back method
* Off-lattice boundary condition on wall by Filippova & Hanel’s Interpolated Bounce-Back method
* Off-lattice boundary setting by importing a STL file
* Outlet boundary condition which suppresses wave reflections [Geier et al., Comput. Math. Appl. (2015), Appendix F]
* Spongezones which suppress wave reflections at outlet boundaries

## Test cases
* Poiseuille flow (Re=100)
<table>
<tr>
<td>Velocity distribution</td>
<td>u profile</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180640617-7e83c0b4-61df-4ed4-ac4f-39554b86affe.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/180640633-b6779f8d-1921-493f-b64f-876f08a873d8.png" width="320px"></td>
</tr>
</table>

* Lid driven cavity flow (Re=100)
<table>
<tr>
<td>Velocity vector</td>
<td>u profile at x=0.5m</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180638527-6905b752-ebff-4695-a5c2-aacec47b16ac.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/180638616-66064f75-6159-4655-a28d-1c7f0be1dcc7.png" width="320px"></td>
</tr>
</table>

* Flow around a cylinder

https://user-images.githubusercontent.com/109857341/180644337-b0e62fda-41a7-487d-9cee-98e37b96f939.mp4
<table>
<tr>
<td>Cx</td>
<td>Cy</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180644297-db37a9b0-177e-4a2a-8390-9ada3c0c96dd.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/180644308-5bb19345-ec3a-4d10-8ad9-4452e091d878.png" width="320px"></td>
</tr>
</table>

* Backward facing step flow
* Flow around a car
