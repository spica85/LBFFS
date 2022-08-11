# LBFFS [English/[Japanese](README_ja.md)]
LBFFS is a flow simulator based on lattice boltzmann method, written in C++ and OpenCL languages to run on a single GPU

## Features
* D3Q19 discrete velocity model
* Recursive regularized collision model
* Large eddy simulation based on Smagorinsky sub-grid scale model with damping near walls
* On-lattice boundary condition on wall by half-way Bounce-Back method
* Off-lattice boundary condition on wall by Filippova & Hanel’s Interpolated Bounce-Back method
* Off-lattice boundary setting by importing a STL file
* Outlet boundary condition which suppresses wave reflections [Geier et al., Comput. Math. Appl. (2015), Appendix F]
* Spongezones which suppress wave reflections at outlet boundaries

## Run script on Google Colaboratory  
[Sample](runScriptOnColab.ipynb)

## Test cases
* Poiseuille flow (Re=100, Laminar)
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

* Lid driven cavity flow (Re=100, Laminar)
<table>
<tr>
<td>Velocity vector</td>
<td>u profile at x=0.5m</td>
</tr>
<tr>
<td><img src="https://user-images.githubusercontent.com/109857341/180638527-6905b752-ebff-4695-a5c2-aacec47b16ac.png" width="320px"></td>
<td><img src="https://user-images.githubusercontent.com/109857341/182006837-1e144cbc-5c16-4bc1-aefe-eba9ca35f386.png" width="320px"></td>
</tr>
</table>

* Flow around a cylinder (Re=100, Laminar)

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

* Flow around a moving cylinder (Laminar)

https://user-images.githubusercontent.com/109857341/184156620-08c4a85d-8176-46f2-8cad-d64d2287c577.mp4

* Backward facing step flow (Re=5500, Turbulence)

https://user-images.githubusercontent.com/109857341/180644458-212d29d3-9d87-4b73-b8d2-fd086c1d4b44.mp4

<img src="https://user-images.githubusercontent.com/109857341/180644496-94171507-2454-4ed6-b495-355ab656610b.png" width="640px">


* Flow around a car (Actual physical properties, Turbulence)

https://user-images.githubusercontent.com/109857341/180644599-89a6945f-214d-449f-8b31-b7e3ab75fc98.mp4

## Licence
[BSD-3-Clause license](LICENSE)


