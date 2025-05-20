! ! WARNING ! !                 Only results in 2D are available               ! ! WARNING ! ! 

Python code conceived for simulating acosutic waves propagation in uniform base flows with the particularity that monopolar source can change its position in time.
It discretizes the Linearized Euler Equations using the Discontinuous Galerkin Method in time.

Some functions are adapted from the MATLAB scripts that accompany the text on Nodal Disontinuous Galerkin methods by Jan S. Hesthaven and Tim Warburton (https://github.com/tcew/nodal-dg/tree/master) (Publisher web-page: http://www.springer.com/mathematics/computational+science+%26+engineering/book/978-0-387-72065-4)

VERSIONS
-------------------------
Python 3.11.11
GMSH 4.13.1 (24 May 2024) 

MESH
-------------------------
Mesh generated using GMSH. Documentation and .exe available at https://gmsh.info/ .
There is 2D meshes as well as a 3D mesh in this repository


SIMULATIONS
-------------------------
Lauch the simulation: $ python .\main.py
Parameters can be defined in 'parameters.py'
Results are displayed using the plt.scatter() function in Python
