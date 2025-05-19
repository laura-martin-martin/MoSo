// Gmsh project created on Mon Apr 28 13:12:46 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {-1, -1, -1, 2, 2, 2};//+
MeshSize {:} = 2;
Physical Volume(25) = {1};