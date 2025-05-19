//Include "parameters.geo";
xmin_ph = -1;
xmax_ph = 1;
ymin_ph = -1;
ymax_ph = 1;
zmin_ph = 2;
zmax_ph = 4;
buf = 1.0;
Np = 4.0;
Np_ph = Np;

Point(1) = {xmin_ph, ymin_ph, zmin_ph, Np};
Point(2) = {xmax_ph, ymin_ph, zmin_ph, Np};
Point(3) = {xmax_ph, ymax_ph, zmin_ph, Np};
Point(4) = {xmin_ph, ymax_ph, zmin_ph, Np};

Point(5) = {xmin_ph, ymin_ph, zmax_ph, Np};
Point(6) = {xmax_ph, ymin_ph, zmax_ph, Np};
Point(7) = {xmax_ph, ymax_ph, zmax_ph, Np};
Point(8) = {xmin_ph, ymax_ph, zmax_ph, Np};


//Point(11) = {xmin_ph-buf, ymin_ph-buf, zmin_ph-buf, Np_ph};
//Point(12) = {xmax_ph+buf, ymin_ph-buf, zmin_ph-buf, Np_ph};
//Point(13) = {xmax_ph+buf, ymax_ph+buf, zmin_ph-buf, Np_ph};
//Point(14) = {xmin_ph-buf, ymax_ph+buf, zmin_ph-buf, Np_ph};

//Point(15) = {xmin_ph-buf, ymin_ph-buf, zmax_ph+buf, Np_ph};
//Point(16) = {xmax_ph+buf, ymin_ph-buf, zmax_ph+buf, Np_ph};
//Point(17) = {xmax_ph+buf, ymax_ph+buf, zmax_ph+buf, Np_ph};
//Point(18) = {xmin_ph-buf, ymax_ph+buf, zmax_ph+buf, Np_ph};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9) = {2, 6};
Line(10) = {7, 3};
Line(11) = {4, 8};
Line(12) = {5, 1};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Line Loop(3) = {1, 9, -5, 12};
Line Loop(4) = {2, -10, -6, -9};
Line Loop(5) = {3, 11, -7, 10};
Line Loop(6) = {4, -12, -8, -11};


Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
//Plane Surface(2) = {2, -1};

// Physical Surface(0) = {1};
//Physical Surface(1) = {2};
//Physical Line(2) = {5, 6, 7, 8};

Surface Loop(25) = {1, -2, -3, -4, -5, -6};
Volume(26) = {25};

Physical Surface(25) = {1};
Physical Volume(26) = {26};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;
Mesh.SecondOrderLinear = 1;
// Mesh.Smoothing = 20;
