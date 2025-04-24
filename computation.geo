Include "parameters.geo";
xmin_ph = -1;
xmax_ph = 1;
ymin_ph = -1;
ymax_ph = 1;
buf = 1.0;
Np = 2.0;
Np_ph = Np;

Point(1) = {xmin_ph, ymin_ph, 0, Np};
Point(2) = {xmax_ph, ymin_ph, 0, Np};
Point(3) = {xmax_ph, ymax_ph, 0, Np};
Point(4) = {xmin_ph, ymax_ph, 0, Np};

Point(5) = {xmin_ph-buf, ymin_ph-buf, 0, Np_ph};
Point(6) = {xmax_ph+buf, ymin_ph-buf, 0, Np_ph};
Point(7) = {xmax_ph+buf, ymax_ph+buf, 0, Np_ph};
Point(8) = {xmin_ph-buf, ymax_ph+buf, 0, Np_ph};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};


Plane Surface(1) = {1};
Plane Surface(2) = {2, -1};

Physical Surface(0) = {1};
Physical Surface(1) = {2};
Physical Line(2) = {5, 6, 7, 8};


Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;
Mesh.SecondOrderLinear = 1;
// Mesh.Smoothing = 20;
