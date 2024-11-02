// Gmsh project created on Fri Sep 27 18:03:25 2024
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 0.5, 0.5, 2*Pi};
//+
Cylinder(2) = {0, 0, 0, 0, 0, 0.25, 0.5, 2*Pi};
//+
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Physical Volume("Volume_Bottom", 1) = {2};
//+
Physical Volume("Volume_Top", 2) = {3};
//+
Physical Surface("interface", 3) = {2};
//+
Physical Surface("surface_bottom", 4) = {5};
//+
Physical Surface("surface_top", 5) = {3};
//+
Physical Surface("surface_side_bottom", 6) = {4};
//+
Physical Surface("surface_side_top", 7) = {1};
//+
Physical Curve("contact_line", 1) = {2};
// Interface
MeshSize {1} = 0.2;
// Bottom and Top Face
MeshSize {3, 2} = 0.2;