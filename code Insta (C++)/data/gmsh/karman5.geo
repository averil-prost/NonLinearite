// Gmsh project created on Fri Dec  3 19:48:19 2021
//+
Point(1) = {0, 0, 0, 0.15};
//+
Point(2) = {0, 1, 0, 0.15};
//+
Point(3) = {2.0, 1, 0, 0.15};
//+
Point(4) = {2.0, 0, 0, 0.15};
//+
Point(5) = {0.4, 0.4, 0, 0.06};
//+
Point(6) = {0.4, 0.6, 0, 0.06};
//+
Point(7) = {0.35, 0.5, 0, 0.05};
//+
Point(8) = {0.25, 0.5, 0, 0.08};
//+
Line(1) = {8, 6};
//+
Line(2) = {6, 7};
//+
Line(3) = {7, 5};
//+
Line(4) = {5, 8};
//+
Line(5) = {2, 1};
//+
Line(6) = {1, 4};
//+
Line(7) = {4, 3};
//+
Line(8) = {3, 2};
//+
Curve Loop(1) = {5, 6, 7, 8};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1, 2};
//+
Field[1] = Box;
Field[1].VIn = 0.07;
Field[1].VOut = 0.2;
Field[1].XMin = 0.4;
Field[1].XMax = 1.5;
Field[1].YMin = 0.25;
Field[1].YMax = 0.75;
Field[1].Thickness = 0.2;
//+
Background Field = 1;
//+
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm = 5;
