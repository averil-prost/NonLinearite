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
Point(5) = {0.4, 0.5, 0, 0.05};
//+
Point(6) = {0.2, 0.5, 0, 0.1};
//+
Point(7) = {0.3, 0.4, 0, 0.07};
//+
Point(8) = {0.3, 0.6, 0, 0.07};
//+
Point(9) = {0.37, 0.57, 0, 0.06};
//+
Point(10) = {0.23, 0.57, 0, 0.09};
//+
Point(11) = {0.23, 0.43, 0, 0.09};
//+
Point(12) = {0.37, 0.43, 0, 0.06};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Line(5) = {5, 9};
//+
Line(6) = {9, 8};
//+
Line(7) = {8, 10};
//+
Line(8) = {10, 6};
//+
Line(9) = {6, 11};
//+
Line(10) = {11, 7};
//+
Line(11) = {7, 12};
//+
Line(12) = {12, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8, 9, 10, 11, 12};
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
