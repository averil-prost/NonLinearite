//+
Point(1) = {0, 0, 0, 0.15};
//+
Point(2) = {2, 0, 0, 0.1};
//+
Point(3) = {2, 1, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.15};
//+
Point(5) = {0.2, 0.6, 0, 0.12};
//+
Point(6) = {0.3, 0.7, 0, 0.12};
//+
Point(7) = {0.4, 0.7, 0, 0.09};
//+
Point(8) = {0.5, 0.6, 0, 0.08};
//+
Point(9) = {0.5, 0.5, 0, 0.08};
//+
Point(10) = {0.4, 0.4, 0, 0.09};
//+
Point(11) = {0.3, 0.4, 0, 0.12};
//+
Point(12) = {0.2, 0.5, 0, 0.12};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 5};
//+
Line(12) = {5, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {12, 5, 6, 7, 8, 9, 10, 11};
//+
Plane Surface(1) = {1, 2};
//+
Field[1] = Box;
//+
Field[1].Thickness = 0;
//+
Field[1].XMax = 1.5;
//+
Field[1].XMin = 0.5;
//+
Field[1].YMax = 0.75;
//+
Field[1].YMin = 0.25;
//+
Field[1].VIn = 0.09;
//+
Field[1].VOut = 0.12;
//+
Field[1].VIn = 0.01;
//+
Field[1].VIn = 0.09;
//+
Background Field = 1;
