//+
Point(1) = {0, 0, 0, 0.15};
//+
Point(2) = {4, 0, 0, 0.15};
//+
Point(3) = {4, 1, 0, 0.15};
//+
Point(4) = {0, 1, 0, 0.15};
//+
Point(5) = {0.4, 0.6, 0, 0.15};
//+
Point(6) = {0.5, 0.7, 0, 0.15};
//+
Point(7) = {0.6, 0.7, 0, 0.15};
//+
Point(8) = {0.7, 0.6, 0, 0.15};
//+
Point(9) = {0.7, 0.5, 0, 0.15};
//+
Point(10) = {0.6, 0.4, 0, 0.15};
//+
Point(11) = {0.5, 0.4, 0, 0.15};
//+
Point(12) = {0.4, 0.5, 0, 0.15};
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
