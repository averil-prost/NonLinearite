Point(1) = {0, 0, 0, 0.05};
Point(2) = {1, 0, 0, 0.05};
Point(3) = {1, 1, 0, 0.05};
Point(4) = {0, 1, 0, 0.05};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {3, 4, 1, 2};
Surface(1) = {1};
