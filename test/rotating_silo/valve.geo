// Rotating valve for silo hopper
hx = 0.5;
hy = 0.25;
hz = 0.5;

Point(1) = {-hx, 0, hz, 1.0};
Point(2) = {hx, 0, hz, 1.0};
Point(3) = {-hx, -hy, hz, 1.0};
Point(4) = {hx, hy, hz, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 1};
//+
Line Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(6) = {5};

Point(5) = {0, -hx, hz, 1.0};
Point(6) = {0, hx, hz, 1.0};
Point(7) = {-hy, -hx, hz, 1.0};
Point(8) = {hy, hx, hz, 1.0};
//+
Line(7) = {8, 5};
//+
Line(8) = {6, 7};
//+
Line(9) = {8, 6};
//+
Line(10) = {7, 5};
//+
Line Loop(11) = {10, -7, 9, 8};
//+
Plane Surface(12) = {11};
