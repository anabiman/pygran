// Setup hopper dims
h = 3;
s = 0.3;
l = 1.0;

Point(1) = {0, 0, h, 1.0};
Point(2) = {-l, 0, h, 1.0};
Point(3) = {l, 0, h, 1.0};
Point(4) = {0, -l, h, 1.0};
Point(5) = {0, l, h, 1.0};

Circle(1) = {2, 1, 5};
Circle(2) = {5, 1, 3};
Circle(3) = {3, 1, 4};
Circle(4) = {4, 1, 2};

Line Loop(5) = {1, 2, 3, 4};

Point(6) = {0, 0, 0, 1.0};
Point(7) = {-s, 0, 0, 1.0};
Point(8) = {s, 0, 0, 1.0};
Point(9) = {0, -s, 0, 1.0};
Point(10) = {0, s, 0, 1.0};

Circle(7) = {7, 6, 10};
Circle(8) = {7, 6, 9};
Circle(9) = {9, 6, 8};
Circle(10) = {8, 6, 10};
Line Loop(11) = {8, 9, 10, -7};

Line(13) = {2, 7};
Line(14) = {10, 5};
Line(15) = {8, 3};
Line(16) = {9, 4};
Line Loop(17) = {14, 2, -15, 10};
Ruled Surface(18) = {17};
Line Loop(19) = {15, 3, -16, 9};
Ruled Surface(20) = {19};
Line Loop(21) = {16, 4, 13, 8};
Ruled Surface(22) = {21};
Line Loop(23) = {14, -1, 13, 7};
Ruled Surface(24) = {23};
