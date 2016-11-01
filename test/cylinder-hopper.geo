lc = 0.5;
h = 20;
hs = 10;
r = 6;
rs = 4;

Point(1) = {-r, 0, 0, lc}; 
Point(2) = {0, -r, 0, lc}; 
Point(3) = {r, 0, 0, lc}; 
Point(4) = {0, r, 0, lc}; 
Point(5) = {0, 0, 0, lc}; 

Point(6) = {-rs, 0, hs, lc}; 
Point(7) = {0, -rs, hs, lc}; 
Point(8) = {rs, 0, hs, lc}; 
Point(9) = {0, rs, hs, lc};
Point(10) = {0, 0, hs, lc};

Circle(1) = {1, 5, 2}; 
Circle(2) = {2, 5, 3}; 
Circle(3) = {3, 5, 4}; 
Circle(4) = {4, 5, 1}; 

Circle(5) = {6, 10, 7}; 
Circle(6) = {7, 10, 8}; 
Circle(7) = {8, 10, 9}; 
Circle(8) = {9, 10, 6}; 

Line Loop(9) = {1, 2, 3, 4}; 
geoExtrude[] = Extrude {0, 0, h} { Line{1,2,3,4}; }; 

Line Loop(10) = {5, 6, 7, 8};
Plane Surface(6) = {10};