n = 40;
radius = 6;
depth = 1;

[X,Y,Z] = cylinder(radius, n);
DT = [];

for i = 1 : 2*(n+1) - 3
   DT = [DT; i, i+1, i+2; i+1 i+2 i+3];
end

n = n + 1;

X = reshape(X, n*2, 1);
Y = reshape(Y, n*2, 1);
Z = reshape(Z, n*2, 1);

n = n - 1;

figure(1)
hold on
shift = 1.0;
ncyl = 20;
C = Z * 0;

trisCoords = [];

for i = 1 : ncyl
    
    if i == ncyl/2
       r = radius - depth;
       x = min(X) : 0.5 * (max(X) - min(X)) / n : max(X);
       y = min(Y) : 0.5 * (max(Y) - min(Y)) / n : max(Y);
       
       [x,y] = meshgrid(x, y);
       nr = length(x);
       
       x = reshape(x, nr^2, 1);
       y = reshape(y, nr^2, 1);
       
       xr = x(x.^2 + y.^2 <= r^2);
       yr = y(x.^2 + y.^2 <= r^2);
       zr = ones(length(xr),1) *  (Z(i) + i * shift);
       
       dt = delaunay(xr, yr);
       
       trimesh(dt, xr, yr, zr, xr .* 0);
       
       trisCoords = [trisCoords; [xr, yr, zr]];
    end
    
    trimesh(DT, X,Y,Z + i*shift, C);
    trisCoords = [trisCoords; [X, Y, Z + i*shift]];
    
end

axis equal;
xlabel('x'); ylabel('y'); zlabel('z');

for i = 1 : 5
    CX = [CX, X];
    CY = [CY, Y];
    CZ = [CZ, Z * i];
end
