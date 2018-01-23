clear;
np = hdf5read('./output/flds.tot.001','densi');
ne = hdf5read('./output/flds.tot.001','dens');

Nx = size(np, 1);
Ny = size(np, 2);

figure(1);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, np);
shading interp;
title ('np');
xlabel ('y');
ylabel ('x');
zlabel ('np');
grid ;

figure(2);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, ne);
shading interp;
title ('ne');
xlabel ('y');
ylabel ('x');
zlabel ('ne');
grid ;