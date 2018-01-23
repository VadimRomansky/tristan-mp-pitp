clear;
Bx = hdf5read('./output/flds.tot.001','bx');
By = hdf5read('./output/flds.tot.001','by');
Bz = hdf5read('./output/flds.tot.001','bz');
Ex = hdf5read('./output/flds.tot.001','ex');
Ey = hdf5read('./output/flds.tot.001','ey');
Ez = hdf5read('./output/flds.tot.001','ez');

Nx = size(Bx, 1);
Ny = size(By, 2);

figure(1);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, Bx);
shading interp;
title ('Bx');
xlabel ('y');
ylabel ('x');
zlabel ('Bx');
grid ;

figure(2);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('y');
ylabel ('x');
zlabel ('By');
grid ;

figure(3);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, Bz);
shading interp;
title ('Bz');
xlabel ('y');
ylabel ('x');
zlabel ('Bz');
grid ;

figure(4);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, Ex);
shading interp;
title ('Ex');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
grid ;

figure(5);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, Ey);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(6);
[X, Y] = meshgrid(1:Ny, 1:Nx);
surf(X, Y, Ez);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;