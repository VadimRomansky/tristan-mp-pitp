clear;
np = hdf5read('./output/flds.tot.001','densi');
ne = hdf5read('./output/flds.tot.001','dens');

Nx = size(np, 1);
Ny = size(np, 2);

figure(1);
plot (1:Nx,np(1:Nx, Ny/2), 'red');
title ('np');
xlabel ('x');
ylabel ('np');
grid ;

figure(2);
plot (1:Nx,ne(1:Nx, Ny/2), 'red');
title ('ne');
xlabel ('x');
ylabel ('ne');
grid ;
