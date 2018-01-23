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
plot (1:Nx,Bx(1:Nx, Ny/2), 'red');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
plot (1:Nx,By(1:Nx, Ny/2), 'red');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
plot (1:Nx,Bz(1:Nx, Ny/2), 'red');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
plot (1:Nx,Ex(1:Nx, Ny/2), 'red');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
plot (1:Nx,Ey(1:Nx, Ny/2), 'red');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
plot (1:Nx,Ez(1:Nx, Ny/2), 'red');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;