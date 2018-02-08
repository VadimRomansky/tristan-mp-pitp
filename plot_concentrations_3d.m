clear;
np = hdf5read('./output/flds.tot.005','densi');
ne = hdf5read('./output/flds.tot.005','dens');

Nx = size(np, 1);
Ny = size(np, 2);

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 20;
mp = 1.67262*10^-24;
me = mp/mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;

omega = 4*pi*n*q*q/me;

rho = c0/(omega*Nskinlength);
c1=0.45;

tau = c1*rho/c0;

densityFactor = 1.0/(rho*rho*rho);

figure(1);
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, np*densityFactor);
shading interp;
title ('np');
xlabel ('y');
ylabel ('x');
zlabel ('np');
grid ;

figure(2);
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, ne*densityFactor);
shading interp;
title ('ne');
xlabel ('y');
ylabel ('x');
zlabel ('ne');
grid ;