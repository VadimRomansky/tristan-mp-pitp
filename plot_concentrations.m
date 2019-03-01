clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fileinfo = hdf5info(full_name);
np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');

Nx = size(np, 1);
Ny = size(np, 2);

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 20;
mp = 1.67262*10^-24;
me = mp/mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;

omega = sqrt(4*pi*n*q*q/me);

rho = c0/(omega*Nskinlength);
c1=0.45;
samplingFactor = 10;
tau = c1*rho/c0;
rho =0.1;

%densityFactor = 1.0/(rho*rho*rho);
rho = rho*samplingFactor;
densityFactor = 1/8;

figure(1);
plot ((1:Nx)*rho,np(1:Nx, fix(Ny/2)+1)*densityFactor, 'red');
title ('np');
xlabel ('x');
ylabel ('np');
grid ;

figure(2);
plot ((1:Nx)*rho,(ne(1:Nx, fix(Ny/2)+1)-np(1:Nx, fix(Ny/2)+1))*densityFactor, 'red');
title ('ne');
xlabel ('x');
ylabel ('ne');
grid ;
