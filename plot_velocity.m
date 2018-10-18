clear;
directory_name = './output1/';
file_name = 'flds1.tot';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fileinfo = hdf5info(full_name);
Vp = hdf5read(full_name,'v4xi');
Ve = hdf5read(full_name,'v4x');

Nx = size(Vp, 1);
Ny = size(Vp, 2);

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
samplingFactor = 5;
tau = c1*rho/c0;
rho =0.1;

%densityFactor = 1.0/(rho*rho*rho);
rho = rho*samplingFactor;

figure(1);
plot ((1:Nx)*rho,Vp(1:Nx, fix(Ny/2)+1), 'red');
title ('Vp');
xlabel ('x');
ylabel ('Vp');
grid ;

figure(2);
plot ((1:Nx)*rho,Ve(1:Nx, fix(Ny/2)+1), 'red');
title ('Ve');
xlabel ('x');
ylabel ('Ve');
grid ;
