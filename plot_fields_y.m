clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.000';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');

Nx = size(Bx, 1);
Ny = size(By, 2);

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

tau = c1*rho/c0;

fieldFactor = me*rho/(q*tau*tau);

samplingFactor = 10;

rho = 0.2;
rho = rho*samplingFactor;

xpoint = fix(Nx/2) + 1;

figure(1);
plot ((1:Ny)*rho,Bx(xpoint, 1:Ny)*fieldFactor, 'red');
title ('Bx');
xlabel ('y');
ylabel ('Bx');
grid ;

figure(2);
plot ((1:Ny)*rho,By(xpoint, 1:Ny)*fieldFactor, 'red');
title ('By');
xlabel ('y');
ylabel ('By');
grid ;

figure(3);
plot ((1:Ny)*rho,Bz(xpoint, 1:Ny)*fieldFactor, 'red');
title ('Bz');
xlabel ('y');
ylabel ('Bz');
grid ;

figure(4);
plot ((1:Ny)*rho,Ex(xpoint, 1:Ny)*fieldFactor, 'red');
title ('Ex');
xlabel ('y');
ylabel ('Ex');
grid ;

figure(5);
plot ((1:Ny)*rho,Ey(xpoint, 1:Ny)*fieldFactor, 'red');
title ('Ey');
xlabel ('y');
ylabel ('Ey');
grid ;

figure(6);
plot ((1:Ny)*rho,Ez(xpoint, 1:Ny)*fieldFactor, 'red');
title ('Ez');
xlabel ('y');
ylabel ('Ez');
grid ;