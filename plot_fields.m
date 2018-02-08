clear;
Bx = hdf5read('./output/flds.tot.005','bx');
By = hdf5read('./output/flds.tot.005','by');
Bz = hdf5read('./output/flds.tot.005','bz');
Ex = hdf5read('./output/flds.tot.005','ex');
Ey = hdf5read('./output/flds.tot.005','ey');
Ez = hdf5read('./output/flds.tot.005','ez');

Nx = size(Bx, 1);
Ny = size(By, 2);

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

fieldFactor = me*rho/(q*tau*tau);

figure(1);
plot ((1:Nx)*rho,Bx(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
plot ((1:Nx)*rho,By(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
plot ((1:Nx)*rho,Bz(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
plot ((1:Nx)*rho,Ex(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
plot ((1:Nx)*rho,Ey(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
plot ((1:Nx)*rho,Ez(1:Nx, fix(Ny/2))*fieldFactor, 'red');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;