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

ypoint = fix(Ny/2) + 1;

Bnorm(1:Nx) = 0;
Bperp(1:Nx) = 0;
theta(1:Nx) = 0;

for i=1:Nx,
        Bnorm(i) = sqrt(Bx(i,ypoint)*Bx(i,ypoint) + By(i,ypoint)*By(i,ypoint) + Bz(i,ypoint)*Bz(i,ypoint));
        Bperp(i) = sqrt(By(i,ypoint)*By(i,ypoint) + Bz(i,ypoint)*Bz(i,ypoint));
        theta(i) = acos(Bx(i,ypoint)/Bnorm(i))*180/pi;
end;

figure(1);
plot ((1:Nx)*rho,Bx(1:Nx, ypoint)*fieldFactor, 'red');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
plot ((1:Nx)*rho,By(1:Nx, ypoint)*fieldFactor, 'red');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
plot ((1:Nx)*rho,Bz(1:Nx, ypoint)*fieldFactor, 'red');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
plot ((1:Nx)*rho,Ex(1:Nx, ypoint)*fieldFactor, 'red');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
plot ((1:Nx)*rho,Ey(1:Nx, ypoint)*fieldFactor, 'red');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
plot ((1:Nx)*rho,Ez(1:Nx, ypoint)*fieldFactor, 'red');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(7);
plot ((1:Nx)*rho,Bnorm(1:Nx)*fieldFactor, 'red');
title ('Bnorm');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(8);
plot ((1:Nx)*rho,Bperp(1:Nx)*fieldFactor, 'red');
title ('Bperp');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(9);
plot ((1:Nx)*rho,theta(1:Nx), 'red');
title ('theta');
xlabel ('x');
ylabel ('theta');
grid ;