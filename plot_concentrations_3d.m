clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.025';
full_name = strcat(directory_name, file_name, file_number);
np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');

Nx = size(np, 1);
Ny = size(np, 2);

offsetx = 10;
offsety = 0;

ne1(1:Nx-2*offsetx,1:Ny-2*offsety)=0;
np1(1:Nx-2*offsetx,1:Ny-2*offsety)=0;

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


for i=1:Nx-2*offsetx,
    for j=1:Ny-2*offsety,
        ne1(i,j)=ne(i+offsetx,j+offsety);
        np1(i,j)=ne(i+offsetx,j+offsety);
    end;
end;

densityFactor = 1.0/(rho*rho*rho);
set(0,'DefaultFigureColormap',feval('jet'));
figure(1);
[X, Y] = meshgrid((1+offsety:Ny-offsety)*rho, (1+offsetx:Nx-offsetx)*rho);
surf(X, Y, np1*densityFactor);
shading interp;
title ('np');
xlabel ('y');
ylabel ('x');
zlabel ('np');
grid ;

figure(2);
[X, Y] = meshgrid((1+offsety:Ny-offsety)*rho, (1+offsetx:Nx-offsetx)*rho);
surf(X, Y, ne1*densityFactor);
shading interp;
title ('ne');
xlabel ('y');
ylabel ('x');
zlabel ('ne');
grid ;