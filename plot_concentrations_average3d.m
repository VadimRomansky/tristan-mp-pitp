clear;
directory_name = './output6/';
file_name = 'flds.tot';
file_number = '.002';
full_name = strcat(directory_name, file_name, file_number);
np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');


Nx = size(np, 1);
Ny = size(np, 2);
Nz = size(np, 3);

npa(1:Nx)=0;
nea(1:Nx)=0;

for i=1:Nx,
    for j=1:Ny,
        for k = 1:Nz,
        npa(i)=npa(i) + np(i,j,k)/(Ny*Nz);
        nea(i)=nea(i) + (ne(i,j,k) - np(i,j,k))/(Ny*Nz);
        end,
    end;
end;

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
plot ((1:Nx),npa(1:Nx), 'red');
title ('np');
xlabel ('x');
ylabel ('np');
grid ;

figure(2);
plot ((1:Nx),nea(1:Nx), 'red');
title ('ne');
xlabel ('x');
ylabel ('ne');
grid ;

dlmwrite('np.dat',np,'delimiter',' ');
dlmwrite('ne.dat',ne,'delimiter',' ');
