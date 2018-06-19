clear;
directory_name = './output1/';
file_name = 'flds';
file_number = '.tot.003';
Nd = 2;
start = 4;

full_name = strcat(directory_name, file_name, num2str(start), file_number);
Np0 = hdf5read(full_name,'densi');
Ne0 = hdf5read(full_name,'dens');
Nx = size(Np0, 1);
Ny = size(Np0, 2);

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
densityFactor = 1.0/(rho*rho*rho);

samplingFactor = 1;

rho = rho*samplingFactor;


Np(1:Nd,1:Nx)=0;
Ne(1:Nd,1:Nx)=0;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    Np0 = hdf5read(full_name,'densi');
    Ne0 = hdf5read(full_name,'dens');
    for i = 1:Nx,
        Np(j,i)=Np0(i, fix(Ny/2)+1);
        Ne(j,i)=Ne0(i, fix(Ny/2)+1);
    end;
end;

Color = {'red','blue'};

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Np(j, 1:Nx)*densityFactor,'color', Color{j});
end;
legend('smooth','no smooth','Location','southeast');
title ('Np');
xlabel ('x');
ylabel ('Np');
grid ;

figure(2);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Ne(j, 1:Nx)*densityFactor,'color', Color{j});
end;
legend('smooth','no smooth','Location','southeast');
title ('Ne');
xlabel ('x');
ylabel ('Ne');
grid ;