clear;
directory_name = './output2/';
file_name = 'flds';
file_number = '.tot.015';
Nd = 5;
start = 0;

Color = {'red','blue','green','black','magenta'};
%LegendTitle = {'2','16','32','128'};
LegendTitle = {'B normal','B quasiparallel', 'anisotropic turbulence bz', 'isotropic turbulence','anisotropic turbulence by'};


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

samplingFactor = 5;

rho = rho*samplingFactor;
rho = samplingFactor/Nskinlength;



Np(1:Nd,1:Nx)=0;
Ne(1:Nd,1:Nx)=0;

ypoint = fix(Ny/4);
for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    Np0 = hdf5read(full_name,'densi');
    Ne0 = hdf5read(full_name,'dens');
    for i = 1:Nx,
        Np(j,i)=Np0(i, ypoint);
        Ne(j,i)=Ne0(i, ypoint)- Np0(i, ypoint);
    end;
end;

N0 = Np(1,fix(Nx-1));

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Np(j, 1:Nx)/N0,'color', Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},LegendTitle{5},'Location','southeast');
title ('Np');
xlabel ('x');
ylabel ('Np');
grid ;

figure(2);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Ne(j, 1:Nx)/N0,'color', Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},LegendTitle{5},'Location','southeast');
title ('Ne');
xlabel ('x');
ylabel ('Ne');
grid ;