clear;
directory_name = './output5/';
file_name = 'flds.tot';
file_number = '.002';
full_name = strcat(directory_name, file_name, file_number);
np2 = hdf5read(full_name,'densi');
ne2 = hdf5read(full_name,'dens');
file_number = '.004';
full_name = strcat(directory_name, file_name, file_number);
np4 = hdf5read(full_name,'densi');
ne4 = hdf5read(full_name,'dens');
file_number = '.006';
full_name = strcat(directory_name, file_name, file_number);
np6 = hdf5read(full_name,'densi');
ne6 = hdf5read(full_name,'dens');
file_number = '.008';
full_name = strcat(directory_name, file_name, file_number);
np8 = hdf5read(full_name,'densi');
ne8 = hdf5read(full_name,'dens');
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
np10 = hdf5read(full_name,'densi');
ne10 = hdf5read(full_name,'dens');

Nx = size(np2, 1);
Ny = size(np2, 2);

npa(1:Nx,1:5)=0;
nea(1:Nx,1:5)=0;

for i=1:Nx,
    for j=1:Ny,
        npa(i,1)=npa(i,1) + np2(i,j)/Ny;
        nea(i,1)=nea(i,1) + (ne2(i,j) - np2(i,j))/Ny;
        npa(i,2)=npa(i,2) + np4(i,j)/Ny;
        nea(i,2)=nea(i,2) + (ne4(i,j) - np4(i,j))/Ny;
        npa(i,3)=npa(i,3) + np6(i,j)/Ny;
        nea(i,3)=nea(i,3) + (ne6(i,j) - np6(i,j))/Ny;
        npa(i,4)=npa(i,4) + np8(i,j)/Ny;
        nea(i,4)=nea(i,4) + (ne8(i,j) - np8(i,j))/Ny;
        npa(i,5)=npa(i,5) + np10(i,j)/Ny;
        nea(i,5)=nea(i,5) + (ne10(i,j) - np10(i,j))/Ny;
    end;
end;

mp = 1.67262177E-24;
q = 4.84*10^-10;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 0.04;
gamma = 1.5;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));
samplingFactor = 20;
timestep = 10000;

omega = sqrt(4*3.14*n*q*q/(gamma*me));
dt = ctristan/(comp*omega);
dx = samplingFactor*c*dt/ctristan;
dt = dt*timestep;

figure(1);
plot ((1:Nx)*dx,npa(1:Nx,1)/ntristan, 'magenta',(1:Nx)*dx,npa(1:Nx,2)/ntristan, 'black',(1:Nx)*dx,npa(1:Nx,3)/ntristan, 'green',(1:Nx)*dx,npa(1:Nx,4)/ntristan, 'blue',(1:Nx)*dx,npa(1:Nx,5)/ntristan, 'red');
title ('density');
xlabel ('x cm');
ylabel ('\rho/{\rho_0}');
legend('t{\omega}_p = 1800','t{\omega}_p = 3600','t{\omega}_p = 5400','t{\omega}_p = 7200','t{\omega}_p = 9000')
grid ;

