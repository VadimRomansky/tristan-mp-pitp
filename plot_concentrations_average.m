clear;
directory_name = './output/';
file_name = 'flds.tot';
number = 2;
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');


Nx = size(np, 1);
Ny = size(np, 2);

npa(1:Nx)=0;
nea(1:Nx)=0;

for i=1:Nx,
    for j=1:Ny,
        npa(i)=npa(i) + np(i,j)/Ny;
        nea(i)=nea(i) + (ne(i,j) - np(i,j))/Ny;
    end;
end;

mp = 1.67262177E-24;
q = 4.84*10^-10;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 0.004;
gamma = 1.5;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));
samplingFactor = 20;
timestep = 10000;

omega = sqrt(4*3.14*n*q*q/(gamma*(me*mp)/(me + mp)));
dt = ctristan/(comp*omega);
dx = samplingFactor*c*dt/ctristan;
dt = dt*timestep;

t = dt*number;
x = 0;

for i = Nx/2:Nx,
    if(npa(Nx - i + 1) > 2*ntristan)
        x = (Nx - i + 1)*dx;
        break;
    end;
end;


figure(1);
plot ((1:Nx),npa(1:Nx)/ntristan, 'red');
title ('np');
xlabel ('x');
ylabel ('np');
grid ;

figure(2);
plot ((1:Nx),nea(1:Nx)/ntristan, 'red');
title ('ne');
xlabel ('x');
ylabel ('ne');
grid ;

%dlmwrite('np.dat',np,'delimiter',' ');
%dlmwrite('ne.dat',ne,'delimiter',' ');
