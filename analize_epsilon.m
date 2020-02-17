clear;

directory_name = './output/';
file_name = 'spect';
file_number = '.013';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

file_name = 'flds.tot';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');

mp = 1.67262177E-24;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 4.0;
gamma = 1.5;
v = c*sqrt(1 - 1/(gamma*gamma));
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));

samplingFactor = 20;

Nx = size(Bx,1);
Ny = size(Bx,2);
Np = size(g, 1);
startx = 10000;
endx = 20000;
fieldStartx = startx/samplingFactor;
fieldEndx = endx/samplingFactor;
magneticDensity = 0;
for i = fieldStartx:fieldEndx,
    for j = 1:Ny,
        magneticDensity = magneticDensity + (Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j))*fieldScale*fieldScale/(8*3.14);
    end;
end;
magneticDensity = magneticDensity/((fieldEndx - fieldStartx + 1)*Ny);

electronDensity = 0;
for i = startx:endx,
    for j = 1:Np,
        electronDensity = electronDensity + me*c*c*n*fe(i,j)*g(j)/(Ny*samplingFactor);
    end;
end;
electronDensity = electronDensity/(endx - startx + 1);

ramPressure = n*gamma*(me + mp)*v*v;

epsilonB = magneticDensity/ramPressure;
epsilonE = electronDensity/ramPressure;