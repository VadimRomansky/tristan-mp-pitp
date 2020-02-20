clear;

directory_name = './output/';
file_name = 'spect';
file_number = '.000';
full_name = strcat(directory_name, file_name, file_number);
fileinfo = h5info(full_name);
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
ntristan = 8; %0.5*ppc0 maybe
sigma = 4.0;
gamma = 1.5;
v = c*sqrt(1 - 1/(gamma*gamma));
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));

samplingFactor = 5;

Nx = size(Bx,1);
Ny = size(Bx,2);
Np = size(g, 1);
startx = 10;
endx = Nx;
fieldStartx = startx;
fieldEndx = endx;
magneticDensity = 0;
dense(1:Nx*samplingFactor) = 0;
densp(1:Nx*samplingFactor) = 0;
for i = fieldStartx:fieldEndx,
    for j = 1:Ny,
        magneticDensity = magneticDensity + (Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j))*fieldScale*fieldScale/(8*3.14);
    end;
end;
magneticDensity = magneticDensity/((fieldEndx - fieldStartx + 1)*Ny);

electronDensity = 0;
for i = startx*samplingFactor:endx*samplingFactor,
    for j = 2:Np,
        electronDensity = electronDensity + me*c*c*(n/ntristan)*fe(i,j)*(g(j)+1)*(g(j) - g(j-1))/(Ny*samplingFactor);
        dense(i) = dense(i) + fe(i,j)/(Ny*samplingFactor);
    end;
end;
electronDensity = electronDensity/(endx*samplingFactor - startx*samplingFactor + 1);

protonDensity = 0;
for i = startx*samplingFactor:endx*samplingFactor,
    for j = 2:Np,
        protonDensity = protonDensity + mp*c*c*(n/ntristan)*fp(i,j)*(g(j)+1)*(g(j) - g(j-1))/(Ny*samplingFactor);
        densp(i) = densp(i) + fe(i,j)*(g(j) - g(j-1))/(Ny*samplingFactor);
    end;
end;
protonDensity = protonDensity/(endx*samplingFactor - startx*samplingFactor + 1);

ramPressure = n*gamma*(me + mp)*v*v;

energyDensity = n*mp*c*c;

epsilonB = magneticDensity/ramPressure;
epsilonE = electronDensity/ramPressure;
epsilonP = protonDensity/ramPressure;

protonEnergyFrac = protonDensity/energyDensity;

newsigma = magneticDensity/energyDensity;