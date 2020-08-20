clear;

directory_name = './output2/';
file_name = 'spect';
file_number = '.010';
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
Te = 4.2E11;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
ntristan = 2; %0.5*ppc0 maybe
sigma = 0.4;
gamma = 1.5;
v = c*sqrt(1 - 1/(gamma*gamma));
beta = v/c;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));

samplingFactor = 20;

gmax = 200;

Nx = size(Bx,1);
Ny = size(Bx,2);
Np = size(g, 1);
startx = 5000;
endx = 10000;
fieldStartx = startx/samplingFactor;
fieldEndx = endx/samplingFactor;
magneticDensity = 0;
dense(1:Nx*samplingFactor) = 0;
densp(1:Nx*samplingFactor) = 0;
for i = fieldStartx:fieldEndx,
    for j = 1:Ny,
        magneticDensity = magneticDensity + (Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j))*fieldScale*fieldScale/(8*3.14);
    end;
end;
magneticDensity = magneticDensity/((fieldEndx - fieldStartx + 1)*Ny);

Fejuttner(1:Np) = 0;
for i = 1:Np,
    exp1 = exp(-(1 + g(i))/thetae);
    beta = sqrt(1 - 1.0/((1 + g(i))^2));
    bes = besselk(2, 1/thetae);
    Fejuttner(i) = (((1+g(i))^2)*beta/(thetae*bes))*exp1;
end;


electronEnergyDensity = 0;
electronHighEnergyDensity = 0;
electronHighEnergyDensity2 = 0;
electronDensity = 0;
for i = startx:endx,
    for j = 1:Np,
        electronEnergyDensity = electronEnergyDensity + me*c*c*(n/ntristan)*fe(i,j)*g(j)*(g(j)+1)/(Ny*samplingFactor);
        electronDensity = electronDensity + (n/ntristan)*fe(i,j)*g(j)/(Ny*samplingFactor);
        if(g(j) > gmax)
            electronHighEnergyDensity = electronHighEnergyDensity + me*c*c*(n/ntristan)*fe(i,j)*g(j)*(g(j)+1)/(Ny*samplingFactor);
        end;
        dense(i) = dense(i) + fe(i,j)*g(j)/(Ny*samplingFactor);
    end;
end;
electronEnergyDensity = electronEnergyDensity/(endx - startx + 1);
electronDensity = electronDensity/(endx - startx + 1);
electronHighEnergyDensity = electronHighEnergyDensity/(endx - startx + 1);

for i = startx:endx,
    for j = 2:Np,
        if(g(j)>gmax)
            a = (n/ntristan)*fe(i,j)*g(j)/(Ny*samplingFactor);
            b = electronDensity*Fejuttner(j)*(g(j) - g(j-1));
            electronHighEnergyDensity2 = electronHighEnergyDensity2 + me*c*c*((n/ntristan)*fe(i,j)*g(j)/(Ny*samplingFactor) - electronDensity*Fejuttner(j)*(g(j) - g(j-1)))*(g(j)+1);
        end;
    end;
end;
electronHighEnergyDensity2 = me*c*c*electronHighEnergyDensity2/(endx - startx + 1);

electronMeanGamma = electronEnergyDensity/(me*c*c*electronDensity);

protonDensity = 0;
for i = startx:endx,
    for j = 1:Np,
        protonDensity = protonDensity + mp*c*c*(n/ntristan)*fp(i,j)*g(j)*(g(j)+1)/(Ny*samplingFactor);
        densp(i) = densp(i) + fp(i,j)*g(j)/(Ny*samplingFactor);
    end;
end;
protonDensity = protonDensity/(endx - startx + 1);

ramPressure = n*gamma*(me + mp)*v*v;

energyDensity = n*gamma*mp*c*c;

epsilonB = magneticDensity/ramPressure;
epsilonE = electronHighEnergyDensity/ramPressure;
epsilonE2 = electronHighEnergyDensity2/ramPressure;
epsilonEtot = electronEnergyDensity/ramPressure;
epsilonP = protonDensity/ramPressure;

epsilonECrumley = (electronHighEnergyDensity/electronEnergyDensity)*me*(electronMeanGamma - 1)/(mp*(gamma - 1));

protonEnergyFrac = protonDensity/energyDensity;

newsigma = magneticDensity/energyDensity;