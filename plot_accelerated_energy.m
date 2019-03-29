clear;
directory_name = './output/';
file_name = 'spect';
density_name = 'flds.tot';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

full_name = strcat(directory_name, density_name, file_number);
fileinfo = hdf5info(full_name);
np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');

Nx = size(fp,1);
Np = size(fp,2);

Nxn = size(np,1);
Nyn = size(np, 2);


startx = 1;
endx = Nx/4;

Fp(1:5,1:Nx)=0;
Fe(1:5,1:Nx)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

me = 0.91*10^-27;
mass_ratio = 16;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

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

samplingFactor = Nx/Nxn;

rho = rho*samplingFactor;
rho = 1/Nskinlength;
rho2 = samplingFactor/Nskinlength;

normp = 0;
norme = 0;
norm = 1;

for i = 1:Np,
    %Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    %Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    Pp(i) = sqrt((g(i)+1)^2 - 1);
    Pe(i) = sqrt((g(i)+1)^2 - 1);
    for j = startx:endx,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i)=Fp(i);
    Fe(i)=Fe(i);
    
end;

normp = Fp(1)*(g(2) - g(1));
norme = Fe(1)*(g(2) - g(1));

for i = 2:Np,
    normp = normp + Fp(i)*(g(i) - g(i-1));
    norme = norme + Fe(i)*(g(i) - g(i-1));
end;

startAccelerationE = 175;
startAccelerationP = 145;
acceleratedProtonEnergyFraction(1:Nx) = 0;
acceleratedElectronEnergyFraction(1:Nx) = 0;
acceleratedProtonFraction(1:Nx)=0;
acceleratedElectronFraction(1:Nx)=0;
totalEnergyP(1:Nx) = 0;
totalEnergyE(1:Nx) = 0;
totalEnergy = 0;
concentrationP(1:Nx) = 0;
concentrationE(1:Nx) = 0;

averageConcentrationP(1:Nxn,2) = 0;
averageConcentrationE(1:Nxn,2) = 0;
for i = 1:Nx,
    acceleratedEnergyP = 0;
    acceleratedEnergyE = 0;
    acceleratedP = 0;
    acceleratedE = 0;
    for j = 2:Np,
        concentrationP(i) = concentrationP(i) + fp(i,j)*(g(j) - g(j-1));
        concentrationE(i) = concentrationE(i) + fe(i,j)*(g(j) - g(j-1));
        totalEnergyP(i) = totalEnergyP(i) +  mass_ratio*fp(i,j)*(1 + g(j))*(g(j) - g(j-1));
        totalEnergyE(i) = totalEnergyE(i) +  fe(i,j)*(1 + g(j))*(g(j) - g(j-1));
        if(j > startAccelerationP)
            acceleratedEnergyP = acceleratedEnergyP + mass_ratio*fp(i,j)*(1 + g(j))*(g(j) - g(j-1));
            acceleratedP = acceleratedP + fp(i,j)*(g(j) - g(j-1));
        end;
        if(j > startAccelerationE)
            acceleratedEnergyE = acceleratedEnergyE + fe(i,j)*(1 + g(j))*(g(j) - g(j-1));
            acceleratedE = acceleratedE + fe(i,j)*(g(j) - g(j-1));
        end;
    end;
    acceleratedProtonEnergyFraction(i) = acceleratedEnergyP/totalEnergyP(i);
    acceleratedElectronEnergyFraction(i) = acceleratedEnergyE/totalEnergyE(i);
    
    acceleratedProtonFraction(i) = acceleratedP/concentrationP(i);
    acceleratedElectronFraction(i) = acceleratedE/concentrationE(i);
    
    totalEnergy = totalEnergy + totalEnergyP(i) + totalEnergyE(i);
end;

for i=1:Nxn,
    averageConcentrationP(i,1) = concentrationP((i-1)*samplingFactor + 1);
    averageConcentrationE(i,1) = concentrationE((i-1)*samplingFactor + 1);
    for j = 1:Nyn,
        averageConcentrationP(i,2) = averageConcentrationP(i,2) + np(i,j);
        averageConcentrationE(i,2) = averageConcentrationE(i,2) + (ne(i,j) - np(i,j));
    end;
end;

figure(1);
plot ((1:Nx)*rho,acceleratedProtonEnergyFraction(1:Nx), 'red');
title ('accelerated proton energy fraction');
xlabel ('x');
ylabel ('f');
grid ;

figure(2);
plot ((1:Nx)*rho,acceleratedElectronEnergyFraction(1:Nx), 'red');
title ('accelerated electron energy fraction');
xlabel ('x');
ylabel ('f');
grid ;

figure(3);
plot ((1:Nx)*rho,totalEnergyP(1:Nx), 'red');
title ('proton energy');
xlabel ('x');
ylabel ('f');
grid ;

figure(4);
plot ((1:Nx)*rho,totalEnergyE(1:Nx), 'red');
title ('electron energy');
xlabel ('x');
ylabel ('f');
grid ;

figure(5);
plot ((1:Nxn)*rho2,averageConcentrationP(1:Nxn,1), 'red', (1:Nxn)*rho2, averageConcentrationP(1:Nxn,2),'blue');
title ('proton n');
xlabel ('x');
ylabel ('f');
grid ;

figure(6);
plot ((1:Nxn)*rho2,averageConcentrationE(1:Nxn,1), 'red', (1:Nxn)*rho2, averageConcentrationE(1:Nxn,2),'blue');
title ('electron n');
xlabel ('x');
ylabel ('f');
grid ;

figure(7);
plot ((1:Nx)*rho,acceleratedProtonFraction(1:Nx), 'red');
title ('accelerated proton fraction');
xlabel ('x');
ylabel ('f');
grid ;

figure(8);
plot ((1:Nx)*rho,acceleratedElectronFraction(1:Nx), 'red');
title ('accelerated electron  fraction');
xlabel ('x');
ylabel ('f');
grid ;
