clear;
directory_name = './output1/';
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
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');

Nx = size(fp,1);
Np = size(fp,2);

Nxn = size(np,1);
Nyn = size(np, 2);


startx = 1;
endx = Nx/4;


Pp(1:Np)=0;
Pe(1:Np)=0;
Fp(1:Np)=0;
Fe(1:Np)=0;

me = 4.67307679E-03;
mass_ratio = 64;
mp = me*mass_ratio;
c = 0.45;
c0 = 0.45;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

g0 = 1.5;
beta0 = sqrt(1 - 1/(g0*g0));

Nskinlength = 10;

q = 4.67307679E-03;
n = 4;

omega = sqrt(4*pi*n*q*q/me);

rho = c0/(omega*Nskinlength);
c1=0.45;

tau = c1*rho/c0;

fieldFactor = me*rho/(q*tau*tau);

samplingFactor = fix(Nx/Nxn);
if(Nxn > Nx)
    samplingFactor = 1;
    Nxn = Nx;
end;

rho = rho*samplingFactor;
rho = 1/Nskinlength;
rho2 = samplingFactor/Nskinlength;

normp = 0;
norme = 0;
norm = 1;

for i = 1:Np,
    Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    %Pp(i) = sqrt((g(i)+1)^2 - 1);
    %Pe(i) = sqrt((g(i)+1)^2 - 1);
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


acceleratedEnergyP = 0;
acceleratedP = 0;
acceleratedEnergyE = 0;
acceleratedE = 0;

concentrationP = 0;
concentrationE = 0;

acceleratedProtonEnergyFraction = 0;
acceleratedElectronEnergyFraction = 0;
acceleratedProtonFraction=0;
acceleratedElectronFraction=0;
totalEnergyP = 0;
totalEnergyE = 0;
totalEnergy = 0;

fieldEnergy = 0;
fieldEnergyFraction = 0;

mag_e_init = 0.02;
mag_e_0 = 0.0002375;

for i = 20:endx/samplingFactor,
    for j = 1:Nyn,
        fieldEnergy = fieldEnergy + ((Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j)))*samplingFactor/Nyn;
    end;
end;

fieldEnergy = fieldEnergy*mag_e_init/mag_e_0;

meangp = 0;
meange = 0;

for j = 2:Np,
    totalEnergyP = totalEnergyP +  mp*c*c*Fp(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
    %totalEnergyP = totalEnergyP +  mp*c*c*Fp(j)*(g(j))*(g(j) - g(j-1));
    meangp = meangp + Fp(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
    meange = meange + Fe(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
    totalEnergyE = totalEnergyE +  me*c*c*Fe(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
    %totalEnergyE = totalEnergyE +  me*c*c*Fe(j)*(g(j))*(g(j) - g(j-1));
    concentrationP = concentrationP + Fp(j)*(g(j) - g(j-1))/(4*pi);
    concentrationE = concentrationE + Fe(j)*(g(j) - g(j-1))/(4*pi);
    if(Pp(j) > 2*g0*beta0*mp*c)
        acceleratedEnergyP = acceleratedEnergyP + mp*c*c*Fp(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
        %acceleratedEnergyP = acceleratedEnergyP + mp*c*c*Fp(j)*(g(j))*(g(j) - g(j-1));
        acceleratedP = acceleratedP + Fp(j)*(g(j) - g(j-1))/(4*pi);
    end;
    if(Pe(j) > 2*g0*beta0*mp*c)
        acceleratedEnergyE = acceleratedEnergyE + me*c*c*Fe(j)*(1 + g(j))*(g(j) - g(j-1))/(4*pi);
        %acceleratedEnergyE = acceleratedEnergyE + me*c*c*Fe(j)*(g(j))*(g(j) - g(j-1));
        acceleratedE = acceleratedE + Fe(j)*(g(j) - g(j-1))/(4*pi);
    end;
end;

meangp = meangp/concentrationP;
meange = meange/concentrationE;
    
    
acceleratedProtonFraction = acceleratedP/concentrationP;
acceleratedElectronFraction = acceleratedE/concentrationE;
    
totalEnergy = fieldEnergy + totalEnergyP + totalEnergyE;
    
acceleratedProtonEnergyFraction = acceleratedEnergyP/totalEnergy;
%acceleratedProtonEnergyFraction = (acceleratedEnergyP/totalEnergyP)*(mp*(meangp - 1))/(mp*(g0 - 1));
acceleratedElectronEnergyFraction = acceleratedEnergyE/totalEnergy;
%acceleratedElectronEnergyFraction = (acceleratedEnergyE/totalEnergyE)*(me*(meange - 1))/(mp*(g0 - 1));
fieldEnergyFraction = fieldEnergy/(totalEnergyP + totalEnergyE);



