clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

Nx = size(fp,1);
Np = size(fp,2);

Fp(1:Np)=0;
Fe(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

me = 0.91*10^-27;
mass_ratio = 100;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

normp = 0;
norme = 0;
norm = 1;

Nx1 = 2000;
Nx2 = 4000;

for i = 1:Np,
    Pp(i) = g(i);
    Pe(i) = g(i);
    for j = Nx1:Nx2,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i) = Fp(i)*(g(i));
    Fe(i) = Fe(i)*(g(i));
end;

figure(1);
plot (Pp(1:Np),Fp(1:Np), 'red');
title ('F_p');
xlabel ('gamma-1');
ylabel ('Fp');
grid ;

figure(2);
plot (Pe(1:Np),Fe(1:Np), 'red');
title ('F_e');
xlabel ('gamma-1');
ylabel ('F_e');
grid ;