clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.015';
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
mass_ratio = 25;
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

for i = 1:Np,
    Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    for j = 1:Nx,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
    
    exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*c*c))/thetae);
    bes = besselk(2, 1/thetae);
    p = Pe(i);
    p3 = (p/(me*c))^3;
    Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
    
    
    
    exp1 = exp(-sqrt(1+Pp(i)*Pp(i)/(mp*mp*c*c))/thetap);
    bes = besselk(2, 1/thetap);
    p = Pp(i);
    p3 = (p/(mp*c))^3;
    Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(i);
end;

normp = (Fp(1)/(Pp(2)^2))*(Pp(2) - Pp(1));
norme = (Fe(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normp = normp + (Fp(i)/(Pp(i)^2))*(Pp(i) - Pp(i-1));
    norme = norme + (Fe(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;

for i = 1:Np,
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
end;

figure(1);
%plot (Pp(1:Np)/(mp*c),Fp(1:Np), 'red');
plot (Pp(1:Np)/(mp*c),Fp(1:Np), 'red',Pp(1:Np)/(mp*c), Fpjuttner(1:Np), 'blue');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
grid ;

figure(2);
%plot (Pe(1:Np)/(me*c),Fe(1:Np), 'red');
plot (Pe(1:Np)/(me*c),Fe(1:Np), 'red',Pe(1:Np)/(me*c), Fejuttner(1:Np), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
grid ;