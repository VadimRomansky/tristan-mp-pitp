clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.001';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

Nx = size(fp,1);
Np = size(fp,2);

samplingFactor = 20;
startFieldX = 5000;
endFieldX = 5500;

startx = startFieldX*samplingFactor;
endx = endFieldX*samplingFactor;

Fp(1:Np)=0;
Fe(1:Np)=0;
Fekappa(1:Np) = 0;

Fpold(1:Np)=0;
Feold(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;
Fekappa(1:Np) = 0;
Fpkappa(1:Np) = 0;

mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;
gam = 1.5;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^11;
Tp = 2.4*10^12;
Pekappa = 1.5*me*c;
Ppkappa = mp*c;
kappa = 4;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 1.0;
fractionp = 1.0;

Apkappa = ((pi*(kappa - 1.5))^(-1.5))*gamma(kappa + 1)/(gamma(kappa - 0.5)*(Ppkappa/(mp*c))^3);
Aekappa = ((pi*(kappa - 1.5))^(-1.5))*gamma(kappa + 1)/(gamma(kappa - 0.5)*(Pekappa/(me*c))^3);


normp = 0;
norme = 0;
norm = 1;

for i = 2:Np,
    %Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    %Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    Pp(i) = sqrt((g(i)+1)^2 - 1);
    Pe(i) = sqrt((g(i)+1)^2 - 1);
    for j = startx:endx,
        Fp(i) = Fp(i) + fp(j,i)*g(i)/(g(i) - g(i-1));
        Fe(i) = Fe(i) + fe(j,i)*g(i)/(g(i) - g(i-1));
        Fpold(i) = Fpold(i) + fp(j,i);
        Feold(i) = Feold(i) + fe(j,i);
    end;
    Fp(i)=Fp(i)*(1/(gam*beta))*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(me/(gam*beta*mp))*(Pe(i)^3)/(1+g(i));
    Fpold(i)=Fpold(i)*(1/(gam*beta))*(Pp(i)^3)/(1+g(i));
    Feold(i)=Feold(i)*(me/(gam*beta*mp))*(Pe(i)^3)/(1+g(i));
    
    %exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*c*c))/thetae);
    exp1 = exp(-sqrt(1+Pe(i)*Pe(i))/thetae);
    bes = besselk(2, 1/thetae);
    p = Pe(i);
    p3 = (p)^3;
    Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
    Fekappa(i) = Aekappa*(1 + ((Pe(i)*me*c/Pekappa)^2)/(kappa-3/2))^(-(kappa + 1))*4*pi*p3*Pe(i);
    
    
    %exp1 = exp(-sqrt(1+Pp(i)*Pp(i)/(mp*mp*c*c))/thetap);
    exp1 = exp(-sqrt(1+Pp(i)*Pp(i))/thetap);
    bes = besselk(2, 1/thetap);
    p = Pp(i);
    p3 = (p)^3;
    Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(i);
    Fpkappa(i) = Apkappa*(1 + (Pp(i)*me*c/Ppkappa)^2)^(-(kappa + 1));
end;

normp = (Fp(1)/(Pp(2)^2))*(Pp(2) - Pp(1));
norme = (Fe(1)/(Pe(2)^2))*(Pe(2) - Pe(1));
normkappae = (Fekappa(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normp = normp + (Fp(i)/(Pp(i)^2))*(Pp(i) - Pp(i-1));
    norme = norme + (Fe(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
    normkappae = normkappae + (Fekappa(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;

for i = 1:Np,
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
    Fekappa(i) = Fekappa(i)/normkappae;
    Fpold(i) = Fpold(i)*norm/normp;
    Feold(i) = Feold(i)*norm/norme;
end;

index1 = 165;
index2 = 175;

s = log(Fe(index1)/Fe(index2))/log(Pe(index2)/Pe(index1));
for i = 170:Np,
    Fekappa(i) = Fekappa(170)*(Pe(i)/Pe(170))^(-s);
end;

figure(1);
loglog(Pp(1:Np),Fp(1:Np),'red');
%loglog(Pp(1:Np),Fp(1:Np),'red',Pp(1:Np), Fpjuttner(1:Np), 'blue');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
%legend('Fp','Maxwell-Juttner fit','Location','southeast');
grid ;

figure(2);
%hold on;
loglog(Pe(1:Np),Fe(1:Np), 'red', Pe(1:Np),Fekappa(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np), 'red',Pe(1:Np), Fejuttner(1:Np), 'blue');
%loglog(Pe(1:Np)*me/(gam*beta*mp),Fe(1:Np), 'red',Pe(1:Np)*me/(gam*beta*mp), Feold(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np)*(me*c), 'red');
%loglog(Pe(1:Np),Fekappa(1:Np)*(me*c), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
%legend('Fe','Maxwell-Juttner fit','Location','southeast');
%legend('new','old','Location','southeast');
grid ;

dlmwrite('Pp3.dat',Pp,'delimiter','\n');
dlmwrite('Pe3.dat',Pe,'delimiter','\n');
dlmwrite('Fp3.dat',Fp,'delimiter','\n');
dlmwrite('Fe3.dat',Fe,'delimiter','\n');