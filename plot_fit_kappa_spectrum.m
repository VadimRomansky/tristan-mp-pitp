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

startx = 1000;
endx = 5000;

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
Temin = 10^11;
Temax = 10^12;
Tp = 2.4*10^12;
Tpmin = 0.5*10^12;
Tpmax = 10^13;
Pekappa = 14*me*c;
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
    p3 = Pe(i)*Pe(i)*Pe(i);
    %Fekappa(i) = Aekappa*(1 + ((Pe(i)*me*c/Pekappa)^2)/(kappa-3/2))^(-(kappa + 1))*4*pi*p3*Pe(i);
    Fekappa(i) = Aekappa*(Pe(i)^4)*(1 + (Pe(i)*me*c/Pekappa)^2)^(-(kappa + 1));
    Fpkappa(i) = Apkappa*(1 + (Pp(i)*mp*c/Ppkappa)^2)^(-(kappa + 1));
    Fpold(i)=Fpold(i)*(1/(gam*beta))*(Pp(i)^3)/(1+g(i));
    Feold(i)=Feold(i)*(me/(gam*beta*mp))*(Pe(i)^3)/(1+g(i));
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
    Fp(i) = Fp(i)/normp;
    Fe(i) = Fe(i)/norme;
    Fekappa(i) = Fekappa(i)/normkappae;
    Fpold(i) = Fpold(i)/normp;
    Feold(i) = Feold(i)/norme;
end;

type kappafit;
fun = @(x)kappafit(x,Pe,Fe, Np);

x0(1:2) = 0;
x0(1) = Pekappa/(me*c);
x0(2) = kappa;

bestx(1) = Pekappa/(me*c);
bestx(2) = kappa;

bestx = fminsearch(fun,x0);

Fekappa(1:Np) = 0;
for i = 1:Np,
    Fekappa(i) = (Pe(i)^4)*(1 + (Pe(i)/bestx(1))^2)^(-(bestx(2) + 1));
end;

normkappae = (Fekappa(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normkappae = normkappae + (Fekappa(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;
for i = 1:Np,
    Fekappa(i) = Fekappa(i)/normkappae;
end;



figure(2);
%hold on;
plot (Pe(1:Np),Fe(1:Np), 'red', Pe(1:Np),Fekappa(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np)*10^5, 'red',Pe(1:Np), Fejuttner(1:Np)*10^5, 'blue', Pe(1:Np),Fekappa(1:Np), 'green');
%loglog(Pe(1:Np)*me/(gam*beta*mp),Fe(1:Np), 'red',Pe(1:Np)*me/(gam*beta*mp), Feold(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np)*(me*c), 'red');
%loglog(Pe(1:Np),Fekappa(1:Np)*(me*c), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
legend('Fe','kappa fit','kappa','Location','southeast');
%legend('new','old','Location','southeast');
grid ;

tempSpectrum(1:Np,2) = 0;
for i = 1:Np,
    tempSpectrum(i,1) = Pe(i);
    tempSpectrum(i,2) = Fejuttner(i);
end;

dlmwrite('Femax.dat',tempSpectrum,'delimiter',' ');

dlmwrite('Pp.dat',Pp,'delimiter','\n');
dlmwrite('Pe.dat',Pe,'delimiter','\n');
dlmwrite('Fp.dat',Fp,'delimiter','\n');
dlmwrite('Fe.dat',Fejuttner,'delimiter','\n');