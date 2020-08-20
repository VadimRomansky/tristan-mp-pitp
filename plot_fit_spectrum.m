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
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
    Fekappa(i) = Fe(i)/normkappae;
    Fpold(i) = Fpold(i)*norm/normp;
    Feold(i) = Feold(i)*norm/norme;
end;

index1 = 120;
index2 = 180;

Tpleft = Tpmin;
Tpright = Tpmax;

for j = 1:20,
    Tp1 = Tpleft + (Tpright - Tpleft)/3;
    Tp2 = Tpleft + (Tpright - Tpleft)*2/3;
    s1 = 0;
    s2 = 0;
    thetap = kB*Tp1/(mp*c*c);
    bes = besselk(2, 1/thetap);
    for i = index1:index2,
        exp1 = exp(-sqrt(1+Pp(i)*Pp(i))/thetap);
        p = Pp(i);
        p3 = (p)^3;
        Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(i);
        s1 = s1 + ((Fpjuttner(i) - Fp(i))^2)*(Pp(i) - Pp(i-1))/(Pp(i)^4);
    end;
    thetap = kB*Tp2/(mp*c*c);
    bes = besselk(2, 1/thetap);
    for i = index1:index2,
        exp1 = exp(-sqrt(1+Pp(i)*Pp(i))/thetap);
        p = Pp(i);
        p3 = (p)^3;
        Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(i);
        s2 = s2 + ((Fpjuttner(i) - Fp(i))^2)*(Pp(i) - Pp(i-1))/(Pp(i)^4);
    end;
    if(s1 < s2)
        Tpright = Tp2;
    else 
        Tpleft = Tp1;
    end;
end;
Tp = (Tpleft + Tpright)/2;
thetap = kB*Tp/(mp*c*c);
bes = besselk(2, 1/thetap);
for i = 2:Np,   
    exp1 = exp(-sqrt(1+Pp(i)*Pp(i))/thetap);
    p = Pp(i);
    p3 = (p)^3;
    Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(i);
end;

figure(1);
%loglog(Pp(1:Np),Fp(1:Np),'--','color','red');
loglog(Pp(1:Np),Fp(1:Np),'red',Pp(1:Np), Fpjuttner(1:Np), 'blue');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
legend('Fp','Maxwell-Juttner fit','Location','southeast');
grid ;





Teleft = Temin;
Teright = Temax;

for j = 1:20,
    Te1 = Teleft + (Teright - Teleft)/3;
    Te2 = Teleft + (Teright - Teleft)*2/3;
    s1 = 0;
    s2 = 0;
    thetae = kB*Te1/(me*c*c);
    bes = besselk(2, 1/thetae);
    for i = index1:index2,
        exp1 = exp(-sqrt(1+Pe(i)*Pe(i))/thetae);
        p = Pe(i);
        p3 = (p)^3;
        Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
        s1 = s1 + ((Fejuttner(i) - Fe(i))^2)*(Pe(i) - Pe(i-1))/(Pe(i)^4);
    end;
    thetae = kB*Te2/(me*c*c);
    bes = besselk(2, 1/thetae);
    for i = index1:index2,
        exp1 = exp(-sqrt(1+Pe(i)*Pe(i))/thetae);
        p = Pe(i);
        p3 = (p)^3;
        Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
        s2 = s2 + ((Fejuttner(i) - Fe(i))^2)*(Pe(i) - Pe(i-1))/(Pe(i)^4);
    end;
    if(s1 < s2)
        Teright = Te2;
    else 
        Teleft = Te1;
    end;
end;
Te = (Teleft + Teright)/2;
thetae = kB*Te/(me*c*c);
bes = besselk(2, 1/thetae);
for i = 2:Np,   
    exp1 = exp(-sqrt(1+Pe(i)*Pe(i))/thetae);
    p = Pe(i);
    p3 = (p)^3;
    Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
end;

figure(2);
%hold on;
%plot (Pe(1:Np),Fe(1:Np), 'red');
loglog(Pe(1:Np),Fe(1:Np)*10^5, 'red',Pe(1:Np), Fejuttner(1:Np)*10^5, 'blue', Pe(1:Np),Fekappa(1:Np), 'green');
%loglog(Pe(1:Np)*me/(gam*beta*mp),Fe(1:Np), 'red',Pe(1:Np)*me/(gam*beta*mp), Feold(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np)*(me*c), 'red');
%loglog(Pe(1:Np),Fekappa(1:Np)*(me*c), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
legend('Fe','Maxwell-Juttner fit','kappa','Location','southeast');
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