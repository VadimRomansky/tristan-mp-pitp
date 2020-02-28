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

startx = 10;
endx = 40000;

Fp(1:Np)=0;
Fe(1:Np)=0;

Fpold(1:Np)=0;
Feold(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;
Fekappa(1:Np) = 0;
Fpkappa(1:Np) = 0;

me = 0.91*10^-27;
mass_ratio = 100;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 5*10^9;
Tp = 3.5*10^10;
Pekappa = 14*me*c;
Ppkappa = mp*c;
kappa = 4;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

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
    Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
    Fpold(i)=Fpold(i)*(Pp(i)^3)/(1+g(i));
    Feold(i)=Feold(i)*(Pe(i)^3)/(1+g(i));
    
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
    Fpold(i) = Fpold(i)*norm/normp;
    Feold(i) = Feold(i)*norm/norme;
end;

figure(1);
loglog(Pp(1:Np),Fp(1:Np),'--','color','red');
%plot (Pp(1:Np),Fp(1:Np), 'red',Pp(1:Np), Fpjuttner(1:Np), 'blue');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
legend('Fp','Location','southeast');
grid ;

figure(2);
%plot (Pe(1:Np),Fe(1:Np), 'red');
%loglog(Pe(1:Np),Fe(1:Np), 'red',Pe(1:Np), Fejuttner(1:Np), 'blue', Pe(1:Np), Fekappa(1:Np),'green');
loglog(Pe(1:Np),Fe(1:Np), 'red',Pe(1:Np), Feold(1:Np), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
legend('new','old','Location','southeast');
grid ;

dlmwrite('Pp.dat',Pp,'delimiter','\n');
dlmwrite('Pe.dat',Pe,'delimiter','\n');
dlmwrite('Fp.dat',Fp,'delimiter','\n');
dlmwrite('Fe.dat',Fe,'delimiter','\n');