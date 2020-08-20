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
endx = 8000;

Fp(1:Np)=0;
Fe(1:Np)=0;


Fpold(1:Np)=0;
Feold(1:Np)=0;

Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;


mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;
gam = 1.5;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 2.6*10^11;
Temin = 10^9;
Temax = 10^12;
Tp = 2.4*10^12;
Tpmin = 0.5*10^12;
Tpmax = 10^14;
Pekappa = 14*me*c;
Ppkappa = mp*c;
kappa = 4;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 1.0;
fractionp = 1.0;

normp = 0;
norme = 0;
norm = 1;

for i = 2:Np,
    for j = startx:endx,
        Fp(i) = Fp(i) + fp(j,i)*g(i)*(g(i) + 1)^2/(g(i) - g(i-1));
        Fe(i) = Fe(i) + fe(j,i)*g(i)*(g(i) + 1)^2/(g(i) - g(i-1));
        Fpold(i) = Fpold(i) + fp(j,i);
        Feold(i) = Feold(i) + fe(j,i);
    end;
end;   

normp = Fp(1)*(g(2) - g(1));
norme = Fe(1)*(g(2) - g(1));

for i = 2:Np,
    normp = normp + Fp(i)*(g(i) - g(i-1))/(g(i) + 1)^2;
    norme = norme + Fe(i)*(g(i) - g(i-1))/(g(i) + 1)^2;
end;

for i = 1:Np,
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
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
        beta = sqrt(1 - 1/(1+g(i))^2);
        exp1 = exp(-(1+g(i))/thetap);
        Fpjuttner(i) = fractionp*((1+g(i))^4)*beta*(1.0/(thetap*bes))*exp1;
        s1 = s1 + ((Fpjuttner(i) - Fp(i))^2)*(g(i) - g(i-1));
    end;
    thetap = kB*Tp2/(mp*c*c);
    bes = besselk(2, 1/thetap);
     for i = index1:index2,
        beta = sqrt(1 - 1/(1+g(i))^2);
        exp1 = exp(-(1+g(i))/thetap);
        Fpjuttner(i) = fractionp*((1+g(i))^4)*beta*(1.0/(thetap*bes))*exp1;
        s2 = s2 + ((Fpjuttner(i) - Fp(i))^2)*(g(i) - g(i-1));
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
    beta = sqrt(1 - 1/(1+g(i))^2);
    exp1 = exp(-(1+g(i))/thetap);
    Fpjuttner(i) = fractionp*((1+g(i))^4)*beta*(1.0/(thetap*bes))*exp1;
end;

figure(1);
%loglog(Pp(1:Np),Fp(1:Np),'--','color','red');
plot(g(1:Np)+1,Fp(1:Np),'red',g(1:Np)+1, Fpjuttner(1:Np), 'blue');
title ('F_p');
xlabel ('\gamma');
ylabel ('Fp');
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
        beta = sqrt(1 - 1/(1+g(i))^2);
        exp1 = exp(-(1+g(i))/thetae);
        Fejuttner(i) = fractione*((1+g(i))^4)*beta*(1.0/(thetae*bes))*exp1;
        s1 = s1 + ((Fejuttner(i) - Fe(i))^2)*(g(i) - g(i-1));
    end;
    thetae = kB*Te2/(me*c*c);
    bes = besselk(2, 1/thetae);
    for i = index1:index2,
        beta = sqrt(1 - 1/(1+g(i))^2);
        exp1 = exp(-(1+g(i))/thetae);
        Fejuttner(i) = fractione*((1+g(i))^4)*beta*(1.0/(thetae*bes))*exp1;
        s2 = s2 + ((Fejuttner(i) - Fe(i))^2)*(g(i) - g(i-1));
    end;
    if(s1 < s2)
        Teright = Te2;
    else 
        Teleft = Te1;
    end;
end;
Te = (Teleft + Teright)/2;
%Te = 1E10;
thetae = kB*Te/(me*c*c);
bes = besselk(2, 1/thetae);
for i = 2:Np,   
    beta = sqrt(1 - 1/(1+g(i))^2);
    exp1 = exp(-(1+g(i))/thetae);
    Fejuttner(i) = fractione*((1+g(i))^4)*beta*(1.0/(thetae*bes))*exp1;
end;

figure(2);
%hold on;
%plot (Pe(1:Np),Fe(1:Np), 'red');
plot(g(1:Np)+1,Fe(1:Np), 'red',g(1:Np)+1, Fejuttner(1:Np), 'blue');
%loglog(Pe(1:Np)*me/(gam*beta*mp),Fe(1:Np), 'red',Pe(1:Np)*me/(gam*beta*mp), Feold(1:Np), 'blue');
%loglog(Pe(1:Np),Fe(1:Np)*(me*c), 'red');
%loglog(Pe(1:Np),Fekappa(1:Np)*(me*c), 'blue');
title ('F_e');
xlabel ('\gamma');
ylabel ('F_e');
legend('Fe','Maxwell-Juttner fit','Location','southeast');
%legend('new','old','Location','southeast');
grid ;