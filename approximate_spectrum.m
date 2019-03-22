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

startx = 1;
endx = fix(Nx/4);

startPowerP = 130;
endPowerP = 155;

startPowerE = 180;
endPowerE = 190;

Fp(1:Np)=0;
Fe(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

Fpa(1:Np)=0;
Fea(1:Np)=0;
gammap = 1;
gammae = 1;
gp(1:Np) = 0;
ge(1:Np) = 0;
ap = 1;
ae = 1;

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
    Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
    
    %exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*c*c))/thetae);
    exp1 = exp(-sqrt(1+Pe(i)*Pe(i))/thetae);
    bes = besselk(2, 1/thetae);
    p = Pe(i);
    p3 = (p)^3;
    Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(i);
    
    
    
    %exp1 = exp(-sqrt(1+Pp(i)*Pp(i)/(mp*mp*c*c))/thetap);
    exp1 = exp(-sqrt(1+Pp(i)*Pp(i))/thetap);
    bes = besselk(2, 1/thetap);
    p = Pp(i);
    p3 = (p)^3;
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

Fpa(startPowerP) = Fp(startPowerP);
Fpa(endPowerP) = Fp(endPowerP);

Fea(startPowerE) = Fe(startPowerE);
Fea(endPowerE) = Fe(endPowerE);

gammap = log(Fpa(startPowerP)/Fpa(endPowerP))/log(Pp(startPowerP)/Pp(endPowerP));
gammae = log(Fea(startPowerE)/Fea(endPowerE))/log(Pe(startPowerE)/Pe(endPowerE));

ap = exp(log(Fpa(startPowerP)) - gammap*log(Pp(startPowerP)));
ae = exp(log(Fea(startPowerE)) - gammae*log(Pe(startPowerE)));

for i = startPowerP:endPowerP,
    Fpa(i) = ap*(Pp(i)^gammap);
end;

for i = startPowerE:endPowerE,
    Fea(i) = ae*(Pe(i)^gammae);
end;

for i = 2:Np,
    gp(i) = log(Fp(i)/Fp(i-1))/log(Pp(i)/Pp(i-1)) - 4;
    ge(i) = log(Fe(i)/Fe(i-1))/log(Pe(i)/Pe(i-1)) - 4;
end;

figure(1);
%plot (Pp(1:Np),Fp(1:Np), 'red');
plot (Pp(1:Np),Fp(1:Np), 'red',Pp(startPowerP:endPowerP), Fpa(startPowerP:endPowerP), 'blue');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
name = strcat('approximation gamma = ',num2str(gammap - 4));
legend('Fp', name,'Location','southeast');
grid ;

figure(2);
plot (Pe(1:Np),Fe(1:Np), 'red',Pe(startPowerE:endPowerE), Fea(startPowerE:endPowerE),'blue');
%plot (Pe(1:Np),Fe(1:Np), 'red',Pe(1:Np), Fejuttner(1:Np), 'blue');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
name = strcat('approximation gamma = ',num2str(gammae - 4));
legend('Fe', name,'Location','southeast');
grid ;

figure(3);
plot (Pp(1:Np),gp(1:Np), 'red');
title ('{\gamma}_p');
xlabel ('p/{m_p c}');
ylabel ('{\gamma}');
grid ;

figure(4);
plot (Pe(1:Np),ge(1:Np), 'red');
title ('{\gamma}_e');
xlabel ('p/{m_e c}');
ylabel ('{\gamma}');
grid ;

outputArrayE(1:Np, 1:2) = 0;
outputArrayP(1:Np, 1:2) = 0;

for i = 1:Np,
    outputArrayE(i, 1) = Pe(i);
    outputArrayE(i, 2) = Fe(i)/(Pe(i)^2);
    outputArrayP(i, 1) = Pp(i);
    outputArrayP(i, 2) = Fp(i)/(Pp(i)^2);
end;

dlmwrite('Fe_perpendisular.dat',outputArrayE, 'delimiter', ' ');
dlmwrite('Fp_perpendicular.dat',outputArrayP, 'delimiter', ' ');