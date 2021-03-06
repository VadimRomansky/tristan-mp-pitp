clear;
directory_name = './output3/';
file_name = 'spect';
file_number = '.010';
Nd = 4;
start = 0;

Color = {'red','blue','green','black','magenta'};
LegendTitle = {'{\theta} = 0', '{\theta} = 10','{\theta} = 20', '{\theta} = 30'};

full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = size(fp,1);
startx = 10000;
endx = 15000;

g(1:Nd,1:Np) = 0;
Fp(1:Nd,1:Np)=0;
Fe(1:Nd,1:Np)=0;
Pp(1:Nd,1:Np)=0;
Pe(1:Nd,1:Np)=0;

Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

me = 0.91*10^-27;
mass_ratio = 100;
mp = me*mass_ratio;

Te = 1.2*10^10;
Tp = 1*10^11;
kB = 1.3806488*10^-16;
c = 2.99792458*10^10;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 1;
fractionp = 1;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    fp = hdf5read(full_name,'specp');
    fe = hdf5read(full_name,'spece');
    gam=hdf5read(full_name,'gamma');
    for i = 2:Np,
        g(j, i) = gam(i);
        Pp(j,i) = sqrt((g(j,i)+1)^2 - 1);
        Pe(j,i) = sqrt((g(j,i)+1)^2 - 1);
        for k = startx:endx,
            Fp(j,i) = Fp(j,i) + fp(k,i)*gam(i)/(gam(i) - gam(i-1));
            Fe(j,i) = Fe(j,i) + fe(k,i)*gam(i)/(gam(i) - gam(i-1));
        end;
        Fp(j,i)=Fp(j,i)*(Pp(j,i)^3)/(1+g(j,i));
        Fe(j,i)=Fe(j,i)*(Pe(j,i)^3)/(1+g(j,i));
    end;
end;

for i = 1:Np,
    %exp1 = exp(-sqrt(1+Pe(i)*Pe(i)/(me*me*c*c))/thetae);
    exp1 = exp(-sqrt(1+Pe(1,i)*Pe(1,i))/thetae);
    bes = besselk(2, 1/thetae);
    p = Pe(1,i);
    p3 = (p)^3;
    Fejuttner(i) = fractione*(1.0/(thetae*bes))*exp1*p3*Pe(1,i);
    
    
    
    %exp1 = exp(-sqrt(1+Pp(i)*Pp(i)/(mp*mp*c*c))/thetap);
    exp1 = exp(-sqrt(1+Pp(1,i)*Pp(1,i))/thetap);
    bes = besselk(2, 1/thetap);
    p = Pp(1,i);
    p3 = (p)^3;
    Fpjuttner(i) = fractionp*(1.0/(thetap*bes))*exp1*p3*Pp(1,i);
end;

norm = 1;

for j = 1:Nd,
    normp = (Fp(j,1)/(Pp(j,2)^2))*(Pp(j,2) - Pp(j,1));
    norme = (Fe(j,1)/(Pe(j,2)^2))*(Pe(j,2) - Pe(j,1));

    for i = 2:Np,
        normp = normp + (Fp(j,i)/(Pp(j,i)^2))*(Pp(j,i) - Pp(j,i-1));
        norme = norme + (Fe(j,i)/(Pe(j,i)^2))*(Pe(j,i) - Pe(j,i-1));
    end;

    for i = 1:Np,
        Fp(j,i) = Fp(j,i)*norm/normp;
        Fe(j,i) = Fe(j,i)*norm/norme;
    end;
end;


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1);
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
for j=1:Nd,
    plot (Pp(j, 1:Np),Fp(j, 1:Np),'color',Color{j});
end;
%plot (Pp(1, 1:Np),Fpjuttner(1:Np),'color',Color{Nd+1});
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},'Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
%plot (Pe(1, 1:Np),Fejuttner(1:Np),'color',Color{Nd+1});
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},'Location','southeast');
grid ;

spectrum(1:Np,1:8) = 0;
for i = 1:Np,
    spectrum(i,1) = Pe(1,i);
    spectrum(i,2) = Fe(1,i);
    
    spectrum(i,3) = Pe(2,i);
    spectrum(i,4) = Fe(2,i);
    
    spectrum(i,5) = Pe(3,i);
    spectrum(i,6) = Fe(3,i);
    
    spectrum(i,7) = Pe(4,i);
    spectrum(i,8) = Fe(4,i);
end;
dlmwrite('spectrum4.dat',spectrum,'delimiter',' ');