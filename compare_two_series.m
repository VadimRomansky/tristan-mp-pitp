clear;
directory_name = './output5/';
file_name = 'spect';
file_number = '.000';
Nd = 8;
Ns = 2;
start = 1;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'1','2','3','4','5','6','7','8','9','10'};
FileNumbers = {'.000','.005','.010','.015','.020','.025','.030','.035','.040','.045'};

full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = fix(size(fp,1)/4);
%Nx = 12500;

g(1:Ns,1:Nd,1:Np) = 0;
Fp(1:Ns,1:Nd,1:Np)=0;
Fe(1:Ns,1:Nd,1:Np)=0;
Pp(1:Ns,1:Nd,1:Np)=0;
Pe(1:Ns,1:Nd,1:Np)=0;

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
fractione = 1.0;
fractionp = 1.0;

for m = 1:Ns,
    for j = 1:Nd,
        full_name = strcat(directory_name, file_name, num2str(m), FileNumbers{j});
        fp = hdf5read(full_name,'specp');
        fe = hdf5read(full_name,'spece');
        gam=hdf5read(full_name,'gamma');
        for i = 1:Np,
            g(m,j, i) = gam(i);
            Pp(m,j,i) = sqrt((g(m,j,i)+1)^2 - 1);
            Pe(m,j,i) = sqrt((g(m,j,i)+1)^2 - 1);
            for k = 1:Nx,
                Fp(m,j,i) = Fp(m,j,i) + fp(k,i);
                Fe(m,j,i) = Fe(m,j,i) + fe(k,i);
            end;
            Fp(m,j,i)=Fp(m,j,i)*(Pp(m,j,i)^3)/(1+g(m,j,i));
            Fe(m,j,i)=Fe(m,j,i)*(Pe(m,j,i)^3)/(1+g(m,j,i));
        end;
    end;
end;

for i = 1:Np,
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

for m = 1:Ns,
for j = 1:Nd,
    normp = (Fp(m,j,1)/(Pp(m,j,2)^2))*(Pp(m,j,2) - Pp(m,j,1));
    norme = (Fe(m,j,1)/(Pe(m,j,2)^2))*(Pe(m,j,2) - Pe(m,j,1));

    for i = 2:Np,
        normp = normp + (Fp(m,j,i)/(Pp(m,j,i)^2))*(Pp(m,j,i) - Pp(m,j,i-1));
        norme = norme + (Fe(m,j,i)/(Pe(m,j,i)^2))*(Pe(m,j,i) - Pe(m,j,i-1));
    end;

    for i = 1:Np,
        Fp(m,j,i) = Fp(m,j,i)*norm/normp;
        Fe(m,j,i) = Fe(m,j,i)*norm/norme;
    end;
end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
for m = 1:Ns,
for j=1:Nd,
    p(1:Np) = 0;
    f(1:Np) = 0;
    for k = 1:Np,
        p(k) = Pp(m,j, k);
        f(k) = Fp(m,j, k);
    end;
    plot (p(1:Np),f(1:Np),'color',Color{j});
end;
end;
%plot (Pp(j, 1:Np),Fpjuttner(1:Np),'color','green');
legend('10rg 1','10rg 2','10rg 3','10rg 4','10rg 5','10rg 6','10rg 7','10rg 8','5rg 1','5rg 2','5rg 3','5rg 4','5rg 5','5rg 6','5rg 7','5rg 8','Location','southeast');
%legend(LegendTitle{1}, LegendTitle{2},'juttner','Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
for m = 1:Ns,
for j=1:Nd,
    p(1:Np) = 0;
    f(1:Np) = 0;
    for k = 1:Np,
        p(k) = Pe(m,j, k);
        f(k) = Fe(m,j, k);
    end;
    plot (p(1:Np),f(1:Np),'color',Color{j});
end;
end;
%plot (Pp(j, 1:Np),Fpjuttner(1:Np),'color','green');
legend('10rg 1','10rg 2','10rg 3','10rg 4','10rg 5','10rg 6','10rg 7','10rg 8','5rg 1','5rg 2','5rg 3','5rg 4','5rg 5','5rg 6','5rg 7','5rg 8', LegendTitle{2},'Location','southeast');
%legend(LegendTitle{1}, LegendTitle{2},'juttner','Location','southeast');
grid ;

spectrum(1:Np,1:4) = 0;
for i = 1:Np,
    spectrum(i,1) = Pe(1,1,i);
    spectrum(i,2) = Fe(1,1,i);
    
    spectrum(i,3) = Pe(1,2,i);
    spectrum(i,4) = Fe(1,2,i);
end;
dlmwrite('spectrum2.dat',spectrum,'delimiter',' ');