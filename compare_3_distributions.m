clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.007';
Nd = 3;
start = 0;

Color = {'red','blue','green'};
LegendTitle = {'mp/me 100','mp/me 50','mp/me 25'};


full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = size(fp,1);
startx(1:Nd) = 1;
endx(1:Nd) = 0;
endx(1) = 3000;
endx(2) = 3000;
endx(3) = 3000;
startx(1) = 1000;
startx(2) = 1000;
startx(3) = 1000;

g(1:Nd,1:Np) = 0;
Fp(1:Nd,1:Np)=0;
Fe(1:Nd,1:Np)=0;
Pp(1:Nd,1:Np)=0;
Pe(1:Nd,1:Np)=0;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    fp = hdf5read(full_name,'specp');
    fe = hdf5read(full_name,'spece');
    gam=hdf5read(full_name,'gamma');
    for i = 1:Np,
        g(j, i) = gam(i);
        Pp(j,i) = sqrt((g(j,i)+1)^2 - 1);
        Pe(j,i) = sqrt((g(j,i)+1)^2 - 1);
        for k = startx(j):endx(j),
            Fp(j,i) = Fp(j,i) + fp(k,i);
            Fe(j,i) = Fe(j,i) + fe(k,i);
        end;
        Fp(j,i)=Fp(j,i)*(Pp(j,i)^3)/(1+g(j,i));
        Fe(j,i)=Fe(j,i)*(Pe(j,i)^3)/(1+g(j,i));
    end;
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
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp p^4');
for j=1:Nd,
    plot (Pp(j, 1:Np),Fp(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3},'Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e p^4');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3},'Location','southeast');
grid ;

factor(1:Nd) = 0;
factor(1) = 1;
factor(2) = sqrt(2);
factor(3) = 2;
figure(3);
hold on;
title ('F_e');
xlabel ('p/{m_e c} * sqrt(me/me0)');
ylabel ('F_e p^4  * sqrt(me/me0)');
for j=1:Nd,
    plot (Pe(j, 1:Np)*factor(j),Fe(j, 1:Np)*factor(j),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3},'Location','southeast');
grid ;

spectrum(1:Np,1:6) = 0;
for i = 1:Np,
    spectrum(i,1) = Pe(1,i);
    spectrum(i,2) = Fe(1,i);
    
    spectrum(i,3) = Pe(2,i);
    spectrum(i,4) = Fe(2,i);
    
    spectrum(i,5) = Pe(3,i);
    spectrum(i,6) = Fe(3,i);
end;
dlmwrite('spectrum3.dat',spectrum,'delimiter',' ');