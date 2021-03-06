clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.020';
Nd = 5;
start = 0;

Color = {'green','red','blue','black','magenta'};
%LegendTitle = {'l = 6 rg','l = 11 rg','l = 22 rg', 'l = 33 rg', 'l = 45 rg'};
LegendTitle = {'1','2','3','4','5'};

full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specprest');
Np = size(fp,2);
Nx = size(fp,1);

samplingfactor = 20;
startx(1:Nd) = 0;
endx(1:Nd) = 0;

startx(1) = 500*samplingfactor;
endx(1) = startx(1) + 320*samplingfactor;

startx(2) = 700*samplingfactor;
endx(2) = startx(2) + 320*samplingfactor;

startx(3) = 900*samplingfactor;
endx(3) = startx(3) + 320*samplingfactor;

startx(4) = 1100*samplingfactor;
endx(4) = startx(4) + 320*samplingfactor;

startx(5) = 1300*samplingfactor;
endx(5) = startx(5) + 320*samplingfactor;

g(1:Nd,1:Np) = 0;
Fp(1:Nd,1:Np)=0;
Fe(1:Nd,1:Np)=0;
Pp(1:Nd,1:Np)=0;
Pe(1:Nd,1:Np)=0;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    fp = hdf5read(full_name,'specprest');
    fe = hdf5read(full_name,'specerest');
    gam=hdf5read(full_name,'gamma');
    for i = 2:Np,
        g(j, i) = gam(i);
        Pp(j,i) = sqrt((g(j,i)+1)^2 - 1);
        Pe(j,i) = sqrt((g(j,i)+1)^2 - 1);
        for k = startx(j):endx(j),
            Fp(j,i) = Fp(j,i) + fp(k,i)*gam(i)/(gam(i) - gam(i-1));
            Fe(j,i) = Fe(j,i) + fe(k,i)*gam(i)/(gam(i) - gam(i-1));
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
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5},'Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e p^4');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5},'Location','southeast');
grid ;

spectrum(1:Np,1:10) = 0;
for i = 1:Np,
    spectrum(i,1) = Pe(1,i);
    spectrum(i,2) = Fe(1,i);
    
    spectrum(i,3) = Pe(2,i);
    spectrum(i,4) = Fe(2,i);
    
    spectrum(i,5) = Pe(3,i);
    spectrum(i,6) = Fe(3,i);
    
    spectrum(i,7) = Pe(4,i);
    spectrum(i,8) = Fe(4,i);
    
    spectrum(i,9) = Pe(5,i);
    spectrum(i,10) = Fe(5,i);
end;
dlmwrite('spectrum5.dat',spectrum,'delimiter',' ');