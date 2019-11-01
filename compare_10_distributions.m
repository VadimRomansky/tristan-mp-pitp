clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.010';
Nd = 10;
start = 0;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'t*{\Omega} = 30','t*{\Omega} = 60','t*{\Omega} = 90', 't*{\Omega} = 120', 't*{\Omega} = 150','t*{\Omega} = 180'};
LegendTitle = {'{\theta} = 20', '{\theta} = 30','{\theta} = 40', '{\theta} = 42', '{\theta} = 44', '{\theta} = 46','{\theta} = 48', '{\theta} = 50', '{\theta} = 55', '{\theta} = 60'};


full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = size(fp,1);
startx = 1;
endx = Nx/4;

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
        for k = startx:endx,
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
set(0, 'DefaultLineLineWidth', 1.5);
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
for j=1:Nd,
    plot (Pp(j, 1:Np),Fp(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest');
grid ;

spectrum(1:Np,1:20) = 0;
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
    
    spectrum(i,11) = Pe(6,i);
    spectrum(i,12) = Fe(6,i);
    
    spectrum(i,13) = Pe(7,i);
    spectrum(i,14) = Fe(7,i);
    
    spectrum(i,15) = Pe(8,i);
    spectrum(i,16) = Fe(8,i);
    
    spectrum(i,17) = Pe(9,i);
    spectrum(i,18) = Fe(9,i);
    
    spectrum(i,19) = Pe(10,i);
    spectrum(i,20) = Fe(10,i);
end;
dlmwrite('spectrum10.dat',spectrum,'delimiter',' ');