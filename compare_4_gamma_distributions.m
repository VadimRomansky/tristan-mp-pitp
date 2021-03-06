clear;
directory_name = './output2/';
file_name = 'spect';
file_number = '.050';
Nd = 4;
start = 0;


full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = size(fp,1);
startx = 4000;
endx = 6000;

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
        Pp(j,i) = gam(i);
        Pe(j,i) = gam(i);
        for k = startx:endx,
            Fp(j,i) = Fp(j,i) + fp(k,i);
            Fe(j,i) = Fe(j,i) + fe(k,i);
        end;
        Fp(j,i)=Fp(j,i)*g(j,i);
        Fe(j,i)=Fe(j,i)*g(j,i);
    end;
end;

norm = 1;

for j = 1:Nd,
    normp = (Fp(j,1)/g(j,2))*(g(j,2) - g(j,1));
    norme = (Fe(j,1)/g(j,2))*(g(j,2) - g(j,1));

    for i = 2:Np,
        normp = normp + (Fp(j,i)/g(j,i))*(g(j,i) - g(j,i-1));
        norme = norme + (Fe(j,i)/g(j,i))*(g(j,i) - g(j,i-1));
    end;

    for i = 1:Np,
        Fp(j,i) = Fp(j,i)*norm/normp;
        Fe(j,i) = Fe(j,i)*norm/norme;
    end;
end;

Color = {'red','blue','green', 'black'};
LegendTitle = {'filter 2','filter 16','filter 32','filter 128'};

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('F_p');
xlabel ('gamma - 1');
ylabel ('Fp*(gamma-1)');
for j=1:Nd,
    plot (Pp(j, 1:Np),Fp(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},'Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('gamma - 1');
ylabel ('F_e*(gamma -1)');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4},'Location','southeast');
grid ;