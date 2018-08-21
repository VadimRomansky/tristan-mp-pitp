clear;
directory_name = './output1/';
file_name = 'spect';
file_number = '.014';
Nd = 10;
start = 0;

full_name = strcat(directory_name, file_name, num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = size(fp,1);

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
        for k = 1:Nx,
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

Color = {'red','blue','green','black','yellow',[.5,.5,.5],'cyan',[.3,.7,.3], 'magenta',[1.0,.5,0],[.75,0.0,.7],[.7,.3,.3]};
LegendTitle = {'filter 16','filter 32','filter 16','filter 32','filter 16','filter 32','filter 16','filter 32','filter 16','filter 32'};

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
for j=1:Nd,
    plot (Pp(j, 1:Np),Fp(j, 1:Np),'color',Color{j});
end;
legend('{\theta} = 0^{\circ}','{\theta} = 10^{\circ}','{\theta} = 20^{\circ}','{\theta} = 30^{\circ}','{\theta} = 40^{\circ}','{\theta} = 50^{\circ}','{\theta} = 60^{\circ}','{\theta} = 70^{\circ}', '{\theta} = 80^{\circ}','{\theta} = 90^{\circ}','Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
for j=1:Nd,
    plot (Pe(j, 1:Np),Fe(j, 1:Np),'color',Color{j});
end;
legend('{\theta} = 0^{\circ}','{\theta} = 10^{\circ}','{\theta} = 20^{\circ}','{\theta} = 30^{\circ}','{\theta} = 40^{\circ}','{\theta} = 50^{\circ}','{\theta} = 60^{\circ}','{\theta} = 70^{\circ}','{\theta} = 80^{\circ}','{\theta} = 90^{\circ}','Location','southeast');
grid ;