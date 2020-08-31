clear;
directory_name = './output/';
file_name = 'flds';
file_number = '.tot.010';
Nd = 2;
start = 0;
Color = {'red','blue'};
LegendTitle = {'out','in'};

full_name = strcat(directory_name, file_name, num2str(start), file_number);
Bx0 = hdf5read(full_name,'bx');
By0 = hdf5read(full_name,'by');
Bz0 = hdf5read(full_name,'bz');
Ex0 = hdf5read(full_name,'ex');
Ey0 = hdf5read(full_name,'ey');
Ez0 = hdf5read(full_name,'ez');
Nx = size(Bx0, 1);
Ny = size(Bx0, 2);

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 20;
mp = 1.67262*10^-24;
me = mp/mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;

omega = sqrt(4*pi*n*q*q/me);

rho = c0/(omega*Nskinlength);
c1=0.45;

tau = c1*rho/c0;

fieldFactor = me*rho/(q*tau*tau);

samplingFactor = 20;

rho = rho*samplingFactor;


Bx(1:Nd,1:Nx)=0;
By(1:Nd,1:Nx)=0;
Bz(1:Nd,1:Nx)=0;
Ex(1:Nd,1:Nx)=0;
Ey(1:Nd,1:Nx)=0;
Ez(1:Nd,1:Nx)=0;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(start + j-1), file_number);
    Bx0 = hdf5read(full_name,'bx');
    By0 = hdf5read(full_name,'by');
    Bz0 = hdf5read(full_name,'bz');
    Ex0 = hdf5read(full_name,'ex');
    Ey0 = hdf5read(full_name,'ey');
    Ez0 = hdf5read(full_name,'ez');
    for i = 1:Nx,
        for k = 1:Ny,
            Bx(j,i)=Bx(j,1) + Bx0(i, k);
            By(j,i)=By(j,i) + By0(i, k);
            Bz(j,i)=Bz(j,i) + Bz0(i, k);
            Ex(j,i)=Ex(j,i) + Ex0(i, k);
            Ey(j,i)=Ey(j,i) + Ey0(i, k);
            Ez(j,i)=Ez(j,i) + Ez0(i, k);
        end;
        Bx(j,i)=Bx(j,1)/Ny;
        By(j,i)=By(j,i)/Ny;
        Bz(j,i)=Bz(j,i)/Ny;
        Ex(j,i)=Ex(j,i)/Ny;
        Ey(j,i)=Ey(j,i)/Ny;
        Ez(j,i)=Ez(j,i)/Ny;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, Bx(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle(1), LegendTitle(2),'Location','southeast');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, By(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle(1), LegendTitle(2),'Location','southeast');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, Bz(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle(1), LegendTitle(2),'Location','southeast');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, Ex(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle(1), LegendTitle(2),'Location','southeast');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, Ey(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2},'Location','southeast');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
hold on;
for j=1:Nd,
    plot ((1:Nx)*samplingFactor, Ez(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend(LegendTitle{1}, LegendTitle{2},'Location','southeast');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;

