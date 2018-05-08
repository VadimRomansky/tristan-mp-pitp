clear;
directory_name = './output1/';
file_name = 'flds';
file_number = '.tot.002';
Nd = 2;

full_name = strcat(directory_name, file_name, num2str(0), file_number);
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

samplingFactor = 1;

rho = rho*samplingFactor;


Bx(1:Nd,1:Nx)=0;
By(1:Nd,1:Nx)=0;
Bz(1:Nd,1:Nx)=0;
Ex(1:Nd,1:Nx)=0;
Ey(1:Nd,1:Nx)=0;
Ez(1:Nd,1:Nx)=0;


for j = 1:Nd,
    full_name = strcat(directory_name, file_name, num2str(j-1), file_number);
    Bx0 = hdf5read(full_name,'bx');
    By0 = hdf5read(full_name,'by');
    Bz0 = hdf5read(full_name,'bz');
    Ex0 = hdf5read(full_name,'ex');
    Ey0 = hdf5read(full_name,'ey');
    Ez0 = hdf5read(full_name,'ez');
    for i = 1:Nx,
        Bx(j,i)=Bx0(i, fix(Ny/2));
        By(j,i)=By0(i, fix(Ny/2));
        Bz(j,i)=Bz0(i, fix(Ny/2));
        Ex(j,i)=Ex0(i, fix(Ny/2));
        Ey(j,i)=Ey0(i, fix(Ny/2));
        Ez(j,i)=Ez0(i, fix(Ny/2));
    end;
end;

Color = {'red','blue'};

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 

figure(1);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Bx(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, By(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Bz(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Ex(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Ey(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
hold on;
for j=1:Nd,
    plot ((1:Nx)*rho, Ez(j, 1:Nx)*fieldFactor,'color', Color{j});
end;
legend('turbulent','regular','Location','southeast');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;

