clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.000';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');
B0 = 0.03030750;
Nx = size(Bx, 1);
Ny = size(By, 2);

Bnorm(1:Nx/2, 1:Ny) = 0;
Ediff(1:Nx/2, 1:Ny) = 0;
%Bperp(1:Nx, 1:Ny) = 0;
g = 1.5;
beta = sqrt(1-1/(g*g));

for i=1:Nx/2,
    for j = 1:Ny,
        Bnorm(i,j) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        Ediff(i,j) = sqrt((Ex(i,j))*(Ex(i,j)) + (Ey(i,j) +beta*Bz(i,j))*(Ey(i,j) +beta*Bz(i,j)) + (Ez(i,j) - beta*By(i,j))*(Ez(i,j) - beta*By(i,j)))/Bnorm(i,j);
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

%for i=1:Nx/2,
%    for j = 1:Nx/2,
%        k = rem(j,Ny) + 1;
%        Bnorm(i,j) = Bnorm(i,k);
%    end;
%end;

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 20;
mp = 1.67262*10^-24;
me = mp/mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;

omega = sqrt(4*pi*n*q*q/me);

rho = c0/(omega*Nskinlength);
rho = 0.2;
c1=0.45;

tau = c1*rho/c0;
samplingFactor = 5;
fieldFactor = me*rho/(q*tau*tau);
rho = rho*samplingFactor;

figure(1);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Bx*fieldFactor);
shading interp;
title ('Bx');
xlabel ('y');
ylabel ('x');
zlabel ('Bx');
grid ;

figure(2);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, By*fieldFactor);
shading interp;
title ('By');
xlabel ('y');
ylabel ('x');
zlabel ('By');
grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Bz*fieldFactor);
shading interp;
title ('Bz');
xlabel ('y');
ylabel ('x');
zlabel ('Bz');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ex/B0);
shading interp;
title ('Ex/B0');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
%grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ey/B0);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ez/B0);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;

figure(7);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx/2)*rho);
surf(X, Y, Bnorm/B0);
shading interp;
title ('B/B_0');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('B/B_0');
grid ;

%dlmwrite('B.dat',Bnorm);

%figure(8);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Bperp*fieldFactor);
%shading interp;
%title ('Bperp');
%xlabel ('y');
%ylabel ('x');
%zlabel ('B');
%grid ;

figure(9);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx/2)*rho);
surf(X, Y, Ediff);
shading interp;
title ('(E-vxB)/B');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('E');
grid ;