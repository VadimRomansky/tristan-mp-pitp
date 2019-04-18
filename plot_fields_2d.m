clear;
directory_name = './output4/';
file_name = 'flds4.tot';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
%Ex = hdf5read(full_name,'ex');
%Ey = hdf5read(full_name,'ey');
%Ez = hdf5read(full_name,'ez');
B0 = 0.03030750;
Nx = size(Bx, 1);
Ny = size(By, 2);

Bnorm(1:Nx/2, 1:Ny) = 0;
%Bperp(1:Nx, 1:Ny) = 0;

for i=1:Nx/2,
    for j = 1:Ny,
        Bnorm(i,j) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
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

%figure(1);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Bx*fieldFactor);
%shading interp;
%title ('Bx');
%xlabel ('y');
%ylabel ('x');
%zlabel ('Bx');
%grid ;

%figure(2);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, By*fieldFactor);
%shading interp;
%title ('By');
%xlabel ('y');
%ylabel ('x');
%zlabel ('By');
%grid ;

%figure(3);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Bz*fieldFactor);
%shading interp;
%title ('Bz');
%xlabel ('y');
%ylabel ('x');
%zlabel ('Bz');
%grid ;

%figure(4);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Ex*fieldFactor);
%shading interp;
%title ('Ex');
%xlabel ('y');
%ylabel ('x');
%zlabel ('Ex');
%grid ;

%figure(5);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Ey*fieldFactor);
%shading interp;
%title ('Ey');
%xlabel ('y');
%ylabel ('x');
%zlabel ('Ey');
%grid ;

%figure(6);
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
%surf(X, Y, Ez*fieldFactor);
%shading interp;
%title ('Ez');
%xlabel ('y');
%ylabel ('x');
%zlabel ('Ez');
%grid ;

figure(7);
hold on;
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx/2)*rho);
%surf(X, Y, Bnorm/B0);
%contourf(X,Y,Bnorm);
imagesc(1:Ny,1:Nx/2,Bnorm);
plot(1:Ny,(1:Ny)*10,'red');
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