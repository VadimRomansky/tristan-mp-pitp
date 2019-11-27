clear;
directory_name = './output5/';
file_name = 'flds.tot';
file_number = '.000';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');
mp = 1.67262177E-24;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 4;
sigma = 4.0;
gamma = 1.5;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/1)*(c*c/(0.45*0.45)));
B1 = sqrt(gamma*n*(1 + me/mp)*me*c*c*sigma)/fieldScale;
B0 = 0.03030750;
Nx = size(Bx, 1);
Ny = size(By, 2);

Bnorm(1:Nx, 1:Ny) = 0;
Ediff(1:Nx, 1:Ny) = 0;

fourierBx(1:Ny, 1:Ny) = 0;
fourierBy(1:Ny, 1:Ny) = 0;
fourierBz(1:Ny, 1:Ny) = 0;

fourierBx2(1:Ny, 1:Ny) = 0;
fourierBy2(1:Ny, 1:Ny) = 0;
fourierBz2(1:Ny, 1:Ny) = 0;
%Bperp(1:Nx, 1:Ny) = 0;
g = 1.5;
beta = sqrt(1-1/(g*g));

for i=1:Nx,
    for j = 1:Ny,
        Bnorm(i,j) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        Ediff(i,j) = sqrt((Ex(i,j))*(Ex(i,j)) + (Ey(i,j) +beta*Bz(i,j))*(Ey(i,j) +beta*Bz(i,j)) + (Ez(i,j) - beta*By(i,j))*(Ez(i,j) - beta*By(i,j)))/Bnorm(i,j);
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

for i=1:Ny,
    for j = 1:Ny,
        fourierBx(i,j) = Bx(i+100,j);
        %fourierBx(i,j) = sin((2*pi/Ny)*(i + j)) + sin((2*pi/Ny)*(i)) + sin((2*pi/Ny)*j) + sin((2*pi/Ny)*(2*i + 2*j));
        fourierBy(i,j) = By(i+100,j);
        fourierBz(i,j) = Bz(i+100,j);
    end;
end;

fourierBx2 = fft2(fourierBx);
fourierBy2 = fft2(fourierBy);
fourierBz2 = fft2(fourierBz);

for i=1:Ny,
    for j = 1:Ny,
        fourierBx(i,j) = abs(fourierBx2(i,j));
        fourierBy(i,j) = abs(fourierBy2(i,j));
        fourierBz(i,j) = abs(fourierBz2(i,j));
    end;
end;

miniFourierBx(1:10,1:10) = 0;
miniFourierBy(1:10,1:10) = 0;
miniFourierBz(1:10,1:10) = 0;
for i=1:10,
    for j = 1:10,
        miniFourierBx(i,j) = abs(fourierBx2(i,j));
        miniFourierBy(i,j) = abs(fourierBy2(i,j));
        miniFourierBz(i,j) = abs(fourierBz2(i,j));
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
surf(X, Y, Bx);
shading interp;
title ('Bx');
xlabel ('y');
ylabel ('x');
zlabel ('Bx');
grid ;

figure(2);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, By);
shading interp;
title ('By');
xlabel ('y');
ylabel ('x');
zlabel ('By');
grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Bz);
shading interp;
title ('Bz');
xlabel ('y');
ylabel ('x');
zlabel ('Bz');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ex);
shading interp;
title ('Ex');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
%grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ey);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ez);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;

figure(7);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
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
[X, Y] = meshgrid((1:Ny)*rho, (1:Nx)*rho);
surf(X, Y, Ediff);
shading interp;
title ('(E-vxB)/B');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('E');
grid ;

figure(10);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:10)*rho, (1:10)*rho);
surf(X, Y, miniFourierBx);
shading interp;
title ('Bx(k)');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('Bx');
grid ;

figure(11);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:10)*rho, (1:10)*rho);
surf(X, Y, miniFourierBy);
shading interp;
title ('Bx(k)');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('By');
grid ;

figure(12);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:10)*rho, (1:10)*rho);
surf(X, Y, miniFourierBz);
shading interp;
title ('Bx(k)');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('Bz');
grid ;

dlmwrite('Bx.dat',Bx,'delimiter',' ');
dlmwrite('By.dat',By,'delimiter',' ');
dlmwrite('Bz.dat',Bz,'delimiter',' ');
dlmwrite('Ex.dat',Ex,'delimiter',' ');
dlmwrite('Ey.dat',Ey,'delimiter',' ');
dlmwrite('Ez.dat',Ex,'delimiter',' ');