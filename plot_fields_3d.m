clear;
directory_name = './output3/';
file_name = 'flds3.tot';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
Bx = double(hdf5read(full_name,'bx'));
By = double(hdf5read(full_name,'by'));
Bz = double(hdf5read(full_name,'bz'));
Ex = double(hdf5read(full_name,'ex'));
Ey = double(hdf5read(full_name,'ey'));
Ez = double(hdf5read(full_name,'ez'));
mp = 1.67262177E-24;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 4.0;
gamma = 1.5;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));
B1 = sqrt(4*3.14*gamma*n*(1 + me/mp)*me*c*c*sigma)/fieldScale;
%Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)
B0 = Bz(1,1);
samplingFactor = 20;

rho = samplingFactor;

startx = 5;

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

Nxm = fix(Nx/15);

miniBx(1:Nxm, 1:Ny) = 0;
miniBy(1:Nxm, 1:Ny) = 0;
miniBz(1:Nxm, 1:Ny) = 0;
miniEx(1:Nxm, 1:Ny) = 0;
miniEy(1:Nxm, 1:Ny) = 0;
miniEz(1:Nxm, 1:Ny) = 0;

for i = 1:Nxm,
    for j = 1:Ny,
        miniBx(i,j) = Bx(i,j);
        miniBy(i,j) = By(i,j);
        miniBz(i,j) = Bz(i,j);
        miniEx(i,j) = Ex(i,j);
        miniEy(i,j) = Ey(i,j);
        miniEz(i,j) = Ez(i,j);
    end;
end;

for i=1:Nx,
    for j = 1:Ny,
        Bnorm(i,j) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        Ediff(i,j) = sqrt((Ex(i,j))*(Ex(i,j)) + (Ey(i,j) +beta*Bz(i,j))*(Ey(i,j) +beta*Bz(i,j)) + (Ez(i,j) - beta*By(i,j))*(Ez(i,j) - beta*By(i,j)))/Bnorm(i,j);
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

for i=1:Ny,
    for j = 1:Ny,
        fourierBx(i,j) = Bx(i+startx,j);
        %fourierBx(i,j) = sin((2*pi/Ny)*(i + j)) + sin((2*pi/Ny)*(i)) + sin((2*pi/Ny)*j) + sin((2*pi/Ny)*(2*i + 2*j));
        fourierBy(i,j) = By(i+startx,j);
        fourierBz(i,j) = Bz(i+startx,j);
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

miniFourierBx(1:Ny/2,1:Ny/2) = 0;
miniFourierBy(1:Ny/2,1:Ny/2) = 0;
miniFourierBz(1:Ny/2,1:Ny/2) = 0;
for i=1:Ny/2,
    for j = 1:Ny/2,
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

Bc2 = miniFourierBx(1,1)*miniFourierBx(1,1)+miniFourierBy(1,1)*miniFourierBy(1,1)+miniFourierBz(1,1)*miniFourierBz(1,1);
Btot2 = 0;
for i = 1:Ny/2,
    for j = 1:Ny/2,
        factor = 0.5;
        if ((i == 1) && (j == 1))
            factor = 1.0;
        end;
        Btot2 = Btot2 + factor*miniFourierBx(i,j)*miniFourierBx(i,j)+miniFourierBy(i,j)*miniFourierBy(i,j)+miniFourierBz(i,j)*miniFourierBz(i,j);
    end;
end;
modeFraction = (Btot2 - Bc2)/Btot2;

figure(1);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniBx);
shading interp;
title ('Bx');
xlabel ('y');
ylabel ('x');
zlabel ('Bx');
grid ;

figure(2);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniBy);
shading interp;
title ('By');
xlabel ('y');
ylabel ('x');
zlabel ('By');
grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniBz);
shading interp;
title ('Bz');
xlabel ('y');
ylabel ('x');
zlabel ('Bz');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniEx);
shading interp;
title ('Ex');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
%grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniEy);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (1:Nxm));
surf(X, Y, miniEz);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;

figure(7);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:Ny), (1:Nx));
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
[X, Y] = meshgrid((1:Ny), (1:Nx));
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
[X, Y] = meshgrid((1:Ny/2), (1:Ny/2));
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
[X, Y] = meshgrid((1:Ny/2), (1:Ny/2));
surf(X, Y, miniFourierBy);
shading interp;
title ('By(k)');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('By');
grid ;

figure(12);
colormap Jet;
caxis ([0 8])
[X, Y] = meshgrid((1:Ny/2), (1:Ny/2));
surf(X, Y, miniFourierBz);
shading interp;
title ('Bz(k)');
xlabel ('y \omega /c');
ylabel ('x \omega /c');
zlabel ('Bz');
grid ;

dlmwrite('Bx.dat',Bx,'delimiter',' ');
dlmwrite('By.dat',By,'delimiter',' ');
dlmwrite('Bz.dat',Bz,'delimiter',' ');
dlmwrite('Ex.dat',Ex,'delimiter',' ');
dlmwrite('Ey.dat',Ey,'delimiter',' ');
dlmwrite('Ez.dat',Ez,'delimiter',' ');