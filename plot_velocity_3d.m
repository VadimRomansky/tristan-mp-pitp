clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.020';
full_name = strcat(directory_name, file_name, file_number);
info = h5info(full_name);
vxi = double(hdf5read(full_name,'v4xi'));
vyi = double(hdf5read(full_name,'v4yi'));
vzi = double(hdf5read(full_name,'v4zi'));
vxe = double(hdf5read(full_name,'v4x'));
vye = double(hdf5read(full_name,'v4y'));
vze = double(hdf5read(full_name,'v4z'));

Bx = double(hdf5read(full_name,'bx'));
By = double(hdf5read(full_name,'by'));
Bz = double(hdf5read(full_name,'bz'));
Ex = double(hdf5read(full_name,'ex'));
Ey = double(hdf5read(full_name,'ey'));
Ez = double(hdf5read(full_name,'ez'));

np = hdf5read(full_name,'densi');
ne = hdf5read(full_name,'dens');

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

samplingFactor = 20;

rho = samplingFactor;

startx = 600;

Nx = size(vxi, 1);
Ny = size(vyi, 2);

for i = 1:Nx,
    for j = 1:Ny,
        u = sqrt(vxi(i,j)*vxi(i,j) + vyi(i,j)*vyi(i,j)+vzi(i,j)*vzi(i,j));
        g = sqrt(1+u*u);
        vxi(i,j) = vxi(i,j)/g;
        vyi(i,j) = vyi(i,j)/g;
        vzi(i,j) = vzi(i,j)/g;
        
        u = sqrt(vxe(i,j)*vxe(i,j) + vye(i,j)*vye(i,j)+vze(i,j)*vze(i,j));
        g = sqrt(1+u*u);
        vxe(i,j) = vxe(i,j)/g;
        vye(i,j) = vye(i,j)/g;
        vze(i,j) = vze(i,j)/g;
    end;
end;

%Nxm = fix(Nx/15);
Nxm = Ny;

minivxi(1:Nxm, 1:Ny) = 0;
minivyi(1:Nxm, 1:Ny) = 0;
minivzi(1:Nxm, 1:Ny) = 0;
minivxe(1:Nxm, 1:Ny) = 0;
minivye(1:Nxm, 1:Ny) = 0;
minivze(1:Nxm, 1:Ny) = 0;

for i = 1:Nxm,
    for j = 1:Ny,
        minivxi(i,j) = vxi(startx + i,j);
        minivyi(i,j) = vyi(startx + i,j);
        minivzi(i,j) = vzi(startx + i,j);
        minivxe(i,j) = vxe(startx + i,j);
        minivye(i,j) = vye(startx + i,j);
        minivze(i,j) = vze(startx + i,j);
    end;
end;

miniBx(1:Nxm, 1:Ny) = 0;
miniBy(1:Nxm, 1:Ny) = 0;
miniBz(1:Nxm, 1:Ny) = 0;
miniEx(1:Nxm, 1:Ny) = 0;
miniEy(1:Nxm, 1:Ny) = 0;
miniEz(1:Nxm, 1:Ny) = 0;
tempEx(1:Nxm, 1:Ny) = 0;
tempEy(1:Nxm, 1:Ny) = 0;
tempEz(1:Nxm, 1:Ny) = 0;

for i = 1:Nxm,
    for j = 1:Ny,
        miniBx(i,j) = Bx(startx + i,j);
        miniBy(i,j) = By(startx + i,j);
        miniBz(i,j) = Bz(startx + i,j);
        miniEx(i,j) = Ex(startx + i,j);
        miniEy(i,j) = Ey(startx + i,j);
        miniEz(i,j) = Ez(startx + i,j);
        tempEx(i,j) = -(minivye(i,j)*miniBz(i,j) - minivze(i,j)*miniBy(i,j));
        tempEy(i,j) = -(minivze(i,j)*miniBx(i,j) - minivxe(i,j)*miniBz(i,j));
        tempEz(i,j) = -(minivxe(i,j)*miniBy(i,j) - minivye(i,j)*miniBx(i,j));
    end;
end;

minine(1:Nxm, 1:Ny) = 0;
mininp(1:Nxm, 1:Ny) = 0;
miniq(1:Nxm, 1:Ny) = 0;
minijx(1:Nxm, 1:Ny) = 0;
minijy(1:Nxm, 1:Ny) = 0;
minijz(1:Nxm, 1:Ny) = 0;
minij2(1:Nxm, 1:Ny) = 0;
totalq = 0;
for i=1:Nxm,
    for j=1:Ny,
        minine(i,j)=ne(startx + i,j) - np(startx + i,j);
        mininp(i,j)=np(startx + i,j);
        miniq(i,j) = mininp(i,j) - minine(i,j);
        totalq = totalq + miniq(i,j);
        minijx(i,j) = (mininp(i,j)*minivxi(i,j) - minine(i,j)*minivxe(i,j));
        minijy(i,j) = (mininp(i,j)*minivyi(i,j) - minine(i,j)*minivye(i,j));
        minijz(i,j) = (mininp(i,j)*minivzi(i,j) - minine(i,j)*minivze(i,j));
        minij2(i,j) = minijx(i,j)*minijx(i,j) + minijy(i,j)*minijy(i,j) + minijz(i,j)*minijz(i,j);
    end;
end;

for i=1:Nxm,
    for j=1:Ny,
        miniq(i,j) = miniq(i,j) - totalq/(Nxm*Ny);
    end;
end;

figure(1);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivxi);
shading interp;
title ('vxi');
xlabel ('y');
ylabel ('x');
zlabel ('vx/c');
grid ;

figure(2);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivyi);
shading interp;
title ('vyi');
xlabel ('y');
ylabel ('x');
zlabel ('vy/c');
grid ;

figure(3);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivzi);
shading interp;
title ('vzi');
xlabel ('y');
ylabel ('x');
zlabel ('vz/c');
grid ;

figure(4);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivxe);
shading interp;
title ('vxe');
xlabel ('y');
ylabel ('x');
zlabel ('vx/c');
%grid ;

figure(5);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivye);
shading interp;
title ('vye');
xlabel ('y');
ylabel ('x');
zlabel ('vy/c');
grid ;

figure(6);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, minivze);
shading interp;
title ('vze');
xlabel ('y');
ylabel ('x');
zlabel ('vz/c');
grid ;

figure(7);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
hold on;
%surf(X, Y, miniEx);
imagesc((1:Ny), (startx + 1:startx + Nxm), miniEx);
quiver(X, Y, miniBy, miniBx);
shading interp;
title ('Ex');
xlabel ('y');
ylabel ('x');
zlabel ('Ex');
grid ;

figure(8);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, miniEy);
shading interp;
title ('Ey');
xlabel ('y');
ylabel ('x');
zlabel ('Ey');
grid ;

figure(9);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
surf(X, Y, miniEz);
shading interp;
title ('Ez');
xlabel ('y');
ylabel ('x');
zlabel ('Ez');
grid ;

% figure(10);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, tempEx);
% shading interp;
% title ('Ex');
% xlabel ('y');
% ylabel ('x');
% zlabel ('Ex');
% grid ;
% 
% figure(11);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, tempEy);
% shading interp;
% title ('Ey');
% xlabel ('y');
% ylabel ('x');
% zlabel ('Ey');
% grid ;
% 
% figure(12);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, tempEz);
% shading interp;
% title ('Ez');
% xlabel ('y');
% ylabel ('x');
% zlabel ('Ez');
% grid ;

% figure(13);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% quiver(X, Y, miniBy, miniBx);
% shading interp;
% title ('Bxy');
% xlabel ('y');
% ylabel ('x');
% grid ;
% 
% figure(14);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% quiver(X, Y, miniEy, miniEx);
% shading interp;
% title ('Exy');
% xlabel ('y');
% ylabel ('x');
% grid ;
% 
% figure(15);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% quiver(X, Y, minivyi, minivxi);
% shading interp;
% title ('Vxyi');
% xlabel ('y');
% ylabel ('x');
% grid ;
% 
% figure(16);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% quiver(X, Y, minivye, minivxe);
% shading interp;
% title ('vxye');
% xlabel ('y');
% ylabel ('x');
% grid ;
% 
% figure(17);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, mininp);
% shading interp;
% title ('np');
% xlabel ('y');
% ylabel ('x');
% zlabel ('np');
% grid ;
% 
% figure(18);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, minine);
% shading interp;
% title ('ne');
% xlabel ('y');
% ylabel ('x');
% zlabel ('ne');
% grid ;
% 
% figure(19);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, miniq);
% shading interp;
% title ('q');
% xlabel ('y');
% ylabel ('x');
% zlabel ('q');
% grid ;
% 
% potential = poicalc(miniq);
% 
% figure(20);
% colormap Jet;
% [X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
% surf(X, Y, potential);
% shading interp;
% title ('phi');
% xlabel ('y');
% ylabel ('x');
% zlabel ('phi');
% grid ;

figure(21);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
hold on;
%surf(X, Y, minij2);

imagesc((1:Ny), (startx + 1:startx + Nxm), minij2);
quiver(X, Y, miniBy, miniBx);
shading interp;
title ('j^2');
xlabel ('y');
ylabel ('x');
zlabel ('j^2');
grid ;

figure(22);
colormap Jet;
[X, Y] = meshgrid((1:Ny), (startx + 1:startx + Nxm));
hold on;
%surf(X, Y, minijz);
imagesc((1:Ny), (startx + 1:startx + Nxm), minijz);
quiver(X, Y, miniBy, miniBx);
shading interp;
title ('jz');
xlabel ('y');
ylabel ('x');
zlabel ('jz');
grid ;


dlmwrite('vxi.dat',vxi,'delimiter',' ');
dlmwrite('vyi.dat',vyi,'delimiter',' ');
dlmwrite('vzi.dat',vzi,'delimiter',' ');
dlmwrite('vxe.dat',vxe,'delimiter',' ');
dlmwrite('vye.dat',vye,'delimiter',' ');
dlmwrite('vze.dat',vze,'delimiter',' ');