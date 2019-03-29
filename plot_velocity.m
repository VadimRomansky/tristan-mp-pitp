clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.020';
full_name = strcat(directory_name, file_name, file_number);
fileinfo = hdf5info(full_name);
Upx = hdf5read(full_name,'v4xi');
Uex = hdf5read(full_name,'v4x');
Upy = hdf5read(full_name,'v4yi');
Uey = hdf5read(full_name,'v4y');
Upz = hdf5read(full_name,'v4zi');
Uez = hdf5read(full_name,'v4z');

Nx = size(Upx, 1);
Ny = size(Upx, 2);

Vpx(1:Nx,1:Ny) = 0;
Vex(1:Nx,1:Ny) = 0;
Vpy(1:Nx,1:Ny) = 0;
Vey(1:Nx,1:Ny) = 0;
Vpz(1:Nx,1:Ny) = 0;
Vez(1:Nx,1:Ny) = 0;

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
samplingFactor = 5;
tau = c1*rho/c0;
rho =0.1;

%densityFactor = 1.0/(rho*rho*rho);
rho = rho*samplingFactor;

ypoint = fix(Ny/2)+1;

for i = 1:Nx,
    for j = 1:Ny,
        g = sqrt(1 + Upx(i,j)*Upx(i,j) + Upy(i,j)*Upy(i,j) + Upz(i,j)*Upz(i,j));
        Vpx(i,j) = Upx(i,j)/g;
        Vpy(i,j) = Upy(i,j)/g;
        Vpz(i,j) = Upz(i,j)/g;
        
        g = sqrt(1 + Uex(i,j)*Uex(i,j) + Uey(i,j)*Uey(i,j) + Uez(i,j)*Uez(i,j));
        Vex(i,j) = Uex(i,j)/g;
        Vey(i,j) = Uey(i,j)/g;
        Vez(i,j) = Uez(i,j)/g;
    end;
end;

figure(1);
plot ((1:Nx)*rho,Vpx(1:Nx, fix(Ny/2)+1), 'red');
title ('Vp');
xlabel ('x');
ylabel ('Vp');
grid ;

figure(2);
plot ((1:Nx)*rho,Vex(1:Nx, fix(Ny/2)+1), 'red');
title ('Ve');
xlabel ('x');
ylabel ('Ve');
grid ;

dlmwrite('Vpx.dat',Vpx(:,ypoint));
dlmwrite('Vpy.dat',Vpy(:,ypoint));
dlmwrite('Vpz.dat',Vpz(:,ypoint));
dlmwrite('Vex.dat',Vex(:,ypoint));
dlmwrite('Vey.dat',Vey(:,ypoint));
dlmwrite('Vez.dat',Vez(:,ypoint));