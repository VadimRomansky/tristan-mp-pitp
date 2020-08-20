clear;
directory_name = './output1/';
file_name = 'flds4.tot';
file_number = '.020';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');

Nx = size(Bx, 1);
Ny = size(By, 2);

mp = 1.67262177E-24;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 0.4;
gamma = 1.5;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));
B1 = sqrt(4*3.14*gamma*n*(1 + me/mp)*me*c*c*sigma)/fieldScale;
%Binit=sqrt(gamma0*ppc0*.5*c**2*(me*(1+me/mi))*sigma)
B0 = Bz(10,10);
samplingFactor = 20;

rho = samplingFactor;

Bnorm(1:Nx) = 0;
Bperp(1:Nx) = 0;
theta(1:Nx) = 0;
Bxa(1:Nx) = 0;
Bya(1:Nx) = 0;
Bza(1:Nx) = 0;
Exa(1:Nx) = 0;
Eya(1:Nx) = 0;
Eza(1:Nx) = 0;

for i=1:Nx,
    for j=1:Ny,
        Bxa(i) = Bxa(i) + Bx(i,j);
        Bya(i) = Bya(i) + By(i,j);
        Bza(i) = Bza(i) + Bz(i,j);
        Exa(i) = Exa(i) + Ex(i,j);
        Eya(i) = Eya(i) + Ey(i,j);
        Eza(i) = Eza(i) + Ez(i,j);
        
        Bnorm(i) = Bnorm(i) + Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j);
        Bperp(i) = Bperp(i) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j);
        theta(i) = theta(i) + acos(Bx(i,j)/sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j)))*180/pi;
    end;
    Bxa(i) = Bxa(i)/Ny;
    Bya(i) = Bya(i)/Ny;
    Bza(i) = Bza(i)/Ny;
    Exa(i) = Exa(i)/Ny;
    Eya(i) = Eya(i)/Ny;
    Eza(i) = Eza(i)/Ny;
    
    Bnorm(i) = sqrt(Bnorm(i)/Ny);
    Bperp(i) = sqrt(Bperp(i)/Ny);
    theta(i) = theta(i)/Ny;
    
end;

figure(1);
plot ((1:Nx),Bxa(1:Nx)/B1, 'red');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
plot ((1:Nx)*rho,Bya(1:Nx)/B1, 'red');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
plot ((1:Nx),Bza(1:Nx)/B1, 'red');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
hold on;
plot ((1:Nx/2),Exa(1:Nx/2)/B1, 'red');
ytemp(1:2) = 0;
ytemp(1) = -1;
ytemp(2) = 1;
xtemp(1:2) = 0;
xtemp(1) = 170;
xtemp(2) = 170;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 220;
xtemp(2) = 220;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 270;
xtemp(2) = 270;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 320;
xtemp(2) = 320;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 370;
xtemp(2) = 370;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 420;
xtemp(2) = 420;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 470;
xtemp(2) = 470;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 520;
xtemp(2) = 520;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 570;
xtemp(2) = 570;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 620;
xtemp(2) = 620;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 670;
xtemp(2) = 670;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 720;
xtemp(2) = 720;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 770;
xtemp(2) = 770;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 820;
xtemp(2) = 820;
plot(xtemp(1:2),ytemp(1:2),'black');
xtemp(1) = 870;
xtemp(2) = 870;
plot(xtemp(1:2),ytemp(1:2),'black');
title ('Ex');
xlabel ('x');
ylabel ('Ex/B0');
grid ;

figure(5);
plot ((1:Nx/2),Eya(1:Nx/2)/B1, 'red');
title ('Ey');
xlabel ('x');
ylabel ('Ey/B0');
grid ;

figure(6);
plot ((1:Nx/2),Eza(1:Nx/2)/B1, 'red');
title ('Ez');
xlabel ('x');
ylabel ('Ez/B0');
grid ;

figure(7);
plot ((1:Nx),Bnorm(1:Nx)/B1, 'red');
title ('Bnorm');
xlabel ('x');
ylabel ('Bnorm');
grid ;

figure(8);
plot ((1:Nx),Bperp(1:Nx)/B1, 'red');
title ('Bperp');
xlabel ('x');
ylabel ('Bperp');
grid ;

figure(9);
plot ((1:Nx),theta(1:Nx), 'red');
title ('theta');
xlabel ('x');
ylabel ('theta');
grid ;