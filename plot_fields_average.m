clear;
directory_name = './output/';
file_name = 'flds.tot';
file_number = '.060';
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
plot ((1:Nx)*rho,Bxa(1:Nx)*fieldScale, 'red');
title ('Bx');
xlabel ('x');
ylabel ('Bx');
grid ;

figure(2);
plot ((1:Nx)*rho,Bya(1:Nx)*fieldScale, 'red');
title ('By');
xlabel ('x');
ylabel ('By');
grid ;

figure(3);
plot ((1:Nx)*rho,Bza(1:Nx)*fieldScale, 'red');
title ('Bz');
xlabel ('x');
ylabel ('Bz');
grid ;

figure(4);
plot ((1:Nx)*rho,Exa(1:Nx)*fieldScale, 'red');
title ('Ex');
xlabel ('x');
ylabel ('Ex');
grid ;

figure(5);
plot ((1:Nx)*rho,Eya(1:Nx)*fieldScale, 'red');
title ('Ey');
xlabel ('x');
ylabel ('Ey');
grid ;

figure(6);
plot ((1:Nx)*rho,Eza(1:Nx)*fieldScale, 'red');
title ('Ez');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(7);
plot ((1:Nx)*rho,Bnorm(1:Nx)*fieldScale, 'red');
title ('Bnorm');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(8);
plot ((1:Nx)*rho,Bperp(1:Nx)*fieldScale, 'red');
title ('Bperp');
xlabel ('x');
ylabel ('Ez');
grid ;

figure(9);
plot ((1:Nx)*rho,theta(1:Nx), 'red');
title ('theta');
xlabel ('x');
ylabel ('theta');
grid ;