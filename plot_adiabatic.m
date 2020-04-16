clear;
directory_name = './output4/';
file_name = 'spect0';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

Nx = size(fp,1);
Np = size(fp,2);



me = 0.91*10^-27;
mass_ratio = 100;
mp = me*mass_ratio;
gam = 1.5;
beta = sqrt(1 - 1/(gam*gam));
c = 2.99792458*10^10;
Te = 5*10^9;
Tp = 3.5*10^10;
Pekappa = 14*me*c;
Ppkappa = mp*c;
kappa = 4;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;




norm = 1;



Prese(1:Nx) = 0;
Presp(1:Nx) = 0;
Pres(1:Nx) = 0;
rhoe(1:Nx) = 0;
rhop(1:Nx) = 0;
rho(1:Nx) = 0;

for i = 1:Nx,
    for j = 2:Np,
        u = sqrt((g(j)+1)*(g(j)+1) - 1);
        v = c*u/sqrt(1 + u*u);
        Prese(i) = Prese(i) + me*c*v*u*fe(i,j)*g(j)/3;
        rhoe(i) = rhoe(i) + me*fe(i,j)*g(j);
        
        Presp(i) = Presp(i) + mp*c*v*u*fp(i,j)*g(j)/3;
        rhop(i) = rhop(i) + mp*fp(i,j)*g(j);
    end;
    Pres(i) = Prese(i) + Presp(i);
    rho(i) = rhoe(i) + rhop(i);
end;
Prese1 = smooth(Prese,0.01,'loess');
rhoe1 = smooth(rhoe,0.01,'loess');
Presp1 = smooth(Presp,0.01,'loess');
rhop1 = smooth(rhop,0.01,'loess');
Pres1 = smooth(Pres,0.01,'loess');
rho1 = smooth(rho,0.01,'loess');

adiabatice(1:Nx) = 0;
adiabaticp(1:Nx) = 0;
adiabatic(1:Nx) = 0;
index = 1000;
for i = 2:Nx,
    adiabatice(i) = log(Prese1(i)/Prese1(i-1))/log(rhoe1(i-1)/rhoe1(i));
    adiabaticp(i) = log(Presp1(i)/Presp1(i-1))/log(rhop1(i-1)/rhop1(i));
    adiabatic(i) = log(Pres1(i)/Pres1(i-1))/log(rho1(i-1)/rho1(i));
end;

figure(1);
plot(1:Nx,adiabatice(1:Nx),'color','red');
title ('\gamma e');
xlabel ('x');
ylabel ('\gamma e');
grid ;

figure(2);
plot(1:Nx,adiabaticp(1:Nx),'color','red');
title ('\gamma p');
xlabel ('x');
ylabel ('\gamma p');
grid ;

figure(3);
plot(1:Nx,adiabatic(1:Nx),'color','red');
title ('\gamma');
xlabel ('x');
ylabel ('\gamma');
grid ;

figure(4);
plot(1:Nx,rhoe1(1:Nx),'color','red');
title ('\rho e');
xlabel ('x');
ylabel ('\rho e');
grid ;

figure(5);
plot(1:Nx,rhop1(1:Nx),'color','red');
title ('\rho p');
xlabel ('x');
ylabel ('\rho p');
grid ;

figure(6);
plot(1:Nx,rho1(1:Nx),'color','red');
title ('\rho');
xlabel ('x');
ylabel ('\rho');
grid ;

figure(7);
plot(1:Nx,Prese1(1:Nx),'color','red');
title ('p e');
xlabel ('x');
ylabel ('p e');
grid ;

figure(8);
plot(1:Nx,Presp1(1:Nx),'color','red');
title ('p p');
xlabel ('x');
ylabel ('p p');
grid ;

figure(9);
plot(1:Nx,Pres1(1:Nx),'color','red');
title ('p');
xlabel ('x');
ylabel ('p');
grid ;