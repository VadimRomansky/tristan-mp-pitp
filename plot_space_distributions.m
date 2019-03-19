clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.005';
full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

Nx = size(fp,1);
Np = size(fp,2);

startx = 1;
endx = Nx/4;

Fp(1:5,1:Nx)=0;
Fe(1:5,1:Nx)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

me = 0.91*10^-27;
mass_ratio = 16;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

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

samplingFactor = 5;

rho = rho*samplingFactor;
rho = samplingFactor/Nskinlength;

normp = 0;
norme = 0;
norm = 1;

for i = 1:Np,
    %Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    %Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    Pp(i) = sqrt((g(i)+1)^2 - 1);
    Pe(i) = sqrt((g(i)+1)^2 - 1);
    
end;
%p1 = fix(Np/2);
%p2 = fix(2*Np/3);
%p3 = fix(3*Np/4);
%p4 = fix(7*Np/8);
%p5 = fix(Np - 1);

pne(1:5) = 0;
pnp(1:5) = 0;

pnp(1) = 135;
pnp(2) = 140;
pnp(3) = 150;
pnp(4) = 155;
pnp(5) = 160;

pne(1) = 130;
pne(2) = 140;
pne(3) = 150;
pne(4) = 155;
pne(5) = 160;

for i = 1:Nx,
    for j = 1:5,
        for k = -2:2,
            Fp(j,i)=Fp(j,i) + fp(i,pnp(j)+k)*(Pp(pnp(j)+k))/(1+g(pnp(j)+k))*(Pp(pnp(j) + k + 1) - Pp(pnp(j) + k));
            Fe(j,i)=Fe(j,i) + fe(i,pne(j)+k)*(Pe(pne(j)+k))/(1+g(pne(j)+k))*(Pe(pne(j) + k + 1) - Pe(pne(j) + k));
        end;
    end;
end;




figure(1);
%plot ((1:Nx)*rho,Fp(1,1:Nx), 'red',(1:Nx)*rho,Fp(2,1:Nx), 'green',(1:Nx)*rho,Fp(3,1:Nx), 'blue',(1:Nx)*rho,Fp(4,1:Nx), 'black',(1:Nx)*rho,Fp(5,1:Nx), 'magenta');
plot ((1:Nx)*rho,Fp(1,1:Nx), 'red',(1:Nx)*rho,Fp(2,1:Nx), 'green',(1:Nx)*rho,smooth(Fp(3,:),11), 'blue',(1:Nx)*rho,smooth(Fp(4,:),11), 'black',(1:Nx)*rho,smooth(Fp(5,:),11), 'magenta');
title ('F_p');
xlabel ('x');
ylabel ('Fp*p^4');
legend(strcat('p/mc = ',num2str(Pp(pnp(1)))), strcat('p/mc = ',num2str(Pp(pnp(2)))),strcat('p/mc = ',num2str(Pp(pnp(3)))),strcat('p/mc = ',num2str(Pp(pnp(4)))),strcat('p/mc = ',num2str(Pp(pnp(5)))),'Location','southeast');
grid ;

figure(2);
plot ((1:Nx)*rho,Fe(1,1:Nx), 'red',(1:Nx)*rho,Fe(2,1:Nx), 'green',(1:Nx)*rho,Fe(3,1:Nx), 'blue',(1:Nx)*rho,Fe(4,1:Nx), 'black',(1:Nx)*rho,Fe(5,1:Nx), 'magenta');
%plot (Pe(1:Np),Fe(1:Np), 'red',Pe(1:Np), Fejuttner(1:Np), 'blue');
title ('F_e');
xlabel ('x');
ylabel ('F_e*p^4');
legend(strcat('p/mc = ',num2str(Pe(pnp(1)))), strcat('p/mc = ',num2str(Pe(pnp(2)))),strcat('p/mc = ',num2str(Pe(pnp(3)))),strcat('p/mc = ',num2str(Pe(pnp(4)))),strcat('p/mc = ',num2str(Pe(pnp(5)))),'Location','southeast');
grid ;