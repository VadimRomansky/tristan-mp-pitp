clear;
shock_x = importdata('mega_shock_x.dat');

c = 2.99792458*10^10;
mp = 1.67*10^-24;
kB = 1.3806488*10^-16;
Nt = size(shock_x,1);
Nn = size(shock_x,2);
dt = 1;
dx = 1;


figure(1);
hold on;
plot (shock_x(1:Nt,1), shock_x(1:Nt,2), 'red');
plot (shock_x(1:Nt,3), shock_x(1:Nt,4), 'blue');
%plot (shock_x(1:Nt,5), shock_x(1:Nt,6), 'red');
plot (shock_x(1:Nt,7), shock_x(1:Nt,8), 'black');
plot (shock_x(1:Nt,9), shock_x(1:Nt,10), 'green');
%plot (shock_x(1:Nt,11), shock_x(1:Nt,12), 'red');
%plot (shock_x(1:Nt,13), shock_x(1:Nt,14), 'red');
plot (shock_x(1:Nt,15), shock_x(1:Nt,16), 'magenta');
%plot (shock_x(1:Nt,17), shock_x(1:Nt,18), 'red');
%plot (shock_x(1:Nt,19), shock_x(1:Nt,20), 'red');
title ('shock');
xlabel ('t s');
ylabel ('x cm');
%legend('reg1','reg2','reg3','turb1','turb2','turb3','turb4','gam1','gam2','gam3');
legend('reg1','reg2','turb1','turb2','gam1');
grid ;
beta(1:Nn/2) = 0;
for i = 1:Nn/2,
    beta(i) = shock_x(Nt,2*i)/(c*shock_x(Nt,2*i-1));
end;
sigma(1:Nn/2) = 0;
sigma(1) = 0.04;
sigma(2) = 0.004;
sigma(3) = 0.0004;
sigma(4) = 0.04;
sigma(5) = 0.04;
sigma(6) = 0.004;
sigma(7) = 0.004;
sigma(8) = 0.04;
sigma(9) = 0.04;
sigma(10) = 0.04;
gamma(1:Nn/2) = 0;
gamma(1) = 1.5;
gamma(2) = 1.5;
gamma(3) = 1.5;
gamma(4) = 1.5;
gamma(5) = 1.5;
gamma(6) = 1.5;
gamma(7) = 1.5;
gamma(8) = 2.0;
gamma(9) = 5.0;
gamma(10) = 10.0;

adiab(1:Nn/2) = 0;
temp(1:Nn/2) = 0;
for i = 1:Nn/2,
    adiab(i) = 2*(beta(i)*beta(i)*sigma(i) - sigma(i) + beta(i) + beta(i)*beta(i))/(2*beta(i) + beta(i)*sigma(i) - sigma(i));
    temp(i) = (mp*c*c/kB)*gamma(i)*beta(i)*(1.0 - (sigma(i)/2)*(1-beta(i))/2);
end;
