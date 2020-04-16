clear;
shock_x = importdata('shock_x.dat');

N = size(shock_x,1);
dt = 1;
dx = 1;

mp = 1.67262177E-24;
q = 4.84*10^-10;
me = mp/100;
c = 2.99792458E10;
n = 1;
ntristan = 2;
sigma = 0.04;
gamma = 5.0;
ctristan = 0.45;
comp = 5;
omp = ctristan/comp;
qtristan = omp*omp*gamma/(ntristan*(1 + me/mp));
metristan = qtristan;
fieldScale = sqrt(4*3.14*(n/ntristan)*(me/metristan)*(c*c/(ctristan*ctristan)));
samplingFactor = 20;
timestep = 10000;

omega = sqrt(4*3.14*n*q*q/(gamma*me));
dt = ctristan/(comp*omega);
dx = samplingFactor*c*dt/ctristan;
dt = dt*timestep;

shock_v(1:N) = 0;

for i = 2:N,
    shock_v(i) = ((shock_x(i,2) - shock_x(i-1,2))/(shock_x(i,1) - shock_x(i-1,1)))*dx/dt;
end;
shock_v(1) = shock_v(2);

figure(1);
plot (shock_x(1:N,1)*dt, shock_x(1:N,2)*dx, 'red');
title ('shock');
xlabel ('t');
ylabel ('x');
grid ;

figure(2);
plot (shock_x(1:N,1)*dt, shock_v(1:N)/c, 'red');
title ('shock');
xlabel ('t');
ylabel ('beta');
grid ;