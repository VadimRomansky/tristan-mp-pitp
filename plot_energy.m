clear;
load energy

N1 = 1;
N2 = size(energy,1);

Nskinlength = 5;

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

energyFactor = rho*rho/tau;

tau = c1/Nskinlength;

figure(1);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (energy(1:N2,1)*tau, energy(1:N2,4)*energyFactor, 'green', energy(1:N2,1)*tau, energy(1:N2,2)*energyFactor, 'blue', energy(1:N2,1)*tau, energy(1:N2,3)*energyFactor, 'black', energy(1:N2,1)*tau, energy(1:N2,5)*energyFactor, 'red');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('particle', 'electric','magnetic', 'full','Location','southwest');
grid ;