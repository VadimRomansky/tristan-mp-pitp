clear;
load energy

N1 = 1;
N2 = size(energy,1);

figure(1);
%plot (general(1:N2,3)*omega_plasma, general(1:N2,4), 'green', general(1:N2,3)*omega_plasma, general(1:N2,5), 'blue', general(1:N2, 3)*omega_plasma, general(1:N2,6), 'black', general(1:N2, 3)*omega_plasma, general(1:N2,7), 'red', general(1:N2, 3)*omega_plasma, general(1:N2,11), 'yellow');
plot (energy(1:N2,1), energy(1:N2,4), 'green', energy(1:N2,1), energy(1:N2,2), 'blue', energy(1:N2,1), energy(1:N2,3), 'black', energy(1:N2,1), energy(1:N2,5), 'red');

xlabel ('{{t {\omega}_p}}');
ylabel ('E/E_0');

%legend('particle', 'electric','magnetic', 'full', 'theoretical','Location','southwest');
legend('particle', 'electric','magnetic', 'full','Location','southwest');
grid ;