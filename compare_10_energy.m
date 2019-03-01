clear;
directory_name = './output10/';
file_name = 'energy';
%file_number = '.010';
Nd = 10;
start = 0;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
%LegendTitle = {'t*{\Omega} = 30','t*{\Omega} = 60','t*{\Omega} = 90', 't*{\Omega} = 120', 't*{\Omega} = 150','t*{\Omega} = 180'};
LegendTitle = {'noturb B in plane', 'noturb B out plane','turb 0.5 iso B in plane', 'turb 0.5 iso B out plane', 'turb 0.5 aniso B in plane', 'turb 0.5 aniso B out plane','turb 0.9 iso B in plane', 'turb 0.9 iso B out plane', 'turb 0.9 aniso B in plane', 'turb 0.9 aniso B out plane'};

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 25;
%mp = 1.67262*10^-24;
%me = mp/mass_ratio;
me = 9.1*10^-28;
mp = me*mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;
gamma = 1.5;

omega = sqrt(4*pi*n*q*q/(me*gamma));

rho = c0/(omega*Nskinlength);
c1=0.45;

tau = c1*rho/c0;

energyFactor = rho*rho/tau;

tau = c1/Nskinlength;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
set(0, 'DefaultLineLineWidth', 1);

figure(1);
hold on;
for i = 1:Nd,   
    full_name = strcat(directory_name, file_name, num2str(start + i-1));
    energy = importdata(full_name);
    N = size(energy,1);
    
    plot (energy(1:N,1)*tau, energy(1:N,5)*energyFactor,'color',Color{i});
end;

xlabel ('{{t {\omega}_p}}');
ylabel ('full E/E_0');

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest')
grid ;

figure(2);
hold on;
for i = 1:Nd,   
    full_name = strcat(directory_name, file_name, num2str(start + i-1));
    energy = importdata(full_name);
    N = size(energy,1);
    
    plot (energy(1:N,1)*tau, energy(1:N,4)*energyFactor,'color',Color{i});
end;

xlabel ('{{t {\omega}_p}}');
ylabel ('particle E/E_0');

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest')
grid ;

figure(3);
hold on;
for i = 1:Nd,   
    full_name = strcat(directory_name, file_name, num2str(start + i-1));
    energy = importdata(full_name);
    N = size(energy,1);
    
    plot (energy(1:N,1)*tau, energy(1:N,3)*energyFactor,'color',Color{i});
end;

xlabel ('{{t {\omega}_p}}');
ylabel ('magnetic E/E_0');

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest')
grid ;

figure(4);
hold on;
for i = 1:Nd,   
    full_name = strcat(directory_name, file_name, num2str(start + i-1));
    energy = importdata(full_name);
    N = size(energy,1);
    
    plot (energy(1:N,1)*tau, energy(1:N,2)*energyFactor,'color',Color{i});
end;

xlabel ('{{t {\omega}_p}}');
ylabel ('electric E/E_0');

legend(LegendTitle{1}, LegendTitle{2}, LegendTitle{3}, LegendTitle{4}, LegendTitle{5}, LegendTitle{6}, LegendTitle{7}, LegendTitle{8}, LegendTitle{9}, LegendTitle{10},'Location','northwest')
grid ;