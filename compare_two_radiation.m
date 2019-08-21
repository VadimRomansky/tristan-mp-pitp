clear;

load radiation0.dat;
load radiation1.dat;

N0 = size(radiation0,1);
N1 = size(radiation1,1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu}');
ylabel ('I_{\nu}');

plot(radiation0(1:N0,1),radiation0(1:N0,2),'red',radiation1(1:N1,1),radiation1(1:N1,2)*1000,'blue');

legend('10^{16}', '10^{17}','Location','southeast');

grid ;