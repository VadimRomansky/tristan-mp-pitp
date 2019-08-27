clear;

load radiation0.dat;
load radiation1.dat;
load radiation2.dat;

N0 = size(radiation0,1);
N1 = size(radiation1,1);
N2 = size(radiation2,1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu}');
ylabel ('I_{\nu} erg/cm^{-3} sr');

plot(radiation0(1:N0,1),radiation0(1:N0,2),'red',radiation1(1:N1,1),radiation1(1:N1,2),'blue',radiation2(1:N2,1),radiation2(1:N2,2),'green');

legend('mp/me = 25','mp/me = 50','mp/me = 100');

grid ;