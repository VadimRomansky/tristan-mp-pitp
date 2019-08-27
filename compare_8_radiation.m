clear;

load radiation0.dat;
load radiation1.dat;
load radiation2.dat;
load radiation3.dat;
load radiation4.dat;
load radiation5.dat;
load radiation6.dat;
load radiation7.dat;

N0 = size(radiation0,1);
N1 = size(radiation1,1);
N2 = size(radiation2,1);
N3 = size(radiation3,1);
N4 = size(radiation4,1);
N5 = size(radiation5,1);
N6 = size(radiation6,1);
N7 = size(radiation7,1);

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('I_{\nu} erg/cm^{-3} sr');

plot(radiation0(1:N0,1),radiation0(1:N0,4),'red',radiation1(1:N1,1),radiation1(1:N1,4),'blue',radiation2(1:N2,1),radiation2(1:N2,4),'green',radiation3(1:N3,1),radiation3(1:N3,4),'black');
plot(radiation4(1:N4,1),radiation4(1:N4,4),'red',radiation5(1:N5,1),radiation5(1:N5,4),'blue',radiation6(1:N6,1),radiation6(1:N6,4),'green',radiation7(1:N7,1),radiation7(1:N7,4),'black');

legend('April','May','June','August','April no turb','May no turb','June no turb','August no turb');

grid ;