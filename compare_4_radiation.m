clear;

load radiation0.dat;
load radiation1.dat;
load radiation2.dat;
load radiation3.dat;

N0 = size(radiation0,1);
N1 = size(radiation1,1);
N2 = size(radiation2,1);
N3 = size(radiation3,1);

augx(1:5) = 0;
augx(1) = 0.335;
augx(2) = 0.625;
augx(3) = 1.46;
augx(4) = 4.92;
augx(5) = 8.57;
augy(1:5) = 0;
augy(1) = 3.29;
augy(2) = 7.77;
augy(3) = 8.53;
augy(4) = 2.42;
augy(5) = 1.06;

augmax = 0.886;
augmaxy = 11.2;

junx(1:4) = 0;
junx(1) = 0.628;
junx(2) = 1.45;
junx(3) = 4.89;
junx(4) = 8.53;
juny(1:4) = 0;
juny(1) = 2.98;
juny(2) = 12.3;
juny(3) = 5.79;
juny(4) = 3.15;

junmaxx = 1.65;
junmaxy = 13.2;

mayx(1:3) = 0;
mayx(1) = 1.46;
mayx(2) = 4.94;
mayx(3) = 8.62;
mayy(1:3) = 0;
mayy(1) = 4.91;
mayy(2) = 12.0;
mayy(3) = 6.67;

maymaxx = 2.96;
maymaxy = 15.2;

aprx(1:4) = 0;
aprx(1) = 1.44;
aprx(2) = 4.91;
aprx(3) = 8.58;
aprx(4) = 22.8;
apry(1:4) = 0;
apry(1) = 0.993;
apry(2) = 13.9;
apry(3) = 17.1;
apry(4) = 5.11;

aprmaxx = 6.50;
aprmaxy = 19.3;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('I_{\nu} erg/Hz cm^{3} sr');

plot(radiation0(1:N0,1),radiation0(1:N0,4),'red',radiation1(1:N1,1),radiation1(1:N1,4),'blue',radiation2(1:N2,1),radiation2(1:N2,4),'green',radiation3(1:N3,1),radiation3(1:N3,4),'black');

legend('April','May','June','August');

grid;

figure(2);
hold on;
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('mJy');

plot(radiation0(1:N0,1),radiation0(1:N0,9),'red',radiation1(1:N1,1),radiation1(1:N1,9),'blue',radiation2(1:N2,1),radiation2(1:N2,9),'green',radiation3(1:N3,1),radiation3(1:N3,9),'black');
plot(aprx(1:4), apry(1:4),'red',mayx(1:3), mayy(1:3),'blue',junx(1:4),juny(1:4),'green',augx(1:5),augy(1:5),'black');

legend('April','May','June','August','April','May','June','August');

grid ;

radiation(1:N0,1:8) = 0;
for i = 1:N0,
    radiation(i,1) = radiation0(i,1);
    radiation(i,2) = radiation0(i,6);
    
    radiation(i,3) = radiation1(i,1);
    radiation(i,4) = radiation1(i,6);
    
    radiation(i,5) = radiation2(i,1);
    radiation(i,6) = radiation2(i,6);
    
    radiation(i,7) = radiation3(i,1);
    radiation(i,8) = radiation3(i,6);
end;

dlmwrite('radiation.dat',radiation,'delimiter',' ');
%writematrix(radiation,'radiation.dat','delimiter',' ');