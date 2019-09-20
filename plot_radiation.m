clear;
radiation = importdata('radiation3.dat');
N = size(radiation,1);

index(1:N) = 0;
for i = 2:N,
    index(i) = (log(radiation(i,2)) - log(radiation(i-1,2)))/(log(radiation(i,1)) - log(radiation(i-1,1)));
end;
index(1) = index(2);

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


factor = 1.0;

figure(1);
hold on;
plot(radiation(1:N,1),radiation(1:N,2),'red');
%plot(radiation(1:N,1),radiation(1:N,7),'green');
plot(augx(1:5),augy(1:5)*factor,'blue');
title ('I_{\nu}');
xlabel ('{\nu} GHz');
ylabel ('I_{\nu} erg/Hz cm^{3} sr');
grid ;

figure(2);
plot(radiation(1:N,1),index(1:N),'red');
title ('{\gamma}');
xlabel ('{\nu}/{\nu}_c');
ylabel ('{\gamma}');
grid ;

figure(3);
plot(radiation(1:N,1),radiation(1:N,3),'red');
title ('t');
xlabel ('{\nu}/{\nu}_c');
ylabel ('t');
grid ;

radiation_wrt(1:N,2) = 0;
for i = 1:N
    radiation_wrt(i,1) = radiation(i,1);
    radiation_wrt(i,2) = radiation(i,2);
end;

dlmwrite('radiation.dat',radiation_wrt,'delimiter',' ');