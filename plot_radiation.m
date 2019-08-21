clear;
radiation = importdata('radiation3.dat');
N = size(radiation,1);

index(1:N) = 0;
for i = 2:N,
    index(i) = (log(radiation(i,2)) - log(radiation(i-1,2)))/(log(radiation(i,1)) - log(radiation(i-1,1)));
end;
index(1) = index(2);

augx(1:4) = 0;
augx(1) = 0.57;
augx(2) = 1.26;
augx(3) = 3.84;
augx(4) = 6.32;
augy(1:4) = 0;
augy(1) = 7.9;
augy(2) =8.7;
augy(3) = 2.45;
augy(4) = 1.07;


factor = 2.5*10^-27;

figure(1);
hold on;
plot(radiation(1:N,1),radiation(1:N,2),'red');
plot(augx(1:4),augy(1:4)*factor,'blue');
title ('I_{\nu}');
xlabel ('{\nu}/{\nu}_c');
ylabel ('I_{\nu} erg/cm^{-2} sr');
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