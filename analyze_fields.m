clear;
directory_name = './output1/';
file_name = 'flds.tot';
file_number = '.000';
full_name = strcat(directory_name, file_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');
B0 = 0.03030750;
Nx = size(Bx, 1);
Ny = size(By, 2);

Bnorm(1:Nx, 1:Ny) = 0;
Ediff(1:Nx, 1:Ny) = 0;

fourierBx(1:Ny, 1:Ny) = 0;
fourierBy(1:Ny, 1:Ny) = 0;
fourierBz(1:Ny, 1:Ny) = 0;

fourierBx2(1:Ny, 1:Ny) = 0;
fourierBy2(1:Ny, 1:Ny) = 0;
fourierBz2(1:Ny, 1:Ny) = 0;

Ndistr = 100;
koef = 10;
Bdistribution(1:Ndistr) = 0;
%Bperp(1:Nx, 1:Ny) = 0;
g = 1.5;
beta = sqrt(1-1/(g*g));

for i=1:Nx,
    for j = 1:Ny,
        Bnorm(i,j) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        Ediff(i,j) = sqrt((Ex(i,j))*(Ex(i,j)) + (Ey(i,j) +beta*Bz(i,j))*(Ey(i,j) +beta*Bz(i,j)) + (Ez(i,j) - beta*By(i,j))*(Ez(i,j) - beta*By(i,j)))/Bnorm(i,j);
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

for i=1:Ny,
    for j = 1:Ny,
        fourierBx(i,j) = Bx(i+100,j);
        %fourierBx(i,j) = sin((2*pi/Ny)*(i + j)) + sin((2*pi/Ny)*(i)) + sin((2*pi/Ny)*j) + sin((2*pi/Ny)*(2*i + 2*j));
        fourierBy(i,j) = By(i+100,j);
        fourierBz(i,j) = Bz(i+100,j);
    end;
end;

fourierBx2 = fft2(fourierBx);
fourierBy2 = fft2(fourierBy);
fourierBz2 = fft2(fourierBz);

miniFourierBx(1:10,1:10) = 0;
miniFourierBy(1:10,1:10) = 0;
miniFourierBz(1:10,1:10) = 0;
for i=1:10,
    for j = 1:10,
        miniFourierBx(i,j) = abs(fourierBx2(i,j));
        miniFourierBy(i,j) = abs(fourierBy2(i,j));
        miniFourierBz(i,j) = abs(fourierBz2(i,j));
    end;
end;

for i = 1:Ny,
    for j = 1:Ny,
        k = fix(koef*Bnorm(i,j)/B0)+1;
        if(k <= Ndistr)
            Bdistribution(k) = Bdistribution(k) + 1;
        end;
    end;
end;

figure(1);
hold on;
title ('distribution');
xlabel ('B/B0');
ylabel ('N');
plot ((1:Ndistr)/koef,Bdistribution(1:Ndistr));
grid ;



