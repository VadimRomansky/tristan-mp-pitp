clear;
directory_name = './output/';
file_name = 'flds.tot';
part_name = 'prtl.tot';
file_number = '.001';
full_name = strcat(directory_name, file_name, file_number);
full_part_name = strcat(directory_name, part_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
fileinfo = hdf5info(full_part_name);
last_number = 10;
a = last_number;

if(a < 10)
        full_name = strcat(directory_name, file_name, '.00', num2str(a));
        full_part_name = strcat(directory_name, part_name, '.00', num2str(a));
    else if (a < 100)
            full_name = strcat(directory_name, file_name, '.0', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.0', num2str(a));  
        else 
            full_name = strcat(directory_name, file_name, '.', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.', num2str(a));
        end;
    end;

gammae = hdf5read(full_part_name, 'gammae');
xe = hdf5read(full_part_name, 'xe');
ye = hdf5read(full_part_name, 'ye');
ze = hdf5read(full_part_name, 'ze');
ve = hdf5read(full_part_name, 've');
ue = hdf5read(full_part_name, 'ue');
we = hdf5read(full_part_name, 'we');
gammai = hdf5read(full_part_name, 'gammai');
xi = hdf5read(full_part_name, 'xi');
yi = hdf5read(full_part_name, 'yi');
zi = hdf5read(full_part_name, 'zi');
%Ex = hdf5read(full_name,'ex');
%Ey = hdf5read(full_name,'ey');
%Ez = hdf5read(full_name,'ez');
B0 = 0.03030750;
Nx = size(Bx, 1);
Ny = size(By, 2);

frameTime = 1.0/last_number;

maxGamma = 1.0;
max_number = 1;
for i = 1: size(gammae,1),
    if (gammae(i) > maxGamma)
        max_number = i;
        maxGamma = gammae(i);
    end;
end;

part_number = max_number;

Bnorm(1:Ny, 1:Nx) = 0;




for i=1:Nx,
    for j = 1:Ny,
        Bnorm(j,i) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

%for i=1:Nx/2,
%    for j = 1:Nx/2,
%        k = rem(j,Ny) + 1;
%        Bnorm(i,j) = Bnorm(i,k);
%    end;
%end;

Nskinlength = 10;

c0 = 2.998*10^10;
mass_ratio = 20;
mp = 1.67262*10^-24;
me = mp/mass_ratio;
q = 4.80320427*10^-10;
n = 10^-4;

omega = sqrt(4*pi*n*q*q/me);

rho = c0/(omega*Nskinlength);
rho = 0.1;
c1=0.45;

tau = c1*rho/c0;
samplingFactor = 1;
fieldFactor = me*rho/(q*tau*tau);
rho = rho*samplingFactor;


%figure(7);
%hold on;
%colormap Jet;
%[X, Y] = meshgrid((1:Ny)*rho, (1:Nx/2)*rho);
%surf(X, Y, Bnorm/B0);
%contourf(X,Y,Bnorm);
%imagesc(1:Ny,1:Nx/2,Bnorm);
%plot(1:Ny,(1:Ny)*10,'red');
%plot(Ny/2, Nx/4, 'ro', 'MarkerSize', 10);
%shading interp;
%title ('B/B_0');
%xlabel ('y \omega /c');
%ylabel ('x \omega /c');
%zlabel ('B/B_0');
%grid ;

figure(1);
%title ('E_x');
xlabel ('x\omega_p/c');
ylabel ('E_x gauss');
grid on;
hold on;
%axis([Xgrid(1) Xgrid(Nx-1) minEx maxEx]);
%fig = plot (Xgrid(1:Nx-1),Ex(1:Nx-1), 'red');
fig = imagesc((1:Nx)*samplingFactor, (1:Ny)*samplingFactor,Bnorm);
fig_part = plot(xe(part_number), ye(part_number), 'ro', 'MarkerSize', 30);
%pos = get(gcf, 'Position');
%width = pos(3);
%height = pos(4);
%mov(1:height, 1:width, 1:1, 1:last_number)=0;
%f = getframe(gcf);
%[mov(:,:,1,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = 1:last_number,
    if(a < 10)
        full_name = strcat(directory_name, file_name, '.00', num2str(a));
        full_part_name = strcat(directory_name, part_name, '.00', num2str(a));
    else if (a < 100)
            full_name = strcat(directory_name, file_name, '.0', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.0', num2str(a));  
        else 
            full_name = strcat(directory_name, file_name, '.', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.', num2str(a));
        end;
    end;
    for i=1:Nx,
        for j = 1:Ny,
            Bnorm(j,i) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        end;
    end;
    Bx = hdf5read(full_name,'bx');
    By = hdf5read(full_name,'by');
    Bz = hdf5read(full_name,'bz');

    gammae = hdf5read(full_part_name, 'gammae');
    xe = hdf5read(full_part_name, 'xe');
    ye = hdf5read(full_part_name, 'ye');
    ze = hdf5read(full_part_name, 'ze');
    
    set(fig, 'CData', Bnorm);
    
    set(fig_part, 'Xdata', xe(part_number));
    set(fig_part, 'Ydata', ye(part_number));
    %f = getframe(gcf);
    %mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    pause(frameTime*10);
end;
%imwrite(mov, map, 'Ex.gif','DelayTime',frameTime,'LoopCount',1);

