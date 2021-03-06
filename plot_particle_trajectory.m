clear;
directory_name = './output4/';
file_name = 'flds.tot';
part_name = 'prtl.tot';
file_number = '.005';
full_name = strcat(directory_name, file_name, file_number);
full_part_name = strcat(directory_name, part_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
Ex = hdf5read(full_name,'ex');
Ey = hdf5read(full_name,'ey');
Ez = hdf5read(full_name,'ez');
fileinfo = hdf5info(full_part_name);
last_number = 50;
a = last_number;
first_number = 1;

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
inde = hdf5read(full_part_name, 'inde');
proce = hdf5read(full_part_name, 'proce');
xe = hdf5read(full_part_name, 'xe');
ye = hdf5read(full_part_name, 'ye');
ze = hdf5read(full_part_name, 'ze');
ve = hdf5read(full_part_name, 've');
ue = hdf5read(full_part_name, 'ue');
we = hdf5read(full_part_name, 'we');
gammai = hdf5read(full_part_name, 'gammai');
indi = hdf5read(full_part_name, 'indi');
xi = hdf5read(full_part_name, 'xi');
yi = hdf5read(full_part_name, 'yi');
zi = hdf5read(full_part_name, 'zi');
%Ex = hdf5read(full_name,'ex');
%Ey = hdf5read(full_name,'ey');
%Ez = hdf5read(full_name,'ez');

Nskinlength = 5;

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
samplingFactor = 20;
fieldFactor = me*rho/(q*tau*tau);
rho = rho*samplingFactor;
rho = 20;

Nx = size(Bx, 1);
Ny = size(By, 2);
endx = Nx/10;

B0 = 0;
for i = 1:endx,
    for j = 1:Ny,
        B0 = B0 + Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j);
    end;
end;
B0 = sqrt(B0/(endx*Ny));

frameTime = 1.0/10;

maxGamma1 = 1.0;
max_number1 = 1;
maxGamma2 = 1.0;
max_number2 = 1;
maxGamma3 = 1.0;
max_number3 = 1;
for i = 1: size(gammae,1),
    if(xe(i) < (endx/2)*samplingFactor)
    if (gammae(i) > maxGamma1)        
        max_number3 = max_number2;
        maxGamma3 = maxGamma2;
        max_number2 = max_number1;
        maxGamma2 = maxGamma1;
        max_number1 = i;
        maxGamma1 = gammae(i);
    else
        if (gammae(i) > maxGamma2) 
            max_number3 = max_number2;
            maxGamma3 = maxGamma2;
            max_number2 = i;
            maxGamma2 = gammae(i);
        else
            if (gammae(i) > maxGamma3)
                max_number3 = i;
                maxGamma3 = gammae(i);
            end
        end
    end
    end
end;

Npart = size(gammae,1);
part_index = inde(max_number1);
part_number = max_number1;

%part_number = 400000;
part_indes = inde(part_number);

part_proc = proce(part_number);

part_number2 = max_number2;
part_number3 = max_number3;


%part_number2 = 6000;
part_index2 = inde(part_number2);
part_proc2 = proce(part_number2);

%part_number3 = 9000;
part_index3 = inde(part_number3);
part_proc3 = proce(part_number3);

Bnorm(1:Ny, 1:endx) = 0;
Enorm(1:Ny, 1:endx) = 0;



for i=1:endx,
    for j = 1:Ny,
        Bnorm(j,i) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        Bnorm(j,i) = sqrt(Ex(i,j)*Ex(i,j) + Ey(i,j)*Ey(i,j) + Ez(i,j)*Ez(i,j));
       % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
    end;
end;

%for i=1:Nx/2,
%    for j = 1:Nx/2,
%        k = rem(j,Ny) + 1;
%        Bnorm(i,j) = Bnorm(i,k);
%    end;
%end;



x1(1:(last_number-first_number + 1)) = 0;
y1(1:(last_number-first_number + 1)) = 0;
g1(1:(last_number-first_number + 1)) = 0;
x2(1:(last_number-first_number + 1)) = 0;
y2(1:(last_number-first_number + 1)) = 0;
g2(1:(last_number-first_number + 1)) = 0;
x3(1:(last_number-first_number + 1)) = 0;
y3(1:(last_number-first_number + 1)) = 0;
g3(1:(last_number-first_number + 1)) = 0;

%figure(1)
figure('Position', [10 10 900 600]);
xlabel ('x');
ylabel ('\gamma');
grid on;
hold on;
for a = first_number:last_number,
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
    Npart = size(gammae,1);
    xe = hdf5read(full_part_name, 'xe');
    ye = hdf5read(full_part_name, 'ye');
    ze = hdf5read(full_part_name, 'ze');
    inde = hdf5read(full_part_name, 'inde');
    proce = hdf5read(full_part_name, 'proce');
    
    for p = 1:Npart,
        if((inde(p) == part_index) && (proce(p)==part_proc))
            part_number = p;
        end;
        if((inde(p) == part_index2) && (proce(p)==part_proc2))
            part_number2 = p;
        end;
        if((inde(p) == part_index3) && (proce(p)==part_proc3))
            part_number3 = p;
        end;
    end;
    if(a == first_number)
        plot(xe(part_number), gammae(part_number), 'ro', 'MarkerSize', 15);
        plot(xe(part_number2), gammae(part_number2), 'ro', 'MarkerSize', 15, 'Color','green');
        plot(xe(part_number3), gammae(part_number3), 'ro', 'MarkerSize', 15, 'Color', 'black');
    end;
    x1(a - first_number+1) = xe(part_number);
    y1(a - first_number+1) = ye(part_number);
    g1(a - first_number+1) = gammae(part_number);
    x2(a - first_number+1) = xe(part_number2);
    y2(a - first_number+1) = ye(part_number2);
    g2(a - first_number+1) = gammae(part_number2);
    x3(a - first_number+1) = xe(part_number3);
    y3(a - first_number+1) = ye(part_number3);
    g3(a - first_number+1) = gammae(part_number3);
end;
plot(x1(1:(last_number-first_number + 1)),g1(1:(last_number-first_number + 1)),'red',x2(1:(last_number-first_number + 1)),g2(1:(last_number-first_number + 1)),'green', x3(1:(last_number-first_number + 1)),g3(1:(last_number-first_number + 1)),'black');
grid;
timeStep = 1;
%figure(2);
figure('Position', [10 50 1200 600]);
%title ('E_x');
xlabel ('Nt');
ylabel ('\gamma');
grid on;
hold on;
plot((1:(last_number-first_number + 1))*timeStep,g1(1:(last_number-first_number + 1)),'red');
plot((1:(last_number-first_number + 1))*timeStep,g2(1:(last_number-first_number + 1)),'green');
plot((1:(last_number-first_number + 1))*timeStep,g3(1:(last_number-first_number + 1)),'black');
fig_part = plot(1, g1(1), 'ro', 'MarkerSize', 10,'Color','red');
fig_part2 = plot(1, g1(1), 'ro', 'MarkerSize', 10,'Color','green');
fig_part3 = plot(1, g1(1), 'ro', 'MarkerSize', 10,'Color','black');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1)=0;
f = getframe(gcf);
[mov(:,:,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = first_number:last_number,
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
   
    set(fig_part, 'Xdata', (a-first_number+1)*timeStep);
    set(fig_part, 'Ydata', g1(a-first_number+1));
    set(fig_part2, 'Xdata', (a-first_number+1)*timeStep);
    set(fig_part2, 'Ydata', g2(a-first_number+1));
    set(fig_part3, 'Xdata', (a-first_number+1)*timeStep);
    set(fig_part3, 'Ydata', g3(a-first_number+1));
    f = getframe(gcf);
    %mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    outname = strcat(directory_name,'g',int2str(a),'.jpg');
    imwrite(rgb2ind(f.cdata, map), map, outname);
end;
imwrite(mov, map, 'T.gif','DelayTime',frameTime*3,'LoopCount',1);


%figure(3);
figure('Position', [10 50 1200 600]);
%title ('E_x');
xlabel ('Nx');
ylabel ('Ny');
grid on;
hold on;
%axis([Xgrid(1) Xgrid(Nx-1) minEx maxEx]);
%fig = plot (Xgrid(1:Nx-1),Ex(1:Nx-1), 'red');
caxis ([0 3])
fig = imagesc((1:endx)*samplingFactor, (1:Ny)*samplingFactor,Bnorm);
fig_part = plot(xe(part_number), ye(part_number), 'ro', 'MarkerSize', 10, 'Color','red','LineWidth',3);
fig_part2 = plot(xe(part_number2), ye(part_number2), 'ro', 'MarkerSize', 10, 'Color','green','LineWidth',3);
fig_part3 = plot(xe(part_number3), ye(part_number3), 'ro', 'MarkerSize', 10, 'Color', 'black','LineWidth',3);

plot(x1(1:(last_number-first_number + 1)),y1(1:(last_number-first_number + 1)),'red');
plot(x2(1:(last_number-first_number + 1)),y2(1:(last_number-first_number + 1)),'green');
plot(x3(1:(last_number-first_number + 1)),y3(1:(last_number-first_number + 1)),'black');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1)=0;
f = getframe(gcf);
[mov(:,:,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = first_number:last_number,
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
    
    Bx = hdf5read(full_name,'bx');
    By = hdf5read(full_name,'by');
    Bz = hdf5read(full_name,'bz');
    
    
    for i=1:endx,
        for j = 1:Ny,
            Bnorm(j,i) = sqrt(Bx(i,j)*Bx(i,j) + By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j))/B0;
        % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        end;
    end;

    
    set(fig, 'CData', Bnorm);
    
    %for p = 1:Npart,
    %    if((inde(p) == part_index) && (proce(p)==part_proc))
    %        part_number = p;
    %    end;
    %    if((inde(p) == part_index2) && (proce(p)==part_proc2))
    %        part_number2 = p;
    %    end;
    %    if((inde(p) == part_index3) && (proce(p)==part_proc3))
    %        part_number3 = p;
    %    end;
    %end;
    
    set(fig_part, 'Xdata', x1(a-first_number+1));
    set(fig_part, 'Ydata', y1(a-first_number+1));
    set(fig_part2, 'Xdata', x2(a-first_number+1));
    set(fig_part2, 'Ydata', y2(a-first_number+1));
    set(fig_part3, 'Xdata', x3(a-first_number+1));
    set(fig_part3, 'Ydata', y3(a-first_number+1));
    f = getframe(gcf);
    %mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    outname = strcat(directory_name,'T',int2str(a),'.jpg');
    imwrite(rgb2ind(f.cdata, map), map, outname);
    pause(1);
end;
%imwrite(mov, map, 'T.gif','DelayTime',frameTime*3,'LoopCount',1);

figure('Position', [10 50 1200 600]);
%title ('E_x');
xlabel ('Nx');
ylabel ('Ny');
grid on;
hold on;
%axis([Xgrid(1) Xgrid(Nx-1) minEx maxEx]);
%fig = plot (Xgrid(1:Nx-1),Ex(1:Nx-1), 'red');
caxis ([0 3])
fig = imagesc((1:endx)*samplingFactor, (1:Ny)*samplingFactor,Enorm);
fig_part = plot(xe(part_number), ye(part_number), 'ro', 'MarkerSize', 10, 'Color','red','LineWidth',3);
fig_part2 = plot(xe(part_number2), ye(part_number2), 'ro', 'MarkerSize', 10, 'Color','green','LineWidth',3);
fig_part3 = plot(xe(part_number3), ye(part_number3), 'ro', 'MarkerSize', 10, 'Color', 'black','LineWidth',3);
plot(x1(1:(last_number-first_number + 1)),y1(1:(last_number-first_number + 1)),'red');
plot(x2(1:(last_number-first_number + 1)),y2(1:(last_number-first_number + 1)),'green');
plot(x3(1:(last_number-first_number + 1)),y3(1:(last_number-first_number + 1)),'black');
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
mov(1:height, 1:width, 1:1)=0;
f = getframe(gcf);
[mov(:,:,1), map]=rgb2ind(f.cdata, colorcube(256));
for a = first_number:last_number,
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
    
    Ex = hdf5read(full_name,'ex');
    Ey = hdf5read(full_name,'ey');
    Ez = hdf5read(full_name,'ez');
    
    for i=1:endx,
        for j = 1:Ny,
            Enorm(j,i) = sqrt(Ex(i,j)*Ex(i,j) + Ey(i,j)*Ey(i,j) + Ez(i,j)*Ez(i,j))/B0;
        % Bperp(i,j) = sqrt(By(i,j)*By(i,j) + Bz(i,j)*Bz(i,j));
        end;
    end;

    
    set(fig, 'CData', Enorm);
    
    %for p = 1:Npart,
    %    if((inde(p) == part_index) && (proce(p)==part_proc))
    %        part_number = p;
    %    end;
    %    if((inde(p) == part_index2) && (proce(p)==part_proc2))
    %        part_number2 = p;
    %    end;
    %    if((inde(p) == part_index3) && (proce(p)==part_proc3))
    %        part_number3 = p;
    %    end;
    %end;
    
    set(fig_part, 'Xdata', x1(a-first_number + 1));
    set(fig_part, 'Ydata', y1(a-first_number + 1));
    set(fig_part2, 'Xdata', x2(a-first_number + 1));
    set(fig_part2, 'Ydata', y2(a-first_number + 1));
    set(fig_part3, 'Xdata', x3(a-first_number + 1));
    set(fig_part3, 'Ydata', y3(a-first_number + 1));
    f = getframe(gcf);
    %mov(:,:,1,a+1)=rgb2ind(f.cdata, map);
    outname = strcat(directory_name,'TE',int2str(a),'.jpg');
    imwrite(rgb2ind(f.cdata, map), map, outname);
end;
%imwrite(mov, map, 'T.gif','DelayTime',frameTime*3,'LoopCount',1);

