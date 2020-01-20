clear;
directory_name = './output/';
file_name = 'spect';
file_number = '.000';
Nd = 10;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'1','2','3','4','5','6','7','8','9','10'};
FileNumbers = {'.002','.004','.006','.008','.010','.012','.014','.016','.018','.020'};

full_name = strcat(directory_name, file_name, file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
Nx = fix(size(fp,1)/4);
%Nx = 12500;
Nx = 50000;

gmax(1:Nd) = 0;

for j = 1:Nd,
    full_name = strcat(directory_name, file_name, FileNumbers{j});
    fe = hdf5read(full_name,'spece');
    gam=hdf5read(full_name,'gamma');
    for i = 1:Np,
        for k = 1:Nx,
            if(fe(k,i) > 0)
                gmax(j) = gam(i);
            end;
        end;
    end;
end;


set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('{\gamma}');
xlabel ('N');
ylabel ('{\gamma}');

plot(1:Nd,gmax(1:Nd),'Color','red');

grid ;