clear;
directory_name = './output/';
file_name = 'spect_ang';
file_number = '.010';
full_name = strcat(directory_name, file_name, file_number);
fileinfo = hdf5info(full_name);
fp = hdf5read(full_name,'specp');
fe = hdf5read(full_name,'spece');
g=hdf5read(full_name,'gamma');

Nx = size(fp,1);
Np = size(fp,2);
Nmu = size(fp,3);
Nphi = size(fp,4);

startx = 1;
endx = Nx/4;

Fp(1:Np,1:Nmu)=0;
Fe(1:Np,1:Nmu)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;
Fejuttner(1:Np)=0;
Fpjuttner(1:Np)=0;

me = 0.91*10^-27;
mass_ratio = 16;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 0.5;
fractionp = 0.5;

normp = 0;
norme = 0;
norm = 1;

for i = 1:Np,
    %Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    %Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    Pp(i) = sqrt((g(i)+1)^2 - 1);
    Pe(i) = sqrt((g(i)+1)^2 - 1);
    for k = 1:Nmu,
        for j = startx:endx,
            for l = 1:Nphi
                a = fp(j,i,k,l);
                if(a > 0)
                    b = 0;
                end
                Fp(i,k) = Fp(i,k) + fp(j,i,k,l);
                Fe(i,k) = Fe(i,k) + fe(j,i,k,l);
            end;
        end;
        Fp(i,k)=Fp(i,k)*(Pp(i))/(1+g(i));
        Fe(i,k)=Fe(i,k)*(Pe(i))/(1+g(i));
    end;
end;

normp = (Fp(1)/(Pp(2)^2))*(Pp(2) - Pp(1));
norme = (Fe(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

%for i = 2:Np,
%    normp = normp + (Fp(i)/(Pp(i)^2))*(Pp(i) - Pp(i-1));
%    norme = norme + (Fe(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
%end;

%for i = 1:Np,
%    Fp(i) = Fp(i)*norm/normp;
%    Fe(i) = Fe(i)*norm/norme;
%end;

figure(1);
[X, Y] = meshgrid(2*(1:Nmu)/Nmu - 1 , Pp(1:Np));
surf(X, Y, Fp);
shading interp;
title ('Fp');
xlabel ('mu');
ylabel ('p/mc');
zlabel ('Fp');
grid ;

figure(2);
[X, Y] = meshgrid(2*(1:Nmu)/Nmu - 1 , Pe(1:Np));
surf(X, Y, Fe);
shading interp;
title ('Fe');
xlabel ('mu');
ylabel ('p/mc');
zlabel ('Fe');
grid ;