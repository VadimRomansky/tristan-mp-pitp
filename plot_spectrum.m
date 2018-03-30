clear;
fp = hdf5read('./output/spect.001','specp');
fe = hdf5read('./output/spect.001','spece');
g=hdf5read('./output/spect.001','gamma');

Nx = size(fp,1);
Np = size(fp,2);

Fp(1:Np)=0;
Fe(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;

mp = 1.67262177*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;
c = 2.99792458*10^10;

normp = 0;
norme = 0;
norm = 1;

for i = 1:Np,
    Pp(i) = sqrt((g(i)+1)^2 - 1)*mp*c;
    Pe(i) = sqrt((g(i)+1)^2 - 1)*me*c;
    for j = 1:Nx,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
end;

normp = (Fp(1)/(Pp(2)^2))*(Pp(2) - Pp(1));
norme = (Fe(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normp = normp + (Fp(i)/(Pp(i)^2))*(Pp(i) - Pp(i-1));
    norme = norme + (Fe(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;

for i = 1:Np,
    Fp(i) = Fp(i)*norm/normp;
    Fe(i) = Fe(i)*norm/norme;
end;

figure(1);
plot (Pp(1:Np)/(mp*c),Fp(1:Np), 'red');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
grid ;

figure(2);
plot (Pe(1:Np)/(me*c),Fe(1:Np), 'red');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
grid ;