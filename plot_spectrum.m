clear;
fp = hdf5read('./output/spect0.001','specp');
fe = hdf5read('./output/spect0.001','spece');
g=hdf5read('./output/spect0.001','gamma');

Nx = size(fp,1);
Np = size(fp,2);

Fp(1:Np)=0;
Fe(1:Np)=0;

Pp(1:Np)=0;
Pe(1:Np)=0;

for i = 1:Np,
    Pp(i) = sqrt((g(i)+1)^2 - 1);
    Pe(i) = sqrt((g(i)+1)^2 - 1);
    for j = 1:Nx,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
    Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
end;

figure(1);
plot (Pp(1:Np),Fp(1:Np), 'red');
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
grid ;

figure(2);
plot (Pe(1:Np),Fe(1:Np), 'red');
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
grid ;