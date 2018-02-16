clear;
fp = hdf5read('./output/spect.023','specp');
fe = hdf5read('./output/spect.023','spece');
g=hdf5read('./output/spect.023','gamma');

Nx = size(fp,1);
Np = size(fp,2);

Fp(1:Np)=0;
Fe(1:Np)=0;

for i = 1:Np,
    for j = 1:Nx,
        Fp(i) = Fp(i) + fp(j,i);
        Fe(i) = Fe(i) + fe(j,i);
    end;
    Fp(i)=Fp(i)*(g(i) + 1)^2;
    Fe(i)=Fe(i)*(g(i) + 1)^2;
end;

figure(1);
plot (1+g(1:Np),Fp(1:Np), 'red');
title ('Fp');
xlabel ('gamma');
ylabel ('Fp*p^4');
grid ;

figure(2);
plot (1+g(1:Np),Fe(1:Np), 'red');
title ('Fe');
xlabel ('gamma');
ylabel ('Fe*p^4');
grid ;