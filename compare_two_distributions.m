clear;
fp0 = hdf5read('./output/spect1.001','specp');
fe0 = hdf5read('./output/spect1.001','spece');
g0=hdf5read('./output/spect1.001','gamma');
fp1 = hdf5read('./output/spect3.001','specp');
fe1 = hdf5read('./output/spect3.001','spece');
g1=hdf5read('./output/spect3.001','gamma');

Nx = size(fp0,1);
Np = size(fp0,2);

Fp(1:Np,1:2)=0;
Fe(1:Np,1:2)=0;
Pp(1:Np,1:2)=0;
Pe(1:Np,1:2)=0;
g(1:Np,1:2) = 0;

Color = {[.7,.3,.3],'red','green','blue','black','yellow',[.5,.5,.5],'cyan',[.3,.7,.3], 'magenta',[1.0,.5,0],[.75,0.0,.7]};

for i = 1:Np,
    g(i,1)=g0(i);
    g(i,2)=g1(i);
    for k = 1:2,
        Pp(i,k) = sqrt((g(i,k)+1)^2 - 1);
        Pe(i,k) = sqrt((g(i,k)+1)^2 - 1);
    end;
    for j = 1:Nx,
        Fp(i,1) = Fp(i,1) + fp0(j,i);
        Fe(i,1) = Fe(i,1) + fe0(j,i);
        Fp(i,2) = Fp(i,2) + fp1(j,i);
        Fe(i,2) = Fe(i,2) + fe1(j,i);
    end;
    for j=1:2,
        Fp(i,j)=Fp(i,j)*(Pp(i,j)^3)/(1+g(i,j));
        Fe(i,j)=Fe(i,j)*(Pe(i,j)^3)/(1+g(i,j));;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
title ('F_p');
xlabel ('p/{m_p c}');
ylabel ('Fp*p^4');
for j=1:2,
    plot (Pp(1:Np,j),Fp(1:Np,j),'color',Color{j});
end;
legend('1','2','Location','southeast');
grid ;

figure(2);
hold on;
title ('F_e');
xlabel ('p/{m_e c}');
ylabel ('F_e*p^4');
for j=1:2,
    plot (Pe(1:Np,j),Fe(1:Np,j),'color',Color{j});
end;
legend('1','2','Location','southeast');
grid ;