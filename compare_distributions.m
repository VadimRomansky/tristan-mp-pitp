clear;
fp0 = hdf5read('./output/spect0.022','specp');
fe0 = hdf5read('./output/spect0.022','spece');
g0=hdf5read('./output/spect0.022','gamma');
fp1 = hdf5read('./output/spect1.022','specp');
fe1 = hdf5read('./output/spect1.022','spece');
g1=hdf5read('./output/spect1.022','gamma');
fp2 = hdf5read('./output/spect2.022','specp');
fe2 = hdf5read('./output/spect2.022','spece');
g2=hdf5read('./output/spect2.022','gamma');
fp3 = hdf5read('./output/spect3.022','specp');
fe3 = hdf5read('./output/spect3.022','spece');
g3=hdf5read('./output/spect3.022','gamma');
fp4 = hdf5read('./output/spect4.022','specp');
fe4 = hdf5read('./output/spect4.022','spece');
g4=hdf5read('./output/spect4.022','gamma');
fp5 = hdf5read('./output/spect5.022','specp');
fe5 = hdf5read('./output/spect5.022','spece');
g5=hdf5read('./output/spect5.022','gamma');
fp6 = hdf5read('./output/spect6.022','specp');
fe6 = hdf5read('./output/spect6.022','spece');
g6=hdf5read('./output/spect6.022','gamma');
fp7 = hdf5read('./output/spect7.022','specp');
fe7 = hdf5read('./output/spect7.022','spece');
g7=hdf5read('./output/spect7.022','gamma');
fp8 = hdf5read('./output/spect8.022','specp');
fe8 = hdf5read('./output/spect8.022','spece');
g8=hdf5read('./output/spect8.022','gamma');
fp9 = hdf5read('./output/spect9.022','specp');
fe9 = hdf5read('./output/spect9.022','spece');
g9=hdf5read('./output/spect9.022','gamma');

Nx = size(fp0,1);
Np = size(fp0,2);

Fp(1:Np,1:10)=0;
Fe(1:Np,1:10)=0;
g(1:Np,1:10) = 0;

Color = {[.7,.3,.3],'red','green','blue','black','yellow',[.5,.5,.5],'cyan',[.3,.7,.3], 'magenta',[1.0,.5,0],[.75,0.0,.7]};

for i = 1:Np,
    g(i,1)=g0(i);
    g(i,2)=g1(i);
    g(i,3)=g2(i);
    g(i,4)=g3(i);
    g(i,5)=g4(i);
    g(i,6)=g5(i);
    g(i,7)=g6(i);
    g(i,8)=g7(i);
    g(i,9)=g8(i);
    g(i,10)=g9(i);
    for j = 1:Nx,
        Fp(i,1) = Fp(i,1) + fp0(j,i);
        Fe(i,1) = Fe(i,1) + fe0(j,i);
        Fp(i,2) = Fp(i,2) + fp1(j,i);
        Fe(i,2) = Fe(i,2) + fe1(j,i);
        Fp(i,3) = Fp(i,3) + fp2(j,i);
        Fe(i,3) = Fe(i,3) + fe2(j,i);
        Fp(i,4) = Fp(i,4) + fp3(j,i);
        Fe(i,4) = Fe(i,4) + fe3(j,i);
        Fp(i,5) = Fp(i,5) + fp4(j,i);
        Fe(i,5) = Fe(i,5) + fe4(j,i);
        Fp(i,6) = Fp(i,6) + fp5(j,i);
        Fe(i,6) = Fe(i,6) + fe5(j,i);
        Fp(i,7) = Fp(i,7) + fp6(j,i);
        Fe(i,7) = Fe(i,7) + fe6(j,i);
        Fp(i,8) = Fp(i,8) + fp7(j,i);
        Fe(i,8) = Fe(i,8) + fe7(j,i);
        Fp(i,9) = Fp(i,9) + fp8(j,i);
        Fe(i,9) = Fe(i,9) + fe8(j,i);
        Fp(i,10) = Fp(i,10) + fp9(j,i);
        Fe(i,10) = Fe(i,10) + fe9(j,i);
    end;
    for j=1:10,
        Fp(i,j)=Fp(i,j)*(g(i,j) + 1)^2;
        Fe(i,j)=Fe(i,j)*(g(i,j) + 1)^2;
    end;
end;

set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',20,'DefaultTextFontName','Times New Roman'); 
figure(1);
hold on;
for j=1:10,
    plot (1+g(1:Np,j),Fp(1:Np,j),'color',Color{j});
end;
legend('{\theta} = 0^{\circ}','{\theta} = 10^{\circ}','{\theta} = 20^{\circ}','{\theta} = 30^{\circ}','{\theta} = 40^{\circ}','{\theta} = 50^{\circ}','{\theta} = 60^{\circ}','{\theta} = 70^{\circ}', '{\theta} = 80^{\circ}','{\theta} = 90^{\circ}','Location','southeast');
grid ;

figure(2);
hold on;
for j=1:10,
    plot (1+g(1:Np,j),Fe(1:Np,j),'color',Color{j});
end;
legend('{\theta} = 0^{\circ}','{\theta} = 10^{\circ}','{\theta} = 20^{\circ}','{\theta} = 30^{\circ}','{\theta} = 40^{\circ}','{\theta} = 50^{\circ}','{\theta} = 70^{\circ}','{\theta} = 90^{\circ}', 'maxwell','maxwell-juttner','Location','southeast');
grid ;