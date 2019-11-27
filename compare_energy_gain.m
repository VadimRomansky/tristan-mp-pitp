clear;
directory_name = './output4/';
file_name = 'spect';
file_number = '.000';
Nd = 10;
Ns = 4;
start = 0;

Color = {'red','blue','green','black','cyan','magenta','yellow',[0.75,0,0.67],[0.5,0.5,0.0],[.98,.5,.44]};
LegendTitle = {'1','2','3','4','5','6','7','8','9','10'};
%FileNumbers = {'.000','.003','.006','.009','.012','.015','.018','.021','.024','.027'};
FileNumbers = {'.000','.010','.020','.030','.040','.050','.060','.070','.080','.090'};
%FileNumbers = {'.000','.001','.002','.003','.004','.005','.006','.007','.008','.009','.010','.011','.012','.013','.014','.015','.016','.017','.018','.019','.020','.021','.022','.023','.024','.025','.026','.027','.028','.029','.030','.031','.032','.033','.034','.035','.036','.037','.038','.039','.040','.041','.042','.043','.044','.045','.046','.047','.048','.049','.050','.051','.052','.053','.054','.055','.056','.057','.058','.059','.060','.061','.062','.063','.064','.065','.066','.067','.068','.069','.070','.071','.072','.073','.074','.075','.076','.077','.078','.079','.080','.081','.082','.083','.084','.085','.086','.087','.088','.089','.090','.091','.092','.093','.094','.095','.096','.097','.098','.099','.100'};
Specifiers = {'-','--',':','-.'};

full_name = strcat(directory_name, file_name,num2str(start), file_number);
fp = hdf5read(full_name,'specp');
Np = size(fp,2);
%Nx = fix(size(fp,1)/4);
Nx = 10000;
%Nx = 12500;

startx = 1;
endx = Nx;

g(1:Np) = 0;
Fp(1:Np)=0;
Fe(1:Np)=0;
Pp(1:Np)=0;
Pe(1:Np)=0;

energyFraction(1:Ns,1:Nd) = 0;
energyFraction1(1:Ns,1:Nd) = 0;
energyFraction2(1:Ns,1:Nd) = 0;
energyFraction3(1:Ns,1:Nd) = 0;
energyFraction4(1:Ns,1:Nd) = 0;
totalEnergy = 0;
highEnergy = 0;
highEnergy1(1:Ns,1:Nd) = 0;
highEnergy2(1:Ns,1:Nd) = 0;
highEnergy3(1:Ns,1:Nd) = 0;
highEnergy4(1:Ns,1:Nd) = 0;
gammaLevel = 100;
gammaLevel1 = 1;
gammaLevel2 = 50;
gammaLevel3 = 200;
gammaLevel4 = 500;
width(1:Ns) = 0;
width(1) = 4000;
width(2) = 2000;
width(3) = 1000;
width(4) = 1000;

me = 0.91*10^-27;
mass_ratio = 25;
mp = me*mass_ratio;
c = 2.99792458*10^10;
Te = 9*10^9;
Tp = 3.5*10^10;
kB = 1.3806488*10^-16;
thetae = kB*Te/(me*c*c);
thetap = kB*Tp/(mp*c*c);
fractione = 1.0;
fractionp = 1.0;

for m = 1:Ns,
    for j = 1:Nd,
        full_name = strcat(directory_name, file_name,num2str(start + m - 1), FileNumbers{j});
        fp = hdf5read(full_name,'specp');
        fe = hdf5read(full_name,'spece');
        gam=hdf5read(full_name,'gamma');
        totalEnergy = 0;
        highEnergy = 0;
        for i = 1:Np,
            Fe(i) = 0;
            Fp(i) = 0;
            g(i) = gam(i);
            Pp(i) = sqrt((g(i)+1)^2 - 1);
            Pe(i) = sqrt((g(i)+1)^2 - 1);
            for k = startx:endx,
                Fp(i) = Fp(i) + fp(k,i);
                Fe(i) = Fe(i) + fe(k,i);
            end;
            totalEnergy = totalEnergy + Fe(i)*(Pe(i));
            if(g(i) > gammaLevel - 1)
                highEnergy = highEnergy + Fe(i)*(Pe(i));
            end;
            if(g(i) > gammaLevel1 - 1)
                highEnergy1(m,j) = highEnergy1(m,j) + Fe(i)*(Pe(i));
            end;
            if(g(i) > gammaLevel2 - 1)
                highEnergy2(m,j) = highEnergy2(m,j) + Fe(i)*(Pe(i));
            end;
            if(g(i) > gammaLevel3 - 1)
                highEnergy3(m,j) = highEnergy3(m,j) + Fe(i)*(Pe(i));
            end;
            if(g(i) > gammaLevel4 - 1)
                highEnergy4(m,j) = highEnergy4(m,j) + Fe(i)*(Pe(i));
            end;
            Fp(i)=Fp(i)*(Pp(i)^3)/(1+g(i));
            Fe(i)=Fe(i)*(Pe(i)^3)/(1+g(i));
        end;
        energyFraction(m,j) = highEnergy/totalEnergy;
        energyFraction1(m,j) = highEnergy1(m,j)/totalEnergy;
        energyFraction2(m,j) = highEnergy2(m,j)/totalEnergy;
        energyFraction3(m,j) = highEnergy3(m,j)/totalEnergy;
        energyFraction4(m,j) = highEnergy4(m,j)/totalEnergy;
    end;
end;

%figure(1);
%hold on;
%for m = 1:Ns,
%    plot(1:Nd,energyFraction(m,1:Nd),'color',Color{m});
%end;
%title ('energy fraction in particles with {\gamma} > {\gamma_0}');
%xlabel ('Nt');
%ylabel ('{\epsilon}');
%legend('20 rg','10 rg','5 rg', '2.5 rg','Location','southwest');
%grid ;

figure(2);
hold on;
for m = 1:Ns,
    plot(1:Nd,energyFraction1(m,1:Nd),'color',Color{m});
    %plot(1:Nd,highEnergy1(m,1:Nd)/width(m),'color',Color{m});
end;
title ('energy fraction in particles with {\gamma} > 15');
xlabel ('Nt');
ylabel ('{\epsilon}');
legend('20 rg','10 rg','5 rg', '2.5 rg','Location','southwest');
grid ;

figure(3);
hold on;
for m = 1:Ns,
    plot(1:Nd,energyFraction2(m,1:Nd),'color',Color{m});
    %plot(1:Nd,highEnergy2(m,1:Nd)/width(m),'color',Color{m});
end;
title ('energy fraction in particles with {\gamma} > 50');
xlabel ('Nt');
ylabel ('{\epsilon}');
legend('20 rg','10 rg','5 rg', '2.5 rg','Location','southwest');
grid ;

figure(4);
hold on;
for m = 1:Ns,
    plot(1:Nd,energyFraction3(m,1:Nd),'color',Color{m});
    %plot(1:Nd,highEnergy3(m,1:Nd)/width(m),'color',Color{m});
end;
title ('energy fraction in particles with {\gamma} > 200');
xlabel ('Nt');
ylabel ('{\epsilon}');
legend('20 rg','10 rg','5 rg', '2.5 rg','Location','southwest');
grid ;

figure(5);
hold on;
for m = 1:Ns,
    plot(1:Nd,energyFraction4(m,1:Nd),'color',Color{m});
    %plot(1:Nd,highEnergy4(m,1:Nd)/width(m),'color',Color{m});
end;
title ('energy fraction in particles with {\gamma} > 500');
xlabel ('Nt');
ylabel ('{\epsilon}');
legend('20 rg','10 rg','5 rg', '2.5 rg','Location','southwest');
grid ;