function diff = kappafit( x, Pe, Fe, Np)
diff = 0;
mp = 1.67*10^-24;
mass_ratio = 100;
me = mp/mass_ratio;
c = 2.99792458*10^10;
Fkappa(1:Np) = 0;
for i = 1:Np,
    Fkappa(i) = (Pe(i)^4)*(1 + (sqrt(1.0 + Pe(i)*Pe(i)) - 1.0)/(x(1)*x(1)*x(2)))^(-(x(2) + 1));
end;

normkappae = (Fkappa(1)/(Pe(2)^2))*(Pe(2) - Pe(1));

for i = 2:Np,
    normkappae = normkappae + (Fkappa(i)/(Pe(i)^2))*(Pe(i) - Pe(i-1));
end;
for i = 1:Np,
    Fkappa(i) = Fkappa(i)/normkappae;
end;

index1 = 60;
index2 = 185;

diff = 0;
for i = index1:index2,
    diff = diff + ((Fkappa(i) - Fe(i))^2)*(Pe(i) - Pe(i-1))/(Pe(i)^2);
end;
end

