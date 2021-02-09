clear;
data = importdata('output/image.dat');
Nr = size(data,1)/4;
Nphi = size(data,2);
image(1:Nphi+1,1:Nr) = 0;
for i = 1:Nr,
    for j = 1:Nphi,
        image(j,i) = data(i,j);
    end;
    image(Nphi+1,i) = data(i,1);
end;


figure(1)
[r,t] = meshgrid(1:Nr,pi/Nphi:2*pi/Nphi:2*pi+pi/Nphi);
x = r.*cos(t);
y = r.*sin(t);
contourf(x,y,image);
hold on
plot([zeros(1,13); 90*cosd(0:30:360)], [zeros(1,13); 90*sind(0:30:360)],'k')
plot(90*((0:0.33:1)'*cosd(0:10:360))', 90*((0:0.33:1)'*sind(0:10:360))','k')
colorbar
set(colorbar,'FontSize',16)
axis equal
%set(gca, 'Box','off', 'XColor','none', 'YColor','none',  'Color','none')
hold off