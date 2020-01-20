clear;
directory_name = './output/';
part_name = 'prtl.tot';
file_number = '.000';
full_part_name = strcat(directory_name, part_name, file_number);

gammae = hdf5read(full_part_name, 'gammae');
inde = hdf5read(full_part_name, 'inde');
proce = hdf5read(full_part_name, 'proce');
xe = hdf5read(full_part_name, 'xe');
ye = hdf5read(full_part_name, 'ye');
ze = hdf5read(full_part_name, 'ze');
ve = hdf5read(full_part_name, 've');
ue = hdf5read(full_part_name, 'ue');
we = hdf5read(full_part_name, 'we');

u = 0;
v = 0;
w = 0;

u2 = 0;
v2 = 0;
w2 = 0;

Np = size(gammae,1);
for i = 1:Np,
    if(xe(i) < 10000)
        u = u + ue(i)/Np;
        u2 = u2 + ue(i)*ue(i)/Np;
        v = v + ve(i)/Np;
        v2 = v2 + ve(i)*ve(i)/Np;
        w = w + we(i)/Np;
        w2 = w2 + we(i)*we(i)/Np;
    end;
end;