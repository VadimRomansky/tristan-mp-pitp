clear;
directory_name = './output/';
file_name = 'flds.tot';
part_name = 'prtl.tot';
file_number = '.001';
full_name = strcat(directory_name, file_name, file_number);
full_part_name = strcat(directory_name, part_name, file_number);
Bx = hdf5read(full_name,'bx');
By = hdf5read(full_name,'by');
Bz = hdf5read(full_name,'bz');
fileinfo = hdf5info(full_part_name);
last_number = 1;
a = last_number;

if(a < 10)
        full_name = strcat(directory_name, file_name, '.00', num2str(a));
        full_part_name = strcat(directory_name, part_name, '.00', num2str(a));
    else if (a < 100)
            full_name = strcat(directory_name, file_name, '.0', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.0', num2str(a));  
        else 
            full_name = strcat(directory_name, file_name, '.', num2str(a));
            full_part_name = strcat(directory_name, part_name, '.', num2str(a));
        end;
    end;

gammae = hdf5read(full_part_name, 'gammae');
inde = hdf5read(full_part_name, 'inde');
proce = hdf5read(full_part_name, 'proce');
xe = hdf5read(full_part_name, 'xe');
ye = hdf5read(full_part_name, 'ye');
ze = hdf5read(full_part_name, 'ze');
ve = hdf5read(full_part_name, 've');
ue = hdf5read(full_part_name, 'ue');
we = hdf5read(full_part_name, 'we');
gammai = hdf5read(full_part_name, 'gammai');
indi = hdf5read(full_part_name, 'indi');
xi = hdf5read(full_part_name, 'xi');
yi = hdf5read(full_part_name, 'yi');
zi = hdf5read(full_part_name, 'zi');
%Ex = hdf5read(full_name,'ex');
%Ey = hdf5read(full_name,'ey');
%Ez = hdf5read(full_name,'ez');
