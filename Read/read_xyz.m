function Geom = read_xyz(filename)
% c
% c   Read xyz-file and save data to Geom
% c

scalingfactor = 1.1;

filename = [filename,'.xyz'];
fid = fopen(filename,'r');

line = fgetl(fid);
it = 0;
while ischar(line)
    if it == 0
        Geom.M = sscanf(line,'%d');
        it = it+1;
        line = fgetl(fid);
        continue
    end
    
    [token,remain] = strtok(line); 
    c =  sscanf(remain,'%f %f %f');
    switch(token)
        %% UFF VdW radius table, Charges
        case  'H', r = 1.443; 
        case  'C', r = 1.9255; 
        case  'O', r = 1.75; 
        case  'N', r = 1.83; 
        case  'P', r = 2.0735; 
        case  'S', r = 2.0175; 
        case 'Mg', r = 1.5105;
        otherwise, line = fgetl(fid); continue;
    end
        
    Geom.R(it,1) = scalingfactor*r;
    Geom.centers(it,:) = c';
    
    it = it + 1;
    line = fgetl(fid);
end
it = it-1;
fclose(fid);
Geom.centers = Geom.centers(1:it,:);
Geom.M = it;

end
