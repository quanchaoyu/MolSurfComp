function Geom = read_PDB(filename)
% c
% c   Read pdb-file and save data to Geom
% c

filename = [filename,'.pdb'];
fid = fopen(filename,'r');

scalingfactor = 1.1;

line = fgetl(fid);
it = 0;
while ischar(line)
	if (strncmp('HETATM',line,6) || strncmp('ATOM',line,4)) 
        it = it + 1;
        switch(line(14))
            %% UFF VdW radius table, Charges
            case  'H', r = 1.443; 
            case  'C', r = 1.9255; 
            case  'O', r = 1.75; 
            case  'N', r = 1.83; 
            case  'P', r = 2.0735; 
            case  'S', r = 2.0175; 
            otherwise, r = 1.75; 
        end
        
        c =  sscanf(line(31:54),'%f %f %f');
        
        Geom.R(it,1) = scalingfactor*r;
        Geom.centers(it,:) = c;
        Geom.T(it,1) = line(14);
        
    end
    
    line = fgetl(fid);
end
fclose(fid);

Geom.M = it;
Geom.centers = Geom.centers(1:it,:);

end