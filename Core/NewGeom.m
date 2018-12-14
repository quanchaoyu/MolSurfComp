function NewGeom(AddSpheres)
%
%   Input: Added Spheres including charges, centers and radii
%   Output new geometry with added spheres
%
global Geom_new;
add = 1;

M_add = 0;
C_add = [];
R_add = [];
    
if add == 1
    M_add = size(AddSpheres,1);
    C_add = AddSpheres(:,1:3);
    R_add = AddSpheres(:,4);
end


Geom_new.M = M_add;
Geom_new.centers = C_add;
Geom_new.R = R_add;

output_Geom_new;
end

function output_Geom_new
global Para;
global Geom_new;

ID = ['./',['AddSpheres_',Para.filename],'.txt'];

fileID = fopen(ID,'wt');

fprintf(fileID,'Add Spheres for Filling Inner Holes\n');
fprintf(fileID,'Number: %d \n',Geom_new.M);
fprintf(fileID,'%10c%10c%10c%10c\n',['R','X','Y','Z']);
formatSpec = '%10f%10f%10f%10f\n';

A = [Geom_new.R,Geom_new.centers];
fprintf(fileID,formatSpec,A');
fprintf(fileID,'\n');

fclose(fileID);

end