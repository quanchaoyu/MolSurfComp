function output_SES_toroidalpatches2(sing,A1,A2,THETA,SAScircle,circle1,circle2,circle1_label,circle2_label,arg_eSES,arg_SAScircle,IJ)
% c
% c   write the file of toroidal SES patches
% c   Remark: double-triangle-shaped or rectangle-shaped
% c

global index_toroidalpatch;

if sing == 1
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,A1,THETA(1),SAScircle,1,circle1,circle1_label,arg_eSES,arg_SAScircle,IJ);
    
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,A2,THETA(2),SAScircle,1,circle2,circle2_label,arg_eSES,arg_SAScircle,IJ);
else
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,[],THETA(3),SAScircle,2,[circle1;circle2],[circle1_label;circle1_label],arg_eSES,arg_SAScircle,IJ);
end


end

function write_patch(index_toroidalpatch,sing,cusp,theta,SAScircle,ncircle_toroidal,circle_toroidal,circle_label,arg_eSES,arg_SAScircle,IJ)
global Para;

ID = ['./',['SEStor_',Para.filename],'.txt'];

if index_toroidalpatch == 1
    fileID = fopen(ID,'wt');
    fprintf(fileID,'# TOROIDAL SES-PATCHES\n');
    fprintf(fileID,'# Cusp Point: singularity coordinates');
    fprintf(fileID,'# Corresponding SAS-Segment: [c,n,r,starting,ending,radian]\n');
    fprintf(fileID,'# 4 or 3 Boundary Segments: [c,n,r,starting,ending,radian,left-hand side]\n');
    fprintf(fileID,'# 2 Boundary Circles: [c,n,r]\n\n');
    
else
    fileID = fopen(ID,'a');
end

fprintf(fileID,'PATCH %d  ',index_toroidalpatch);
if arg_eSES == 1
    fprintf(fileID,'Exterior ');
else
    fprintf(fileID,'Interior ');
end

if sing == 1
    fprintf(fileID,'Singular\n');
    fprintf(fileID,'Cusp Point: %10f %10f %10f\n',cusp'); % the coordinate of the cusp
else
    fprintf(fileID,'Non-Singular\n'); 
end

fprintf(fileID,'Corresponding Two Atoms\n');
fprintf(fileID,'%10d %10d\n',IJ');

fprintf(fileID,'Corresponding SAS-Circle\n');
fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f\n',SAScircle');% the corresponding SAS circle  

fprintf(fileID,'Circles: %d\n',ncircle_toroidal);
fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10d\n',circle_toroidal');

fprintf(fileID,'\n');

fclose(fileID);

end