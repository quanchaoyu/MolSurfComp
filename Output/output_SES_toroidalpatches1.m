function output_SES_toroidalpatches1(sing,A1,A2,THETA,SASsegment,segment1,segment2,segment1_label,segment2_label,segment0,segment0_label,arg_eSES,arg_SAScircle,IJ)
% c
% c   write the file of toroidal SES patches
% c   Remark: double-triangle-shaped or rectangle-shaped
% c

global index_toroidalpatch;

if sing == 1 % singular patch
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,A1,THETA(1),SASsegment,3,segment1,segment1_label,arg_eSES,arg_SAScircle,IJ);
    
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,A2,THETA(2),SASsegment,3,segment2,segment2_label,arg_eSES,arg_SAScircle,IJ);
else
    index_toroidalpatch = index_toroidalpatch+1;
    write_patch(index_toroidalpatch,sing,[],THETA(3),SASsegment,4,segment0,segment0_label,arg_eSES,arg_SAScircle,IJ); % segment_label = [1,atom,nsegment,0,0] or [2,index_I,atom_i,atom_j,-1 or 0 or 1]
end


end

function write_patch(index_toroidalpatch,sing,cusp,theta,SASsegment,nsegment_toroidal,segment_toroidal,segment_label,arg_eSES,arg_SAScircle,IJ)
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

fprintf(fileID,'Corresponding SAS-Segment\n');
fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f\n',SASsegment');% the corresponding SAS segment  

fprintf(fileID,'Segments: %d\n',nsegment_toroidal);
%% [c,n,r,s,e,angle], n points from ci to cj
fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10d\n',segment_toroidal');

fprintf(fileID,'\n');

fclose(fileID);

end