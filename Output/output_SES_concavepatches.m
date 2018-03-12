function output_SES_concavepatches(atoms,NB,index_I_segment,index_I_circle,index_I,index_concavepatch,arg_eSES,c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize)
% c
% c   write the file of external or complete concave SES patches
% c   Remark: the patch has the boundary composed of loops and circles
% c

maxnum = max(loopsize)*patchesize;

segment_concave = zeros(maxnum,15); % [c,n,r,startpoint,endpoint,angle,right-hand]
segment_label = zeros(maxnum,5); % [type = 2,index_I,atom_i,atom_j,-1, 0 or 1] or [type = 3,index_I,index_I1,angle_middlepoint,0]
nsegment_concave = 0;

circle_concave = zeros(maxnum,8); % [c,n,r,right-hand]
circle_label = zeros(maxnum,3); % [type = 3,index_I,index_I1]
ncircle_concave = 0;

loops_concave = zeros(patchesize,2); % loops_concave(i,:) = [index_startingsegment,index_endingsegment]

for i = 1:patchesize
    loopsize_concave = 0;
    nloops_concave = 0;
    if patches(i) > 0
        nloops_concave = nloops_concave+1;
        for j = 1:loopsize(patches(i))
            nsegment_concave = nsegment_concave+1;
            loopsize_concave = loopsize_concave+1;

            ij = loops(patches(i),j);

            c = segment(ij,1:3);
            n = segment(ij,4:6);
            r = segment(ij,7);
            startpoint = segment(ij,8:10);
            u=(startpoint-c)/r;
            v=[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];
            angle = segment(ij,11);
            endpoint = r*cos(angle)*u+r*sin(angle)*v+c;
            
            segment_concave(nsegment_concave,:) = [c,n,r,startpoint,endpoint,angle,-1];
           
        end

        loops_concave(i,1:2) = [nsegment_concave-loopsize_concave+1,nsegment_concave]; 
    else
        ncircle_concave = ncircle_concave+1;
        c = circle(-patches(i),1:3);
        n = circle(-patches(i),4:6);
        r = circle(-patches(i),7);
        
        circle_concave(ncircle_concave,:) = [c,n,r,-1];
    end
end

segment_concave = segment_concave(1:nsegment_concave,:);

circle_concave = circle_concave(1:ncircle_concave,:);

sing = 0;
if patchesize > 1 || loopsize(patches(1))>3
   sing = 1; 
end

if arg_eSES == 1
    write_eSESpatch(sing,index_concavepatch,c_sphere,r_sphere,segment_concave,nsegment_concave,loops_concave,nloops_concave,circle_concave,ncircle_concave);
else
    write_cSESpatch(sing,index_concavepatch,c_sphere,r_sphere,segment_concave,nsegment_concave,loops_concave,nloops_concave,circle_concave,ncircle_concave);
end

end


function write_eSESpatch(sing,index_concavepatch,c_sphere,r_sphere,segment_concave,nsegment_concave,loops_concave,nloops_concave,circle_concave,ncircle_concave)
global Para;

ID = ['./',['eSESconcave_',Para.filename],'.txt'];

if index_concavepatch == 1
    fileID = fopen(ID,'wt');
    fprintf(fileID,'# CONCAVE eSES-PATCHES\n');
    fprintf(fileID,'# SAS Circular Segment: [c,n,r,starting,ending,radian]\n');
    fprintf(fileID,'# SAS Free Circule: [c,n,r]\n\n');
else
    fileID = fopen(ID,'a');
end
    fprintf(fileID,'PATCH ');
    fprintf(fileID,'%d  ',index_concavepatch);% index of the concave patch
    
    if sing
        fprintf(fileID,'Singular\n');% if the patch is on exterior SES, extSES = 1; otherwise, extSES = 0
    else
        fprintf(fileID,'Non-Singular\n');
    end

    fprintf(fileID,'CenterRadius\n');
    fprintf(fileID,'%10f %10f %10f %10f\n',[c_sphere,r_sphere]);% center and radius of the probe
    
    fprintf(fileID,'Loops: %d\n',nloops_concave);
    if nloops_concave > 0
        fprintf(fileID,'%4d %4d\n',loops_concave'); % the interior of the patch is on the right-hand side of each loop
    end
    
    fprintf(fileID,'Segments: %d\n',nsegment_concave);
    if nsegment_concave > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10d \n',segment_concave');
    end
    
    fprintf(fileID,'Circles\n');
    fprintf(fileID,'%10d\n',ncircle_concave);
    if ncircle_concave > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10d \n',circle_concave');
    end
    
    fprintf(fileID,'\n');
    fclose(fileID);

end

function write_cSESpatch(sing,index_concavepatch,c_sphere,r_sphere,segment_concave,nsegment_concave,loops_concave,nloops_concave,circle_concave,ncircle_concave)
global Para;

ID = ['./',['cSESconcave_',Para.filename],'.txt'];

if index_concavepatch == 1
    fileID = fopen(ID,'wt');
    fprintf(fileID,'# CONCAVE cSES-PATCHES\n');
    fprintf(fileID,'# Circular Segment: [c,n,r,starting,ending,radian]\n');
    fprintf(fileID,'# Free Circule: [c,n,r]\n\n');
else
    fileID = fopen(ID,'a');
end
    fprintf(fileID,'PATCH ');
    fprintf(fileID,'%d  ',index_concavepatch);% index of the concave patch
    
    if sing
        fprintf(fileID,'Singular\n');% if the patch is on exterior SES, extSES = 1; otherwise, extSES = 0
    else
        fprintf(fileID,'Non-Singular\n');
    end

    fprintf(fileID,'CenterRadius\n');
    fprintf(fileID,'%10f %10f %10f %10f\n',[c_sphere,r_sphere]);% center and radius of the probe
    
    fprintf(fileID,'Loops: %d\n',nloops_concave);
    if nloops_concave > 0
        fprintf(fileID,'%4d %4d\n',loops_concave'); % the interior of the patch is on the right-hand side of each loop
    end
    
    fprintf(fileID,'Segments: %d\n',nsegment_concave);
    if nsegment_concave > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10d \n',segment_concave');
    end
    
    fprintf(fileID,'Circles\n');
    fprintf(fileID,'%10d\n',ncircle_concave);
    if ncircle_concave > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10d \n',circle_concave');
    end
    
    fprintf(fileID,'\n');
    fclose(fileID);

end