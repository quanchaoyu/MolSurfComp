function output_SES_convexpatches(index_atom,index_convexpatch,extSES,c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize,circleindex)
% c
% c   write the file of convex SES patches
% c   Remark: the patch has the boundary composed of loops and circles
% c

global Para;
Rp = Para.radius_probe;

maxnum = max(loopsize)*patchesize;

segment_convex = zeros(maxnum,15); % [c,n,r,startpoint,endpoint,angle,right-hand], right-hand: 1 or 0
nsegment_convex = 0;

circle_convex = zeros(maxnum,8); % [c,n,r,right-hand]
ncircle_convex = 0;

loops_convex = zeros(patchesize,2); % loops_convex(i,:) = [index_startingsegment,index_endingsegment]


for i = 1:patchesize
    loopsize_convex = 0;
    nloops_convex = 0;
    if patches(i) > 0
        nloops_convex = nloops_convex+1;
        for j = 1:loopsize(patches(i))
            nsegment_convex = nsegment_convex+1;
            loopsize_convex = loopsize_convex+1;

            ij = loops(patches(i),j);

            c = segment(ij,1:3);
            n = segment(ij,4:6);
            r = segment(ij,7);
            startpoint = segment(ij,8:10);
            u=(startpoint-c)/r;
            v=[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];
            angle = segment(ij,11);
            endpoint = r*cos(angle)*u+r*sin(angle)*v+c;

            c1 = c_sphere+(c-c_sphere)*(r_sphere-Rp)/r_sphere;
            r1 = r*(r_sphere-Rp)/r_sphere;
            startpoint1 = c_sphere+(startpoint-c_sphere)*(r_sphere-Rp)/r_sphere;
            endpoint1 = c_sphere+(endpoint-c_sphere)*(r_sphere-Rp)/r_sphere;

            segment_convex(nsegment_convex,:) = [c1,n,r1,startpoint1,endpoint1,angle,1];

        end

        loops_convex(i,1:2) = [nsegment_convex-loopsize_convex+1,nsegment_convex]; 
    else
        ncircle_convex = ncircle_convex+1;
        c = circle(-patches(i),1:3);
        n = circle(-patches(i),4:6);
        r = circle(-patches(i),7);
        
        c2 = c_sphere+(c-c_sphere)*(r_sphere-Rp)/r_sphere;
        r2 = r*(r_sphere-Rp)/r_sphere;
        
        circle_convex(ncircle_convex,:) = [c2,n,r2,1];
    end
end

segment_convex = segment_convex(1:nsegment_convex,:);

circle_convex = circle_convex(1:ncircle_convex,:);

write_patch(index_convexpatch,extSES,c_sphere,r_sphere,segment_convex,nsegment_convex,loops_convex,nloops_convex,circle_convex,ncircle_convex)

end

function write_patch(index_convexpatch,extSES,c_sphere,r_sphere,segment_convex,nsegment_convex,loops_convex,nloops_convex,circle_convex,ncircle_convex)
global Para;
Rp = Para.radius_probe;

ID = ['./',['SESconvex_',Para.filename],'.txt'];

if index_convexpatch == 1
    fileID = fopen(ID,'wt');
    fprintf(fileID,'# CONVEX SES-PATCHES\n');
    fprintf(fileID,'# Circular Segment: [c,n,r,starting,ending,radian]\n');
    fprintf(fileID,'# Free Circule: [c,n,r]\n\n');
else
    fileID = fopen(ID,'a');
end

    fprintf(fileID,'\nPATCH ');
    fprintf(fileID,'%d  ',index_convexpatch);% index of the convex patch
    
    if extSES
        fprintf(fileID,'Exterior\n');% if the patch is on exterior SES, extSES = 1; otherwise, extSES = 0
    else
        fprintf(fileID,'Interior\n');
    end
    
    fprintf(fileID,'CenterRadius\n');
    fprintf(fileID,'%10f %10f %10f %10f\n',[c_sphere,r_sphere]');% center and radius of the VdW-ball
    
    fprintf(fileID,'Loops: %d\n',nloops_convex);
    if nloops_convex > 0
        fprintf(fileID,'%4d %4d\n',loops_convex'); % the interior of the patch is on the right-hand side of each loop
    end
    
    fprintf(fileID,'Segments: %d\n',nsegment_convex);
    if nsegment_convex > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10f %10d\n',[segment_convex]');
        % the common segments between two neighbor patches have same label
    end
    
    fprintf(fileID,'Circles: %d\n',ncircle_convex);
    if ncircle_convex > 0
        fprintf(fileID,'%10f %10f %10f %10f %10f %10f %10f %10d \n',[circle_convex]'); 
        % the common circles between two neighbor patches have same label
    end
    
    fclose(fileID);

end