function data_Seg_Pat
%c
%c   Compute Segments, Loops & Patches on the SAS
%c
global Para Geom Inter;
M = Geom.M;
C = Geom.centers;
R = Geom.R;
Rp = Para.radius_probe;

inter = Inter.M_int; %  intersection matrix
Row = Inter.num_int; %  intersection number for each SAS-ball

global DataI DataSeg DataCir DataLoop DataPat;
I = DataI.I;
s = DataI.nI;
Iijk = DataI.Iijk;
In = DataI.In;
direction = DataI.direction;
I_circle = DataI.I_circle;
I_circle_num = DataI.I_circle_num;

circle = DataCir.circle;
ncircle = DataCir.ncircle;
circleindex = DataCir.circleindex;
ncircleindex = DataCir.ncircleindex;

nsegment = 0; % number of SAS circular segments
segment = zeros(s*3/2,5);% segment(k) -- [i,j,x1,x2,clockwise (1) or backclockwise (-1)], x1 and x2 are indices of starting and ending points of a segment
ncrasegment = zeros(s*3/2,8);% ncrasegment(k) -- [normal vector, center, radius, radian]


smax = max(In); % the max number of intersection points (segments) on an atom; the max number of loops on an atom
satom = zeros(M,smax); % satom(i,:) represents all the segments on the ith atom 
nsatom = zeros(M,1); % nsatom(i) represents the number of segments on the ith atom

lmax = s*3/2+ncircle; % the max number of loops on the SAS
loops = zeros(lmax,smax);% loops(i,:) records all the segments on the i-th loop
loopsize = zeros(lmax,1); % loopsize(i) is the number of segments on the i-th loop
nloops = 0; % number of loops

pmax = 10; % the max number of loops on a spherical patch
patches=zeros(lmax,pmax); % patches(i,:) records the loops on the boundary of the i-th patch
patchesize=zeros(lmax,1); 
npatches=0; % number of patches

loops_index = zeros(M,2); % loops_index(i,:) = [starting loop-index, ending loop-index] on the i-th sphere
patches_index = zeros(M,2); % patches_index(i,:) = [starting patch-index, ending patch-index] on the i-th sphere
patch_atom=zeros(lmax,1); % patch_atom(i) is the indice of the SAS-sphere where the i-th loop lies



%% Loops and Patches Construction for MSAS (i.e. cSAS)
global DataAV;
DataAV.Acsas = 0;
DataAV.Acses = 0;
DataAV.Aesas = 0;
DataAV.Aeses = 0;
DataAV.Vcsas = 0;
DataAV.Vcses = 0;
DataAV.Vesas = 0;
DataAV.Veses = 0;


Area_sphpat = zeros(lmax,1);% Area_sphpat(i) records the area of the i-th spherical SAS-patch


global tree;
max_nloops = 20;
for j = 1:(2*max_nloops+1)
    tree0(j) = treenode(0,[],0,0,0); 
end
tree = tree0;

global Rj;
Rj=zeros(nsegment,1);
segment0 = zeros(nsegment,11);

for i=1:M

    for row=1:Row(i)
        j=inter(i,row);
        
        if j>i && I_circle_num(i,row)>0
            A=circlecenter(C(i,:),C(j,:),R(i)+Rp,R(j)+Rp);
            rij=sqrt((R(i)+Rp)^2-norm(C(i,:)-A)^2);
            nij=(C(j,:)-C(i,:))/norm(C(j,:)-C(i,:));
            point=I_circle(i,row,1:I_circle_num(i,row)); % collect all points on the intersection circle
            spoint=I_circle(i,row,1); % start point
            
            index=Iijk(spoint,:);
            if i==index(1)&& j==index(2)
                direct=direction(spoint,1);
            elseif i==index(2)&&j==index(3)
                direct=direction(spoint,2);
            elseif i==index(1)&&j==index(3)
                direct=direction(spoint,3);
            end
            
            alpha1=zeros(1,I_circle_num(i,row));
            for k=2:I_circle_num(i,row)
                alpha1(k)=alpha(direct,I(spoint,:)-A,I(point(k),:)-A,nij);
            end
            [alpha2,order]=sort(alpha1);
            point=point(order);
                        
            for k=1:I_circle_num(i,row)/2
                nsegment=nsegment+1;
                segment(nsegment,:)=[i,j,point(2*k-1),point(2*k),direct];% j>i; direct=1 means clockwise,direct=-1 means counterclockwise
                ncrasegment(nsegment,:)=[nij,A,rij,alpha2(2*k)-alpha2(2*k-1)];

                nsatom(i)=nsatom(i)+1; % number of segments on the i-th atom
                satom(i,nsatom(i))=nsegment; % segments on the i-th atom
                nsatom(j)=nsatom(j)+1;
                satom(j,nsatom(j))=nsegment;
            end
        end
    end
    
    if nsatom(i)>0
        %c construct loops on the i-th SAS-ball
        [loops_i,loopsize_i,nloops_i]=loopconstruct(i,satom(i,:),nsatom(i),segment,smax);
        loops(nloops+1:nloops+nloops_i,:)=loops_i;
        loopsize(nloops+1:nloops+nloops_i,1)=loopsize_i;
        nloops=nloops+nloops_i;
        loops_index(i,:)=[nloops-nloops_i+1,nloops]; % the loops on the i-th sphere.
    elseif ncircleindex(i)>0
        loops_i=[];
        loopsize_i=[];
        nloops_i=0;
    end
    
    if nsatom(i)>0||ncircleindex(i)>0
        
        [patches_i,patchesize_i,npatches_i]=patchesconstruct(i,C,R,segment,ncrasegment,pmax,I,circle,circleindex(i,:),ncircleindex(i),loops_i,loopsize_i,nloops_i,nloops);
        
        patches(npatches+1:npatches+npatches_i,:)=patches_i;
        patchesize(npatches+1:npatches+npatches_i,1)=patchesize_i;
        npatches=npatches+npatches_i;
        patches_index(i,:)=[npatches-npatches_i+1,npatches];
        patch_atom(npatches-npatches_i+1:npatches,1)=ones(npatches_i,1).*i;        
        
        % modify loops_i,segment,circle so that visu_sphereicalpatch can be used
        % segment0 : [c,n,r,spoint,angle] segment0,
        [loops_i0,segment0,circle0] = mod_seg_loop_cir(i,nloops_i,loops_i,loopsize_i,circleindex,ncircleindex,satom,nsatom,C,R,segment,nsegment,ncrasegment,segment0);
        
        for j = 1:npatches_i
            
            Area = area_spherical(C(i,:),R(i)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j));
            Area_sphpat(npatches-npatches_i+j) = Area;
            
            DataAV.Acsas = DataAV.Acsas+Area;
            DataAV.Acses = DataAV.Acses+(R(i)/(R(i)+Rp))^2*Area;
            
            DataAV.Vcsas=DataAV.Vcsas+Area*(R(i)+Rp)/3;
            DataAV.Vcses=DataAV.Vcses+(R(i)/(R(i)+Rp))^2*Area*R(i)/3;
            
        end
        
    end

end


DataSeg.segment = segment;
DataSeg.nsegment = nsegment;
DataSeg.ncrasegment = ncrasegment;
DataSeg.satom = satom;
DataSeg.nsatom = nsatom;

DataLoop.loops = loops;
DataLoop.nloops = nloops;
DataLoop.loopsize = loopsize;
DataLoop.loops_index = loops_index;

DataPat.patches = patches;
DataPat.npatches = npatches;
DataPat.patchesize = patchesize;
DataPat.patches_index = patches_index;
DataPat.patch_atom = patch_atom; % patch index to atom index

DataAV.Area_sphpat = Area_sphpat;
end

function [loops,loopsize,nloops]=loopconstruct(i,sa,ns,segment,smax)
%c
%c  Construct all loops on the i-th sphere
%c  Remark: segment(nsegment,1:5) is composed of [i,j,p1,p2,direction]
%c
s=zeros(ns,2);
for k=1:ns
    if segment(sa(k),5)*(2*(segment(sa(k),1)==i)-1) == 1 %%%modify the direction of the loop, to ensure the lefthand is interior
        s(k,:)=[segment(sa(k),4),segment(sa(k),3)];
    else
        s(k,:)=segment(sa(k),3:4);
    end
end

true=zeros(ns,1);
number=0;
nloops=0;
loops=zeros(floor(ns/2),smax);
loopsize=zeros(floor(ns/2),1);

for k1=1:ns
    
   loop=zeros(1,smax);
   loopsi=0;
   if true(k1)==0 && number<ns %%% ns means the number of segments on the ith atom
       nloops=nloops+1;
       pointer=k1;
       loopsi=loopsi+1;
       number=number+1; %%% number represents the already used segments.
       loop(loopsi)=sa(pointer);
       true(k1)=1;
       
       while ( s(k1,1)~=s(pointer,2) )
           
           for k2=k1+1:ns
               if true(k2)==0 && s(k2,1)==s(pointer,2)
                   pointer=k2;
                   loopsi=loopsi+1;
                   number=number+1;
                   loop(loopsi)=sa(pointer);
                   
                   true(pointer)=1;
               end
           end
           
       end
       
       loops(nloops,:)=loop;
       loopsize(nloops)=loopsi;
   end

end

loops=loops(1:nloops,:);
loopsize=loopsize(1:nloops);

end

function [patches,patchesize,npatches]=patchesconstruct(i,C,R,segment,ncrasegment,pmax,I,circle,circleindex,ncircleindex,loops_i,loopsize_i,nloops_i,nloops)
%c
%c  Construct all spherical patches on the i-th sphere
%c  Use the binary tree method to build boundaries of patches by classifying loops
%c

global Para;
Rp = Para.radius_probe;

patches=zeros(nloops_i+ncircleindex,pmax); %% it records the boundary loops of spherical patches.
patchesize=zeros(nloops_i+ncircleindex,1);
npatches=0;

global tree;
S0=[1:1:nloops_i,-1:-1:-ncircleindex];
tree(1)=treenode(1,S0,S0(1),nloops_i,ncircleindex);
ntree=1;
j=1;

if ncircleindex >0 && nloops_i>= 2
    
end

for s=1:2*(nloops_i+ncircleindex)+1
    
    if j>ntree%%this means that the set has been completed divided
        break;
    end
    
    if j<=ntree && tree(j).activenode==1
        S1=[];
        S2=[];
        k1=0;
        k2=0;
        t1=1;
        t2=1;
        left_n1=0;
        left_n2=0;
        right_n1=0;
        right_n2=0;
        
        S=tree(j).set;
        k=tree(j).activeelement;% k denotes the fixed loop
        
        if length(S)==1
            tree(j).activenode=0;
            j=j+1;
            continue;
        end
        
        if k>0
            for s1=1:tree(j).n1
                point=I(segment(loops_i(S(s1),1),3),:);
                
                
                ci=C(i,:);
                ri=R(i)+Rp; %%% ????

                if k==S(s1) || interiorloop(point,ci,C,ri,i,loops_i(k,:),loopsize_i(k,:),I,segment,ncrasegment) 
                    S1=[S1,S(s1)];
                    left_n1=left_n1+1;
                    if t1==1 && S(s1)>k
                        t1=0;
                        k1=S(s1);
                    end
                else
                    S2=[S2,S(s1)];
                    right_n1=right_n1+1;
                    if t2==1&& S(s1)>k
                        t2=0;
                        k2=S(s1);
                    end
                end 
            end

            for s2=(tree(j).n1+1):(tree(j).n1+tree(j).n2)
                    [vector1,~]=orthogonalvectors(circle(circleindex(-S(s2)),6:8));
                    point=circle(circleindex(-S(s2)),3:5)+circle(circleindex(-S(s2)),9)*vector1;
                    ci=C(i,:);
                    ri=R(i)+Rp;
                    
                    if interiorloop(point,ci,C,ri,i,loops_i(k,:),loopsize_i(k,:),I,segment,ncrasegment)
                        S1=[S1,S(s2)];
                        left_n2=left_n2+1;
                        if t1==1
                            t1=0;
                            k1=S(s2);
                        end
                    else
                        S2=[S2,S(s2)];
                        right_n2=right_n2+1;
                        if t2==1
                            t2=0;
                            k2=S(s2);
                        end
                    end 
            end
            
            
            if length(S2)==0
                
                if t1==1
                    tree(j).activenode=0;
                    j=j+1;
                else
                    tree(j).activeelement=k1;
                end
                
            else
                tree(j).activenode=0;
                tree(j).leftnode=ntree+1;
                tree(j).rightnode=ntree+2;
                
                ntree=ntree+1;
                tree(ntree)=treenode(1,S1,k1,left_n1,left_n2);
                if t1==1
                    tree(ntree).activenode=0;
                end
                
                ntree=ntree+1;
                tree(ntree)=treenode(1,S2,k2,right_n1,right_n2);
                if t2==1
                    tree(ntree).activenode=0;
                end
                
                j=j+1;
            end
            
        else
            tree(j).activenode=0;
            j=j+1;
        end
        %%% Consider the next treenode.
    else
        j=j+1;
    end
    
end


for j=1:ntree
    if tree(j).activenode==0 && tree(j).leftnode==0
        npatches=npatches+1;
        patchesize(npatches,1)=tree(j).n1+tree(j).n2;
        patches(npatches,1:patchesize(npatches,1))=tree(j).set;
    end
end



patches=patches(1:npatches,:);
patchesize=patchesize(1:npatches,1);


end

function alpha=alpha(direct,u,v,n)
%c
%c  Compute the angle between two neighbor edges e and f, countorclockwise direction
%c
t=sign(det([u;v;n]));%direct=1 means clockwise
if direct*t>0
    alpha=acos(u*v'/(norm(u)*norm(v)));
else
    alpha=2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end

function O=circlecenter(c1,c2,r1,r2)
%c
%c  Compute the center of an intersection circle
%c
d=sqrt((c1(1)-c2(1))^2+(c1(2)-c2(2))^2+(c1(3)-c2(3))^2);
t=(r1^2-r2^2+d^2)/(2*d);
O=c1+(c2-c1)*t/d;
end