function construct_concavepat(k1,atoms,I,I_probe,N_probe,segment,N_segment,circle_interior,N_circle,direct)
%c
%c   Construct the concave patch associated with the k1-th probe, based on its data structure
%c   Remark: segment=[c,n,r,a,j1,j2], n depends on direct (pointing x if direct=1,otherwise x1), [j1,j2] are indices of I_probe;
%c           circle=[c,n,r], n depends on direct
%c

epsilon = 10^-10;
global Para;
Rp = Para.radius_probe;

global map_same;
nsamepoint=0;
samepoint=zeros(N_probe,10);%% the same point set
samepointsize=zeros(N_probe,1);%% the size of the i-th same point set
map_same=zeros(N_probe,1);
for i=1:N_probe
   if i==1
       nsamepoint=1;
       samepointsize(1)=1;
       samepoint(1,1)=i;
       map_same(i)=nsamepoint;
   else
       flag=1;
       for j=1:nsamepoint
          if norm(I_probe(samepoint(j,1),:)-I_probe(i,:))<epsilon
              samepointsize(j)=samepointsize(j)+1;
              samepoint(j,samepointsize(j))=i;
              map_same(i)=j;
              flag=0;
              break;
          end
       end
       
       if flag==1
          nsamepoint=nsamepoint+1;
          samepointsize(nsamepoint)=samepointsize(nsamepoint)+1;
          samepoint(nsamepoint,1)=i;
          map_same(i)=nsamepoint;
       end
   end
end

%%% loops,loopsize,nloops respectively denote the loops on the sphere, the size of each loop, the number of loops

loops=[];
loopsize=[];
nloops=0;
[loops_concave,loopsize_concave,nloops_concave] = loopconstruct_concave(segment,N_segment);
loops = [loops;loops_concave];
loopsize = [loopsize;loopsize_concave];
nloops = nloops+nloops_concave;

% compute the spherical patches
[patches,patchesize,npatches] = patchesconstruct_concave(k1,I,I_probe,segment,circle_interior,N_circle,loops_concave,loopsize_concave,nloops_concave,direct);

segment0 = zeros(N_segment,11);%%%[c,n,r,spoint,angle]
for i = 1:N_segment
    if direct == 1
        segment0(i,1:7) = segment(i,1:7);
        segment0(i,8:10) = I_probe(segment(i,9),:);
        segment0(i,11) = segment(i,8);
    else
        segment0(i,1:3) = segment(i,1:3);
        segment0(i,4:6) = -segment(i,4:6);
        segment0(i,7) = segment(i,7);
        segment0(i,8:10) = I_probe(segment(i,10),:);
        segment0(i,11) = segment(i,8);
    end
end

circle0 = circle_interior;
if direct == -1
    for i = 1:N_circle 
        circle0(i,4:6) = -circle_interior(i,4:6);
    end
    for i = 1:nloops
        loops(i,1:loopsize(i)) = flipud(loops(i,1:loopsize(i))')';
    end
end

global arg_eSAS;
global index_externalconcavepatch;
global index_completeconcavepatch;
global index_I_circle;
global index_I_segment;
global NB; % number of segments on the boundary of triangle-shaped patch

global Para;
global Figs;
global Ext;

%global V_eSES V_cSES;
for i = 1:npatches
    % mesh and visulize each concave patch

   index_I = k1;
   if arg_eSAS == 1
       index_externalconcavepatch = index_externalconcavepatch+1;
       index_concavepatch = index_externalconcavepatch;
   else
       index_completeconcavepatch = index_completeconcavepatch+1;
       index_concavepatch = index_completeconcavepatch;
   end
   
   if Para.out_MolStrc
        output_SES_concavepatches(atoms,NB,index_I_segment,index_I_circle,index_I,index_concavepatch,arg_eSAS,I(k1,:),Rp,loops,loopsize,segment0,circle0,patches(i,:),patchesize(i))
   end
   
   if arg_eSAS
       if Para.arg_meshing(1)
           if Para.arg_viz && Para.arg_ext
                figure(Figs.ext)
                hold on;
           end
           mesh_sphpat(I(k1,:),Rp,loops,loopsize,segment0,circle0,patches(i,:),patchesize(i));
       end
       CompAreaVol_concave(I(k1,:),Rp,loops,loopsize,segment0,circle0,patches(i,:),patchesize(i));
   end
   
   if arg_eSAS == 0
       if Ext.I(k1)==0 && Para.arg_meshing(1) 
           if Para.arg_viz && Para.arg_int
                figure(Figs.int)
                hold on;
           end
           mesh_sphpat(I(k1,:),Rp,loops,loopsize,segment0,circle0,patches(i,:),patchesize(i));
       end
       CompAreaVol_concave(I(k1,:),Rp,loops,loopsize,segment0,circle0,patches(i,:),patchesize(i));
   end
           
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Construct all spherical patches on the i-th sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [patches,patchesize,npatches]=patchesconstruct_concave(i,I,I_probe,segment,circle,N_circle,loops,loopsize,nloops,direct)
global Para;
Rp = Para.radius_probe;

patches=zeros(nloops+N_circle,nloops+N_circle); %% it records the boundary loops of spherical patches.
patchesize=zeros(nloops+N_circle,1);
npatches=0;

global tree;
% for j=1:(2*nloops+1)
%    tree(j)=treenode(0,[],0,0,0); 
% end
S0=[1:1:nloops,-1:-1:-N_circle];
tree(1)=treenode(1,S0,S0(1),nloops,N_circle);
ntree=1;
j=1;

for s=1:2*(nloops+N_circle)+1
    
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
        k=tree(j).activeelement;%%k denotes the fixed loop
        
        if length(S)==1
            tree(j).activenode=0;
            j=j+1;
            continue;
        end
        
        if k>0
            for s1=1:tree(j).n1
                point=I_probe(segment(loops(S(s1),1),9),:);%%point denotes the starting point of the loop
                if k==S(s1)||interiorloop_concave(point,I(i,:),Rp,loops(k,:),loopsize(k),I_probe,segment,direct)
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
                    [vector1,vector2]=orthogonalvectors(circle(-S(s2),4:6));
                    point=circle(-S(s2),1:3)+circle(-S(s2),7)*vector1;%%%[c,n,r]
                    
                    if interiorloop_concave(point,I(i,:),Rp,loops(k,:),loopsize(k),I_probe,segment,direct)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Construct all loops on the i-th sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [loops,loopsize,nloops]=loopconstruct_concave(segment,N_segment)
global map_same;
%%%segment(nsegment,1:5) is composed of [i,j,p1,p2,direction]
s = segment(:,9:10);

if N_segment ~= floor(N_segment)
   pause 
end

true = zeros(N_segment,1);%%records if the segment has been ensured to be in a loop
ntrue = 0;

loops = zeros(floor(N_segment/2),N_segment);
loopsize = zeros(floor(N_segment/2),1);
nloops = 0;%%%the number of loops

for k1 = 1:N_segment
    
   loop = zeros(1,N_segment);
   loopsi = 0;%% loopsize of the loop that we are handling with now.
   if true(k1) == 0 && ntrue < N_segment %%% ns means the number of segments on the ith atom
       nloops = nloops+1;
       pointer = k1;

       ntrue = ntrue+1; %%% ntrue of the already used segments.
       true(k1) = 1;       
       loopsi = loopsi+1;
       loop(loopsi) = pointer;
       
       times = 0;
       while ( map_same(s(k1,1)) ~= map_same(s(pointer,2)) )
           
           for k2 = k1+1:N_segment
               if true(k2) == 0 && map_same(s(k2,1)) == map_same(s(pointer,2))
                   pointer = k2;

                   loopsi = loopsi+1;
                   loop(loopsi) = pointer;
                   
                   true(pointer) = 1;                   
                   ntrue = ntrue+1;
               end
           end
           times = times+1;
           
           if times > 100
               disp('Fail to construct a loop!')
               pause
           end
       end
       
       loops(nloops,:) = loop;
       loopsize(nloops) = loopsi;
   end

end

loops = loops(1:nloops,:);
loopsize = loopsize(1:nloops);

end