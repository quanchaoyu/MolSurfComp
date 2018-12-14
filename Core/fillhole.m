function fillhole
%
%   fill all possible inner holes with virtual atoms
%
DataGlob;

global neighbor;
global Nneighbor;

%% compute inner cavities (holes)
num_innerpatch = npatches-sum(ext_patch); % number of inner patches
max_innercavity = floor((num_innerpatch)/3); % max number of inner holes (cavities)

innerpatch = zeros(num_innerpatch,1);
inner = zeros(npatches,1); % if inner(i) is nonzero, it denotes the index of innercavity where innerpatch(i) lies
innercavity = zeros(max_innercavity,100); % innercavity(i,:) represents the set of patches consistuting the i-th inner cavity
size_innercavity = zeros(max_innercavity,1);
index_inner = 0; % index of inner cavities

num = 0;
for i = 1:npatches
    if ext_patch(i) == 0
        num = num+1;
        innerpatch(num) = i;
    end
end

for i = 1:num_innerpatch
    j = innerpatch(i); % the j-th patch
    if inner(j) == 0
        inner(j) = 1;
        
        % compute the num_inner inner cavity
        index_inner = index_inner+1;
        
        ninnerset = 1;
        innerset = j;
        for i0 = 1:num_innerpatch
            if i0 <= ninnerset
                j0 = innerset(i0);
                for k = 1:Nneighbor(j0)
                    n = neighbor(j0,k); % the n-th patch, which is neighboring to the j0-th path
                    if inner(n) == 0 
                        inner(n) = index_inner;
                        ninnerset = ninnerset+1;
                        innerset = [innerset,n];
                    end
                end
            else
                break;
            end
        end
        
        size_innercavity(index_inner) = ninnerset;
        innercavity(index_inner,1:ninnerset) = innerset;
    end
end

innercavity = innercavity(1:index_inner,1:max(size_innercavity));
size_innercavity = size_innercavity(1:index_inner,1);

%% fill all inner hols
inner_point = zeros(s,1);
inner_atom = zeros(M,1);
AddSpheres = [];
for i1 = 1:index_inner
    inner_atom = zeros(M,1);
    
    polyhedron = []; % polyhedron(i,;) = [index_I,ci,cj,cj], one triangle on the boundary of the polyhedron
    SASball_related = [];
    for j1 = 1:size_innercavity(i1,1) 
        j = innercavity(i1,j1); % j-th spherical patch
        
        for k=1:patchesize(j,1)
            if patches(j,k)>0 %  patches(j,k)>0 denotes the index of loop on the i-th sphere
                
               ln = loops_index(patch_atom(j),1)+patches(j,k)-1; %  the loop number
               for s0 = 1:loopsize(ln)
                   sn = loops(ln,s0); %  the segment number
                   
                   Ip1 = seg(sn,3); % intersection point
                   Ip2 = seg(sn,4);
                   
                   %%%%%%
                   
                   if inner_point(Ip1) == 0
                       inner_point(Ip1) = 1;
                       polyhedron = [polyhedron;Ip1,Iijk(Ip1,:)];
                       if inner_atom(Iijk(Ip1,1)) == 0
                           inner_atom(Iijk(Ip1,1)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip1,1)];
                       elseif inner_atom(Iijk(Ip1,2)) == 0
                           inner_atom(Iijk(Ip1,2)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip1,2)];
                       elseif inner_atom(Iijk(Ip1,3)) == 0
                           inner_atom(Iijk(Ip1,3)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip1,3)];
                       end
                       
                   end
                   if inner_point(Ip2) == 0
                       inner_point(Ip2) = 1;
                       polyhedron = [polyhedron;Ip2,Iijk(Ip2,:)];
                       if inner_atom(Iijk(Ip2,1)) == 0
                            inner_atom(Iijk(Ip2,1)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip2,1)];
                       elseif inner_atom(Iijk(Ip2,2)) == 0
                           inner_atom(Iijk(Ip2,2)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip2,2)];
                       elseif inner_atom(Iijk(Ip2,3)) == 0
                           inner_atom(Iijk(Ip2,3)) = 1;
                            SASball_related = [SASball_related,Iijk(Ip2,3)];
                       end
                   end
               end

            else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere

               cn = -circleindex(patch_atom(j),-patches(j,k)); %% the negative circle number
               if inner_atom(circle(-cn,1)) == 0
                   SASball_related = [SASball_related,circle(-cn,1)];
               elseif inner_atom(circle(-cn,2)) == 0
                   SASball_related = [SASball_related,circle(-cn,2)];
               end
               
            end
        end
    end
    add = add_spheres(SASball_related,polyhedron,innercavity(i1,1:size_innercavity(i1)),size_innercavity(i1),patch_atom,I,C,R,M,hightvalue,seg,ncrasegment,circle,ext_circle,ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index); % add(i,:) = [c,r]
    AddSpheres = [AddSpheres;add];
end

if size(AddSpheres,1)
    NewGeom(AddSpheres); % output new geometry with added spheres
end

global Para Figs;
if Para.viz_holes && size(AddSpheres,1)
    Figs.fillhole =figure('Name','Added spheres filling inner holes');
    figure(Figs.fillhole)
    for i = 1:size(AddSpheres,1)
        visusphere(AddSpheres(i,1:3),AddSpheres(i,4),0.2);
        hold on;
    end
    visu_innersegment;
    FigSetting;
end


end