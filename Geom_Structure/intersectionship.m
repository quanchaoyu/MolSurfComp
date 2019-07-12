function [M_int,num_int] = intersectionship(Geom,Rp,nmax_centers,num_centers,ind_smallboxes,nonemptyboxes,num_nonemptyboxes,smallboxes)
% c
% c   compute the intersectionship between two VdW atoms
% c   M_int: intersection matrix recording the intersected atoms with each atom
% c   num_int: number of intersected atoms with the j-th atom in the i-th box
% c

nmax_int = 400;Geom.M; % max number of intersected SAS balls with an arbitrary SAS ball
% M_int = zeros(num_nonemptyboxes,nmax_centers,nmax_int);
% num_int = zeros(num_nonemptyboxes,nmax_centers);
M_int = zeros(Geom.M,nmax_int); %   matrix of intersection
num_int = zeros(Geom.M,1);

for i = 1:num_nonemptyboxes
    ind = nonemptyboxes(i,:); % index of the i-th nonempty box
    for j = 1:num_centers(ind(1),ind(2),ind(3))
        m = ind_smallboxes(ind(1),ind(2),ind(3),j);
        [INT,num_INT] = INT_SAS(j,ind,Geom,Rp,num_centers,ind_smallboxes,smallboxes,nmax_int);
        
        size(M_int)
        M_int(m,:) = INT;
        num_int(m,1) = num_INT;
%         M_int(i,j,:) = INT;
%         num_int(i,j) = num_INT;
    end
end

nmax_int = max(num_int);
M_int = M_int(:,1:nmax_int);

end

function [INT,num_INT] = INT_SAS(j,ind,Geom,Rp,num_centers,ind_smallboxes,smallboxes,nmax_int)
% c
% c   compute the intersected VdW atoms with j-th atom in box "ind"
% c
ind_xmin = max(1,ind(1)-2);
ind_ymin = max(1,ind(2)-2);
ind_zmin = max(1,ind(3)-2);

nx = smallboxes(1,2)-smallboxes(1,1); % number of boxes along with X axis
ny = smallboxes(2,2)-smallboxes(2,1); % number of boxes along with Y axis
nz = smallboxes(3,2)-smallboxes(3,1); % number of boxes along with Z axis

ind_xmax = min(nx,ind(1)+2);
ind_ymax = min(ny,ind(2)+2);
ind_zmax = min(nz,ind(3)+2);

m = ind_smallboxes(ind(1),ind(2),ind(3),j);
c = Geom.centers(m,:);
r = Geom.R(m);

INT = zeros(1,nmax_int);
num_INT = 0;
for i0 = ind_xmin:ind_xmax
    for j0 = ind_ymin:ind_ymax
        for k0 = ind_zmin:ind_zmax
            for n = 1:num_centers(i0,j0,k0)
                m0 = ind_smallboxes(i0,j0,k0,n);
                if r+Geom.R(m0)+2*Rp-norm(c-Geom.centers(m0,:)) > 0  && m0 ~= m % check if two SAS balls intersect
                    num_INT =  num_INT+1;
                    INT(num_INT) = m0;
                end
            end
        end
    end
end

end
