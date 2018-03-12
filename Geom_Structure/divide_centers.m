function [nmax_centers,num_centers,ind_smallboxes,nonemptyboxes,num_nonemptyboxes] = divide_centers(Geom,smallboxes,rmax)
% c   
% c   divide all centers into different small boxes
% c
% c   num_centers: the number of centers in a small box
% c   ind_smallboxes: the indices of centers inside a small box
% c
% c   nonemptyboxes: collect all nonempty boxes each with index (i,j,k)
% c   num_nonemptyboxes: the number of nonempty boxes
% c

nmax_centers = Geom.M; % max number of possible centers in a small box
nx = smallboxes(1,2)-smallboxes(1,1); % number of boxes along with X axis
ny = smallboxes(2,2)-smallboxes(2,1); % number of boxes along with Y axis
nz = smallboxes(3,2)-smallboxes(3,1); % number of boxes along with Z axis
num_centers = zeros(nx,ny,nz);
ind_smallboxes = zeros(nx,ny,nz,nmax_centers);
for m = 1:Geom.M
    i0 = floor(Geom.centers(m,1)/rmax);
    j0 = floor(Geom.centers(m,2)/rmax);
    k0 = floor(Geom.centers(m,3)/rmax);
    i = i0-smallboxes(1,1)+1;
    j = j0-smallboxes(2,1)+1;
    k = k0-smallboxes(3,1)+1;
    
    num_centers(i,j,k) = num_centers(i,j,k)+1;
    num = num_centers(i,j,k);
    ind_smallboxes(i,j,k,num) = m;
end
nmax_centers = max(max(max(num_centers)));
ind_smallboxes = ind_smallboxes(:,:,:,1:nmax_centers);

% c
% c   collect the indices of nonempty boxes
% c
num_boxes = (smallboxes(1,2)-smallboxes(1,1))*(smallboxes(2,2)-smallboxes(2,1))*(smallboxes(3,2)-smallboxes(3,1));
num_nonemptyboxes = 0;
nonemptyboxes = zeros(num_boxes,3);
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            if num_centers(i,j,k)>0
                num_nonemptyboxes = num_nonemptyboxes+1;
                nonemptyboxes(num_nonemptyboxes,:) = [i,j,k];
            end
        end
    end
end

nonemptyboxes = nonemptyboxes(1:num_nonemptyboxes,:);

end