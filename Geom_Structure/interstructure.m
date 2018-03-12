function intersection = interstructure(Geom,Rp)
% c
% c   create the intersection structure between SAS balls, with complexity of O(MlogM)
% c
% c   divide_centers: calculate the structure containing 0) max number of 
% c         centers in a small box 1) number of centers in a small box 2)
% c         indices of centers in a small box 3) collection of nonempty
% c         boxes 4) number of non empty boxes
% c
% c   M_int: intersection matrix recording the indice of VdW atoms which
% c         intersects the j-th atom with center in the box "ind"=(i,j)
% c   num_int: number of intersected atoms with j-th atom in the box "ind"
% c         =(i,j)
% c         
% c

rmax = max(Geom.R)+Rp; % the max radius of an SAS ball, 2.5 is the maximum UFF radius

% create a box enclosing all centers of atoms, composed of small boxes with side rmax
smallboxes = BOX(Geom,rmax); 

% divide centers into different small boxes
[nmax_centers,num_centers,ind_smallboxes,nonemptyboxes,num_nonemptyboxes] = divide_centers(Geom,smallboxes,rmax); 

% compute the intersectionship between VdW atoms
[M_int,num_int] = intersectionship(Geom,Rp,nmax_centers,num_centers,ind_smallboxes,nonemptyboxes,num_nonemptyboxes,smallboxes);

intersection.M_int = M_int; % n,i(n-th nonempty box,i-th center) --> [m_1,m_2,...,m_n] 
intersection.num_int = num_int;
% intersection.nonemptyboxes = nonemptyboxes; % n --> ind(indices of a nonempty box)
% intersection.ind_smallboxes = ind_smallboxes; % ind --> [m_1,m_2,...,m_n]
% intersection.num_centers = num_centers; % ind --> number of centers

end