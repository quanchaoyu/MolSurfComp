function add = add_spheres(SASball,polyhedron,cavity,npatch,patch_atom,I,C,R,M,hightvalue,segment,ncrasegment,circle,ext_circle,ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index)
%c
%c compute idealized SAS-spheres filling an inner cavity
%c 

global Para;
Rp = Para.radius_probe;

x_min = min(I(polyhedron(:,1),1));
x_max = max(I(polyhedron(:,1),1));
y_min = min(I(polyhedron(:,1),2));
y_max = max(I(polyhedron(:,1),2));
z_min = min(I(polyhedron(:,1),3));
z_max = max(I(polyhedron(:,1),3));
%delta = min(hightvalue(polyhedron(:,1)));

d_cavity = dimeter_cavity(I(polyhedron(:,1),:)); % dimeter of the inner cavity
d_cavity2surf = dist_cavity2surf(polyhedron(:,1),I); % dist from the cavity to other cavity or eSAS

parameter = 1.05;

add = [];

if d_cavity*parameter < d_cavity2surf/parameter
    delta = min(max(d_cavity2surf/parameter,0.4),2); % the radius of each added sphere
    add = [I(polyhedron(1,1),:),delta]; % in this case, one sphere is enough to fill the inner cavity
else
    delta = min(max(d_cavity2surf/(2*parameter),0.4),2);
    l = delta*2/sqrt(3)/parameter; % the length of a small cube
    nx = floor((x_max-x_min)/l)+1;
    ny = floor((y_max-y_min)/l)+1;
    nz = floor((z_max-z_min)/l)+1;
    
    for i = 0:nx
        for j = 0:ny
            for k = 0:nz
                x = x_min+i*l; % center (x,y,z), radius delta
                y = y_min+j*l;
                z = z_min+k*l;
                
                P = [x,y,z];
                
                %% Case 1: P is in the inner cavity
                CASE = 1;
                for t = 1:length(SASball)
                    ind = SASball(t);
                    if norm(P-C(ind,:)) < R(ind)+Rp % P is in a related SAS ball
                        CASE = 0; % it is not in Case 1
                        break;
                    end
                    if t == 1
                        dist0 = norm(P-C(ind,:))-(R(ind)+Rp);
                        ind0 = ind;
                    elseif dist0 > norm(P-C(ind,:))-(R(ind)+Rp)
                        dist0 = norm(P-C(ind,:))-(R(ind)+Rp);
                        ind0 = ind;
                    end
                end
                
                if length(SASball) == 0
                    
                end
                
                if CASE == 1
                    Xp=C(ind0,:)+(R(ind0)+Rp)*(P-C(ind0,:))/norm(P-C(ind0,:));
                    % check if the closest point Xp is on an inner patch
                    %plot3(P(1),P(2),P(3),'r.');
                    %plot3(Xp(1),Xp(2),Xp(3),'b.');
                    %hold on;
                    if interiorpatch_inner(Xp,ind0,C,R,I,Rp,segment,ncrasegment,ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index) == 0
                        CASE = 0;
                    end
                end

                
                if CASE == 1
                    add = [add;x,y,z,delta];
                    continue;
                end
                
                %% Case 2: P has a closest point on an inner patch, dist<delta
                % CASE = 0;
                for t = 1:length(SASball)
                    ind = SASball(t);
                    if norm(P-C(ind,:))<R(ind)+Rp % P is in the ind-ball
                        Xp=C(ind,:)+(R(ind)+Rp)*(P-C(ind,:))/norm(P-C(ind,:));
                        if norm(Xp-P) < delta
                            if interiorpatch_inner(Xp,ind,C,R,I,Rp,segment,ncrasegment,ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index) == 1
                                %plot3(P(1),P(2),P(3),'r.');
                                %plot3(Xp(1),Xp(2),Xp(3),'b.');
                                %hold on;
                                CASE = 2; % P has a closest point on an inner patch
                                break;
                            end
                        end
                    end
                end
                if CASE == 2
                    add = [add;x,y,z,delta];
                    continue;
                end
                
                %% Case 3: P has a closest point on an inner segment, dist<delta
                % CASE = 0;
                
                for t = 1:length(SASball)
                    ind = SASball(t);
                    if norm(P-C(ind,:))<R(ind)+Rp % P is in the ind-ball
                        index_s = patches_index(ind,1); % start and end indices of patches ont the i-th SAS-ball
                        index_e = patches_index(ind,2);
                        loops_i = loops(loops_index(ind,1):loops_index(ind,2),:);
                        loopsize_i = loopsize(loops_index(ind,1):loops_index(ind,2),:);

                        for patch_j = index_s:index_e
                            if ext_patch(patch_j) == 0 % the patch is an inner patch
                                for loop_k = 1:patchesize(patch_j,1)
                                    if patches(patch_j,loop_k) > 0
                                        for tt = 1:loopsize_i(patches(patch_j,loop_k),1)
                                            sn = loops_i(patches(patch_j,loop_k),tt); % segment index
                                            seg = [segment(sn,1:2),ncrasegment(sn,4:6),ncrasegment(sn,1:3),I(segment(sn,3),:),ncrasegment(sn,8),segment(sn,5)];
                                            if check_doublecone(P,seg,C,R,Rp)
                                                Xij = X(seg(1),seg(2),P,C,R,Rp);
                                                if norm(Xij-P) < delta
                                                    %plot3(P(1),P(2),P(3),'r.');
                                                    %hold on;
                                                    CASE = 3; % P has a closest point in an inner segment
                                                    break;
                                                end

                                            end
                                        end

                                    else
                                        cn = -patches(patch_j,loop_k); % circle index
                                        ind_i = circle(cn,1);
                                        ind_j = circle(cn,2);
                                        if ext_circle(cn) == 0 && check_Bij(P,ind_i,ind_j,C,R,Rp)
                                            Xij = X(ind_i,ind_j,P,C,R,Rp);
                                            if norm(Xij-P) < delta
                                                CASE = 3;
                                            end
                                        end
                                    end

                                    if CASE == 3
                                        break;
                                    end
                                end
                            end
                            if CASE == 3
                                break;
                            end
                        end

                    end
                end

                if CASE == 3
                    add = [add;x,y,z,delta];
                    continue;
                end
                
                %% Case 4: P has an inner intersection point as the closest point
                for s = 1:size(polyhedron,1)
                    if norm(P-I(polyhedron(s,1),:))<delta
                        %plot3(P(1),P(2),P(3),'r.');
                        %hold on;
                        add = [add;x,y,z,delta];
                    end
                end

            end
        end
    end
end

end

%% check if P is inside doublecone corresponding to seg
function true = check_doublecone(P,seg,C,R,Rp)
i = seg(1);
j = seg(2);
c = seg(3:5);
n = seg(6:8);
startpoint = seg(9:11);
angle = seg(12);
direct = seg(13);

Pi = C(i,:)+(P-C(i,:))*(R(i)+Rp)/(norm(C(i,:)-P));
Pj = C(j,:)+(P-C(j,:))*(R(j)+Rp)/(norm(C(j,:)-P));
if (R(i)+Rp)-norm(C(i,:)-Pj)>=0 && (R(j)+Rp)-norm(C(j,:)-Pi)>=0
    Xij = X(i,j,P,C,R,Rp);
    u = startpoint-c;
    v = Xij-c;
    if alpha(direct,u,v,n) <= angle
        true = 1;
    else
        true = 0;
    end
else
    true = 0;
end

end

function true = check_Bij(P,i,j,C,R,Rp)
Pi = C(i,:)+(P-C(i,:))*(R(i)+Rp)/(norm(C(i,:)-P));
Pj = C(j,:)+(P-C(j,:))*(R(j)+Rp)/(norm(C(j,:)-P));
if (R(i)+Rp)-norm(C(i,:)-Pj)>=0 && (R(j)+Rp)-norm(C(j,:)-Pi)>=0
    true = 1;
else
    true = 0;
end

end

function alpha=alpha(direct,u,v,n)
%c
%c  Compute the angle between two neighbor vectors u and v
%c
t=sign(det([u;v;n]));%direct=1 means clockwise
if direct*t>0
    alpha=acos(u*v'/(norm(u)*norm(v)));
else
    alpha=2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the nearest point of  P in Bij %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xp=X(i,j,P,C,R,Rp)
d=norm(C(i,:)-C(j,:));
t=((R(i)+Rp)^2-(R(j)+Rp)^2+d^2)/(2*d);
n=(C(j,:)-C(i,:))/norm(C(j,:)-C(i,:));
O=C(i,:)+t*n;
P1=P+((O-P)*n')*n;
Xp=O+sqrt((R(i)+Rp)^2-t^2)*(P1-O)/norm(P1-O);
end

%% compute the dist from the inner cavity to the eSAS or other inner cavities
function d = dist_cavity2surf(index,I)
global Ext; 
ext_I = Ext.I;
outcavity = ones(size(I,1),1);

for i = 1:size(index,1)
    outcavity(index(i)) = 0;
end

%first = 1;
ind_ext = find(ext_I == 1);


for i = 1:size(index,1)
    x = I(index(i),:);
    
    D = I - ones(size(I,1),1)*x;
    S = sum(D.^2,2);
    d = min(S(ind_ext))^0.5;
end
%     for j = 1:size(I,1)
%         if outcavity(j) == 1 && ext_I(j) == 1
%             y = I(j,:);
%             if first == 1
%                 d = norm(x-y);
%                 first = 0;
%             elseif norm(x-y) < d
%                 d = norm(x-y);
%             end
%         end
%     end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  check if Xp is on an inner patch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function true = interiorpatch_inner(Xp,i,C,R,I,Rp,segment,ncrasegment,ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index)

point = Xp;
ci = C(i,:);
ri = R(i)+Rp;

index_s = patches_index(i,1); % start and end indices of patches ont the i-th SAS-ball
index_e = patches_index(i,2);

loops_i = loops(loops_index(i,1):loops_index(i,2),:);
loopsize_i = loopsize(loops_index(i,1):loops_index(i,2),:);

true = 0;
for j = index_s:index_e
    if ext_patch(j) == 0 % the patch is an inner patch
        true0 = 1;
        for k = 1:patchesize(j,1)
            if patches(j,k) > 0
                if interiorloop(point,ci,C,ri,i,loops_i(patches(j,k),:),loopsize_i(patches(j,k),:),I,segment,ncrasegment) == 0 % Xp is outside a loop
                    true0 = 0;
                    break;
                end
            end
        end
        
        if true0 == 1 % there exists an inner patch containing Xp
            true = 1;
            break;
        end
    end
end

end

%% compute the dimeter of the inner cavity
function d = dimeter_cavity(inter) % inter is a n*3 matrix
d = 0;
for i = 1:size(inter,1)
    for j = i+1:size(inter,1)
        d0 = norm(inter(i,:)-inter(j,:));
        if d0 > d
            d = d0;
        end
    end
end
end
