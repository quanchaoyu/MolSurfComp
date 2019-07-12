%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute f(p) at point P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fsas,Xp,fses,patchtype,grad_fsas,laplace_fsas] = fsas(P,index_ball,...
    C,R,M,probe,I,inter,Row,Iijk,Ii,In,Csas_I,segment,ncrasegment,circle,...
    circleindex,ncircleindex,Csas_circle,Csas_patch,patches,patchesize,...
    patches_index,loops,loopsize,loops_index)
%c 
%c this function run N times (N is the number of grid points)
%c
%c compute the signed distance value from p to eSAS
%c Rk: P lies in the label-Ball
%c

S=Sp(P,C,R,probe,inter,Row,index_ball); % S consists of Bi containig P

CASE = 0;
%% CASE 1: P is out of the eSAS
if index_ball == 0 % P is out of every SAS-ball
    [fses,ind]=min((((ones(M,1)*P-C).^2)*ones(3,1)).^(1/2)-R');
    Xp=C(ind,:)+(R(ind)+probe)*(P-C(ind,:))/norm(P-C(ind,:)); % Xp is on the ball ind
    
    if interiorpatch_inner(Xp,ind,C,R,I,probe,segment,ncrasegment,...
            Csas_patch,patches,patchesize,patches_index,loops,loopsize,...
            loops_index,circle,circleindex,ncircleindex) == 0
        patchtype=[1,ind,0,0];
        fsas = fses-probe;
        grad_fsas = (P-C(ind,:))/norm(P-C(ind,:));
        laplace_fsas = 2/norm(P-C(ind,:));
        return;
    else
        CASE = 4;
    end
end

%% CASE 2: P is inside the eSAS and has a closest point on an external SAS patch
if CASE ~= 4
    for i=1:length(S)
        Pi=C(S(i),:)+(R(S(i))+probe)*(P-C(S(i),:))/norm(P-C(S(i),:));
        Xp = Pi;
        ind = S(i);
        if interiorpatch_external(Xp,ind,C,R,I,probe,segment,ncrasegment,...
                Csas_patch,patches,patchesize,patches_index,loops,loopsize,...
                loops_index,circle,circleindex,ncircleindex) == 1
            fses = -norm(P-Pi)+probe;
            fsas = fses-probe;
            patchtype=[2,S(i),0,0];
            grad_fsas = -(C(ind,:)-P)/norm(C(ind,:)-P);%(C(ind,:)-P)/norm(C(ind,:)-P);%????
            laplace_fsas = 2/norm(P-C(ind,:));
            return;
        end
    end
end

%% CASE 3: P is inside the eSAS and has a closest point on an external SAS segment
if CASE ~= 4
    for i = 1:length(S)
        ind = S(i);
        % start and end indices of patches ont the i-th SAS-ball
        index_s = patches_index(ind,1); 
        index_e = patches_index(ind,2);
        if index_s == 0 % the ind-th SAS-sphere is completely covered by others
            continue;
        end

        if loops_index(ind,1) ~= 0 % there exists a loop (not circle) on the sphere
            loops_i = loops(loops_index(ind,1):loops_index(ind,2),:);
            loopsize_i = loopsize(loops_index(ind,1):loops_index(ind,2),:);
        else
            loops_i = [];
            loopsize_i = 0;
        end
        %loops_i = loops(loops_index(ind,1):loops_index(ind,2),:);
        %loopsize_i = loopsize(loops_index(ind,1):loops_index(ind,2),:);

        for patch_j = index_s:index_e
            if Csas_patch(patch_j) == 1 % the patch is an external patch
                for loop_k = 1:patchesize(patch_j,1)
                    if patches(patch_j,loop_k) > 0
                        %loopsize_i
                        for tt = 1:loopsize_i(patches(patch_j,loop_k),1)
                            sn = loops_i(patches(patch_j,loop_k),tt); % segment index
                            seg = [segment(sn,1:2),ncrasegment(sn,4:6),...
                                ncrasegment(sn,1:3),I(segment(sn,3),:),...
                                ncrasegment(sn,8),segment(sn,5)];
                            if check_doublecone(P,seg,C,R,probe)
                                [Xij,scale] = X(seg(1),seg(2),P,C,R,probe);
                                Xp = Xij;
                                fsas = -norm(P-Xp);
                                fses = fsas+probe;
                                patchtype = [3,seg(1:2),0];
                                grad_fsas = (Xp-P)/norm(Xp-P);
                                laplace_fsas = (-2+scale)/norm(Xp-P);
                                return;
                            end
                        end
                    else
                         % circle index
                        cn = circleindex(ind,-patches(patch_j,loop_k));
                        ind_i = circle(cn,1);
                        ind_j = circle(cn,2);
                        if Csas_circle(cn) == 1 && check_Bij(P,ind_i,ind_j,C,R,probe)
                            [Xij,scale] = X(ind_i,ind_j,P,C,R,probe);
                            Xp = Xij;
                            fsas = -norm(P-Xp);
                            fses = fsas+probe;
                            patchtype = [3,ind_i,ind_j,0];
                            grad_fsas = (Xp-P)/norm(Xp-P);
                            laplace_fsas = (-2+scale)/norm(Xp-P);
                            return;
                        end
                    end

                end
            end
        end
    end
end

%% CASE 4: P has an external SAS intersection point as the closest point
%Case I: P is inside a tetrahedron corresponding to an intersection point
for i=1:length(S)
    for j=1:In(S(i))
        k=Ii(S(i),j); % index of intersection point
        if Csas_I(k) == 1 % external intersection point
            Ip=I(k,:);
            i0=Iijk(k,1);j0=Iijk(k,2);k0=Iijk(k,3);
            if tetrahedron([P,1],[C(i0,:),1],[C(j0,:),1],[C(k0,:),1],...
                    [Ip,1]) % P belongs to the tetrahedron corresponding to I(k,:)
                Xp=Ip;
                fses=-norm(Ip-P)+probe;
                fsas = fses-probe;
                patchtype=[4,i0,j0,k0];
                grad_fsas = (Xp-P)/norm(Xp-P);
                laplace_fsas = -2/norm(Xp-P);
                return;
            end
        end
    end
end

%Case II: P is outside each tetrahedron corresponding to an intersection point
I_external = I(find(Csas_I)',:);
s_external=sum(Csas_I); % number of external intersection points

if isempty(I_external) == 0
    [fses,j]=max(-(((ones(s_external,1)*P-I_external).^2)*ones(3,1)).^(1/2)+probe);
    Xp=I_external(j,:);
    fsas = fses-probe;
    patchtype=[4.1,Iijk(j,1),Iijk(j,2),Iijk(j,3)];
    grad_fsas = (Xp-P)/norm(Xp-P);
    laplace_fsas = -2/norm(Xp-P);
else
    [fses,ind]=min((((ones(M,1)*P-C).^2)*ones(3,1)).^(1/2)-R');
    Xp=C(ind,:)+(R(ind)+probe)*(P-C(ind,:))/norm(P-C(ind,:)); % Xp is on the ball ind
    fsas = fses-probe;
    patchtype=[2,ind,0,0];
    grad_fsas = -(C(ind,:)-P)/norm(C(ind,:)-P);% ??????
    laplace_fsas = 2/norm(P-C(ind,:));
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if P is inside double-cone region corresponding to seg or circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function true = check_doublecone(P,seg,C,R,probe)
i = seg(1);
j = seg(2);
c = seg(3:5);
n = seg(6:8);
startpoint = seg(9:11);
angle = seg(12);
direct = seg(13);

Pi = C(i,:)+(P-C(i,:))*(R(i)+probe)/(norm(C(i,:)-P));
Pj = C(j,:)+(P-C(j,:))*(R(j)+probe)/(norm(C(j,:)-P));
if (R(i)+probe)-norm(C(i,:)-Pj)>=0 && (R(j)+probe)-norm(C(j,:)-Pi)>=0
    Xij = X(i,j,P,C,R,probe);
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

function true = check_Bij(P,i,j,C,R,probe)
Pi = C(i,:)+(P-C(i,:))*(R(i)+probe)/(norm(C(i,:)-P));
Pj = C(j,:)+(P-C(j,:))*(R(j)+probe)/(norm(C(j,:)-P));
if (R(i)+probe)-norm(C(i,:)-Pj)>=0 && (R(j)+probe)-norm(C(j,:)-Pi)>=0
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
%%  check if Xp is on a spherical patch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function true = interiorpatch_inner(Xp,i,C,R,I,probe,segment,ncrasegment,...
    Csas_patch,patches,patchesize,patches_index,loops,loopsize,loops_index,...
    circle,circleindex,ncircleindex)

point = Xp;
ci = C(i,:);
ri = R(i)+probe;

index_s = patches_index(i,1); % start and end indices of patches ont the i-th SAS-ball
index_e = patches_index(i,2);

if index_s == 0 %% the i-th SAS-sphere is completely covered by others
    true = 0;
    return;
end

if loops_index(i,1) ~= 0 % there exists a loop (not circle) on the sphere
    loops_i = loops(loops_index(i,1):loops_index(i,2),:);
    loopsize_i = loopsize(loops_index(i,1):loops_index(i,2),:);
else
    loops_i = [];
    loopsize_i = 0;
end

% for j = 1:ncircleindex(i)
%     ind_circle = circleindex(i,j); % index of circle
%     
%     ind1 = circle(ind_circle,1);
%     ind2 = circle(ind_circle,2);
%     A = circle(ind_circle,3:5); % center
%     n = circle(ind_circle,6:8); % normal vector from ind1 to ind2
%     
%     t = (ind1 == i);
%     if t*((Xp-A)*n') > 0 % Xp is eaten
%         true = 0;
%         return;
%     end
%     
% end

% if loops_index(i,1) == 0
%     true = 0;
%     return;
% end
% 
% loops_i = loops(loops_index(i,1):loops_index(i,2),:);
% loopsize_i = loopsize(loops_index(i,1):loops_index(i,2),:);

true = 0;
for j = index_s:index_e
    if Csas_patch(j) == 0 % the patch is an inner patch
        true0 = 1;
        for k = 1:patchesize(j,1)
            if patches(j,k) > 0
                if interiorloop_dist(point,ci,C,ri,i,R,proble,loops_i(patches(j,k),:),...
                        loopsize_i(patches(j,k),:),I,segment,ncrasegment) == 0 
                    % Xp is outside a loop
                    true0 = 0;
                    break;
                end
            else
                ind_circle = circleindex(i,-patches(j,k));

                ind1 = circle(ind_circle,1);
                ind2 = circle(ind_circle,2);
                A = circle(ind_circle,3:5); % center
                n = circle(ind_circle,6:8); % normal vector from ind1 to ind2

                t = 2*(ind1 == i)-1;
                if t*((Xp-A)*n') > 0 % Xp is eaten
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

function true = interiorpatch_external(Xp,i,C,R,I,probe,segment,ncrasegment,...
    Csas_patch,patches,patchesize,patches_index,loops,loopsize,loops_index,...
    circle,circleindex,ncircleindex)

point = Xp;
ci = C(i,:);
ri = R(i)+probe;

index_s = patches_index(i,1); % start and end indices of patches on the i-th SAS-ball
index_e = patches_index(i,2);

if index_s == 0 %% the i-th SAS-sphere is completely covered by others
    true = 0;
    return;
end

if loops_index(i,1) ~= 0 % there exists a loop (not circle) on the sphere
    loops_i = loops(loops_index(i,1):loops_index(i,2),:);
    loopsize_i = loopsize(loops_index(i,1):loops_index(i,2),:);
else
    loops_i = [];
    loopsize_i = 0;
end

% for j = 1:ncircleindex(i)
%     ind_circle = circleindex(i,j); % index of circle
%     
%     ind1 = circle(ind_circle,1);
%     ind2 = circle(ind_circle,2);
%     A = circle(ind_circle,3:5); % center
%     n = circle(ind_circle,6:8); % normal vector from ind1 to ind2
%     
%     t = (ind1 == i);
%     if t*((Xp-A)*n') > 0 % Xp is eaten
%         true = 0;
%         return;
%     end
%     
% end

true = 0;
for j = index_s:index_e
    if Csas_patch(j) == 1 % the patch is an external patch
        
        true0 = 1;
        for k = 1:patchesize(j,1)
            if patches(j,k) > 0
                if interiorloop_dist(point,ci,C,ri,i,R,probe,loops_i(patches(j,k),:),...
                        loopsize_i(patches(j,k),:),I,segment,ncrasegment) == 0 
                    % Xp is outside a loop
                    true0 = 0;
                    break;
                end
            else
                ind_circle = circleindex(i,-patches(j,k));
                 
                ind1 = circle(ind_circle,1);
                ind2 = circle(ind_circle,2);
                A = circle(ind_circle,3:5); % center
                n = circle(ind_circle,6:8); % normal vector from ind1 to ind2

                t = 2*(ind1 == i)-1;
                if t*((Xp-A)*n') > 0 % Xp is eaten, i.e. Xp is outside a circle
                    true0 = 0;
                    break;
                end
            end
        end
        
        if true0 == 1 % there exists an external patch containing Xp
            true = 1;
            break;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  return all indices i s.t. Bi contains P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = Sp(P,C,R,probe,inter,Row,index_ball)
S=zeros(100);%kmax
t=1;
S(t)=index_ball;
if index_ball == 0
    S = [];
    return;
end
if norm(P-C(index_ball,:))>R(index_ball)+probe
    fprintf('error');%%error case
end

for row=1:Row(index_ball)
    j=inter(index_ball,row);
    if norm(P-C(j,:))<(R(j)+probe)
        t=t+1;
        S(t)=j;
    end
end

S=S(1:t);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  return all the indices (i,j) s.t. Bij contains P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=Kp(P,C,R,probe,inter,Row,S)
K=zeros(length(S)*100,2);%%kmax=100
t=0;

for i=1:length(S)
    Pi=C(S(i),:)+(P-C(S(i),:))*(R(S(i))+probe)/(norm(C(S(i),:)-P));
    %Consider the situation P belongs to the intersection of two balls
    for j=(i+1):length(S)
        Pj=C(S(j),:)+(P-C(S(j),:))*(R(S(j))+probe)/(norm(C(S(j),:)-P));
        if (R(S(i))+probe)-norm(C(S(i),:)-Pj)>=0&&...
                (R(S(j))+probe)-norm(C(S(j),:)-Pi)>=0
            t=t+1;
            K(t,:)=[S(i),S(j)];
        end
    end
    %Consider the situation P only belongs to one ball S(i)
    for row=1:Row(S(i))
        j=inter(S(i),row);
        if  ~any(S==j)&&(R(j)+probe)-norm(C(j,:)-Pi)>=0 %%triangle section
            t=t+1;
            K(t,:)=[S(i),j];
        end
    end
end

K=K(1:t,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the nearest point of  P in Bij %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xp,scale]=X(i,j,P,C,R,probe)
d=norm(C(i,:)-C(j,:));
t=((R(i)+probe)^2-(R(j)+probe)^2+d^2)/(2*d);
n=(C(j,:)-C(i,:))/norm(C(j,:)-C(i,:));
O=C(i,:)+t*n;
P1=P+((O-P)*n')*n;

if norm(P1-O) ~= 0
    Xp=O+sqrt((R(i)+probe)^2-t^2)*(P1-O)/norm(P1-O);
    scale = sqrt((R(i)+probe)^2-t^2)/norm(P1-O);
else
    [v,~] = orthogonalvectors(n);
    Xp=O+sqrt((R(i)+probe)^2-t^2)*v;
    scale = sqrt((R(i)+probe)^2-t^2)*0;%10^8;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if P is in the tetrahedron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t=tetrahedron(P,X1,X2,X3,X4)
D0=sign(det([X1;X2;X3;X4]));
D1=sign(det([P;X2;X3;X4]));
D2=sign(det([X1;P;X3;X4]));
D3=sign(det([X1;X2;P;X4]));
D4=sign(det([X1;X2;X3;P]));
t=0;
if D0==D1&&D0==D2&&D0==D3&&D0==D4
    t=1;
end
end
