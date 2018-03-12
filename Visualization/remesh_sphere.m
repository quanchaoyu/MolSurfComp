function [T_new,Nt_new,P,Np] = remesh_sphere(c,r,T,Nt,P,Np)

% c
% c   refine the trianglar patch into 2^(2k) smaller triangles
% c

% global T_sphere;
% global Nt_sphere;
% global P_sphere;
% global Np_sphere;
% 
% readsphere;
% 
% P = P_sphere.*(r/1.2);
% P = P+ones(Np_sphere,1)*c;
% T = T_sphere;
% Nt = Nt_sphere;
% Np = Np_sphere;


ind_mid = zeros(Np,Np);

for i = 1:Nt
    S = sort(T(i,:));
    
    midpoint = 0.5*(P(S(1),:)+P(S(2),:));
    midpoint = midpoint/norm(midpoint)*r;
    if ind_mid(S(1),S(2)) == 0
        ind_mid(S(1),S(2)) = Np+1;
        Np = Np+1;
        P(Np,:) = midpoint;
    end
    
    midpoint = 0.5*(P(S(2),:)+P(S(3),:));
    midpoint = midpoint/norm(midpoint)*r;
    if ind_mid(S(2),S(3)) == 0
        ind_mid(S(2),S(3)) = Np+1;
        Np = Np+1;
        P(Np,:) = midpoint;
    end
    
    midpoint = 0.5*(P(S(1),:)+P(S(3),:));
    midpoint = midpoint/norm(midpoint)*r;
    if ind_mid(S(1),S(3)) == 0
        ind_mid(S(1),S(3)) = Np+1;
        Np = Np+1;
        P(Np,:) = midpoint;
    end
    
    
end

T_new = zeros(4*Nt,3);
Nt_new = 4*Nt;
for i = 1:Nt
    S = sort(T(i,:));
    ind = [ind_mid(S(1),S(2)),ind_mid(S(2),S(3)),ind_mid(S(1),S(3))];
    T_new(4*i-3,:) = [S(1),ind(1),ind(3)];
    T_new(4*i-2,:) = [S(2),ind(1),ind(2)];
    T_new(4*i-1,:) = [S(3),ind(2),ind(3)];
    T_new(4*i,:) = ind;
end


end

function readsphere
global T_sphere;
global Nt_sphere;
global P_sphere;
global Np_sphere;

T_sphere=textread('T_sphere.txt');
Nt_sphere=textread('Nt_sphere.txt');
P_sphere=textread('Pses_sphere.txt');
Np_sphere=textread('Np_sphere.txt');

end