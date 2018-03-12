%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meshsphere(c,r)
%figure('Position',[700 100 600 600])
global nsphere;
global T_sphere;
global Nt_sphere;
global P_sphere;
global Np_sphere;

nsphere = 0;

if nsphere == 0
    readsphere;
end
nsphere=nsphere+1;
%{
rid=[];
for i=1:Nt
    for j=i+1:Nt
       if (T(i,1)==T(j,2)&&T(i,2)==T(j,3)&&T(i,3)==T(j,1))||(T(i,1)==T(j,3)&&T(i,2)==T(j,1)&&T(i,3)==T(j,2))
           rid=[rid,j];
       end
    end
end
rid
T(rid,:)=[];
Nt=Nt-length(rid);
%}

P=P_sphere.*(r/1.2);
P=P+ones(Np_sphere,1)*c;
T=T_sphere;
Nt=Nt_sphere;
Np = Np_sphere;


figure
clf
axis off
zoom on
axis equal
view([145,10])%[138,37]
light('position',[1,1,1]);
lighting gouraud;
%}
k = 3;
for i = 1:k
    [T,Nt,P,Np] =remesh_sphere(c,r,T,Nt,P,Np);
end

P=real(P);
%line([P(E(1:Ne,1),1)';P(E(1:Ne,2),1)'],[P(E(1:Ne,1),2)';P(E(1:Ne,2),2)'],[P(E(1:Ne,1),3)';P(E(1:Ne,2),3)'],'LineStyle','-');%,'Marker','.'
X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];
% for i = 1:Nt_sphere
%     patch(X(:,i),Y(:,i),Z(:,i),'y','LineStyle','none','FaceAlpha',0.1);%'DiffuseStrength',0.8,
%     refined_triangle(c,r,X(:,i),Y(:,i),Z(:,i),1);
%     hold on;
% end


patch(X,Y,Z,'y','FaceAlpha',1);%'DiffuseStrength',0.8,
hold on;
%{
for i=1:Np
    str=strcat(num2str(i));
    %text(P(i,1),P(i,2),P(i,3),str);
end

hold on;
%}

global Triangle;
global Vertices;
global NormalVects;
Triangle = T;
Vertices = P;
for i = 1:Nt

    V1 = P(T(i,2),:)-P(T(i,1),:);
    V2 = P(T(i,3),:)-P(T(i,1),:);
    
    nv = cross(V1,V2);
    nv = nv/norm(nv);
    
    s = sign(nv*(P(T(i,1),:)-c)');
    
    if s == -1
        nv = -nv;
    end
    
    NormalVects(i,:) = nv;
end

writeSTL('sphere')
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

function refined_triangle(c,r,X,Y,Z,k)
% c
% c   refine the trianglar patch into 2^(2k) smaller triangles
% c

color = 'y';

n = 2^k;
A = [X(1),Y(1),Z(1)];
B = [X(2),Y(2),Z(2)];
C = [X(3),Y(3),Z(3)];

fa = 1;

for i = 1:n
    for j = 1:i
        if i < n
            lambda1 = [n-i,i-j+1,j-1]/n;
            lambda2 = [n-i,i-j,j]/n;
            lambda_up = [n-i+1,i-j,j-1]/n;
            lambda_down = [n-i-1,i-j+1,j]/n;
            A_up = lambda_up(1)*A+lambda_up(2)*B+lambda_up(3)*C;
            A_down = lambda_down(1)*A+lambda_down(2)*B+lambda_down(3)*C;
            B0 = lambda1(1)*A+lambda1(2)*B+lambda1(3)*C;
            C0 = lambda2(1)*A+lambda2(2)*B+lambda2(3)*C;
            
            A_up = (A_up-c)/norm(A_up-c,2)*r+c;
            A_down = (A_down-c)/norm(A_down-c,2)*r+c;
            B0 = (B0-c)/norm(B0-c,2)*r+c;
            C0 = (C0-c)/norm(C0-c,2)*r+c;
            
            x = [A_up(1);B0(1);C0(1)];
            y = [A_up(2);B0(2);C0(2)];
            z = [A_up(3);B0(3);C0(3)];
            patch(x,y,z,color,'FaceAlpha',fa);
            x = [A_down(1);B0(1);C0(1)];
            y = [A_down(2);B0(2);C0(2)];
            z = [A_down(3);B0(3);C0(3)];
            patch(x,y,z,color,'FaceAlpha',fa);
        else
            lambda1 = [n-i,i-j+1,j-1]/n;
            lambda2 = [n-i,i-j,j]/n;
            lambda_up = [n-i+1,i-j,j-1]/n;
            A_up = lambda_up(1)*A+lambda_up(2)*B+lambda_up(3)*C;
            B0 = lambda1(1)*A+lambda1(2)*B+lambda1(3)*C;
            C0 = lambda2(1)*A+lambda2(2)*B+lambda2(3)*C;
            
            A_up = (A_up-c)/norm(A_up-c,2)*r+c;
            B0 = (B0-c)/norm(B0-c,2)*r+c;
            C0 = (C0-c)/norm(C0-c,2)*r+c;
            
            x = [A_up(1);B0(1);C0(1)];
            y = [A_up(2);B0(2);C0(2)];
            z = [A_up(3);B0(3);C0(3)];
            patch(x,y,z,color,'FaceAlpha',fa);
        end
        
    end
end

end

function writeSTL(filename)

global Triangle;
global Vertices;
global NormalVects;

T = Triangle(:,1:3);
x = abs(floor((min(min(Vertices)))))+1;
V = Vertices(:,1:3)+x;
NV = NormalVects;

ID = ['./',filename,'.stl'];

N = size(T,1);

fileID = fopen(ID,'wt');
    fprintf(fileID,'%-16s \n',['solid ',filename]);
    fprintf(fileID,'\n');
    
    for i = 1:N
        v1 = V(T(i,1),:);
        v2 = V(T(i,2),:);
        v3 = V(T(i,3),:);
        
        fprintf(fileID,'%13s','facet normal '); fprintf(fileID,'%d %d %d \n',NV(i,:));
        fprintf(fileID,'%10s \n','outer loop');
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v1);
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v2);
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v3);
        fprintf(fileID,'%7s \n','endloop');
        fprintf(fileID,'%8s \n','endfacet');
        
        fprintf(fileID,'\n');
        
%         line(v1(1)+[0:100:1]*NV(i,1),v1(2)+[0:100:1]*NV(i,2),v1(3)+[0:100:1]*NV(i,3))
%         hold on;
    end
    
    fprintf(fileID,'%8s \n','endsolid');
fclose(fileID);

%{
for i = 1:N
    v1 = Vertices(T(i,1),1:3);
    v2 = Vertices(T(i,2),1:3);
    v3 = Vertices(T(i,3),1:3);

    line(v1(1)+[0:0.01:1]*NV(i,1),v1(2)+[0:0.01:1]*NV(i,2),v1(3)+[0:0.01:1]*NV(i,3))
    hold on;
end
%}
    
end

