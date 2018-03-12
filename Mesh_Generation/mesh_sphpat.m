function mesh_sphpat(c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize)
% c   
% c   Use Advancing-front Method to Mesh a (Convex or Concave) Spherical Patch
% c
% c   c_sphere -- the sphere center, 
% c   r_sphere -- sphere radius
% c   loops, loopsize -- the loops forming the spherical patch
% c   patches, patchesize -- the patch composed of indicies (positive for a loop, negative for a circle) of loops and circles
% c
% c   segment -- records the infomation of all cirulcar arcs on the boundary of the spherical patch
% c              [c,n,r,spoint,angle] center, normal vector outwards the sphere, radius
% c              start point of the segment, angle of the segment
% c
% c   circle -- records the infomation of all circles on the boundary of the spherical patch
% c             [c,n,r],pointing outside the sphere
% c
global Vertices;
global Triangle;
global NormalVects;
        
global Para Figs;
Rp = Para.radius_probe;

d = Para.arg_meshing(2);
tolerance = 0.8*min(Para.arg_meshing(2),r_sphere);

global active;
global nactive;
nactive=patchesize-1;

N=floor(4*pi*r_sphere^2/(sqrt(3)*tolerance^2/4));%  N denotes the max number of triangles on the spherical patch
T=[];
Nt=0;
P=[];
Np=0;

if r_sphere~= 0
    
    k = patches(1);% the k-th loop or the -k-th circle
    if k>0
        [Ae,Nae,P,Np]=loop_division(c_sphere,r_sphere,loops(k,:),loopsize(k),segment,P,Np,d,Rp);
    else
        [Ae,Nae,P,Np]=circle_division(c_sphere,r_sphere,circle(-k,:),P,Np,d,Rp);
    end
    
    active0(1)=activefront(0,[],[]);
    for i=2:nactive
        active0(i)=activefront(0,[],[]);
    end
    active=active0;

    for i=1:nactive
        k=patches(i+1);
        if k>0
            [Ae0,Nae0,P,Np]=loop_division(c_sphere,r_sphere,loops(k,:),loopsize(k),segment,P,Np,d,Rp);
        else
            [Ae0,Nae0,P,Np]=circle_division(c_sphere,r_sphere,circle(-k,:),P,Np,d,Rp);
        end
        
        active(i)=activefront(0,Ae0,Nae0);
    end

    [T,Nt,~,~,P,Np] = advancing_front_approach(c_sphere,r_sphere,N,T,Nt,Ae,Nae,P,Np,d,tolerance);
    
    if Para.out_STL
        sv=size(Vertices,1);
        st=size(Triangle,1);
        Vertices=[Vertices(1:sv,:);P(1:Np,:),zeros(Np,1)];
        Triangle=[Triangle(1:st,:);T(1:Nt,:)+sv,ones(Nt,1)*508];
        
        if r_sphere == Rp
            NV = compute_NV(T,P,c_sphere,-1);
            NormalVects = [NormalVects(1:st,:);NV];
        else
            NV = compute_NV(T,P,c_sphere,1);
            NormalVects = [NormalVects(1:st,:);NV];
        end
    end

    if Nt > 0 % there might exist an isolated point as a patch
        if r_sphere == Rp
            if Para.arg_viz && Para.viz_SES(1)
                visuspherical1(c_sphere,T,Nt,P,Np); % visu concave SES patch
            end
            
        else
            
            if Para.arg_viz && Para.viz_SES(1)
                visuspherical2(c_sphere,r_sphere,T,Nt,P,Np); % visu convex SES patch
            end

            if Para.arg_viz && Para.viz_SAS(1)
                visuspherical3(c_sphere,r_sphere,T,Nt,P,Np); % visu convex SAS patch
            end
        end
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Divide the loop into several edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ae,Nae,P,Np]=loop_division(c_sphere,r_sphere,loop,loopsize,segment,P,Np,d,Rp)
%%%[c,n,r,spoint,epoint,angle]=segment,n outwards the sphere,angle is clockwise 
global rj;
global Rj;
Ae=[];
Nae=0;
Np0=Np;

for i=1:loopsize
   s=loop(i);%%%the segment number
   c=segment(s,1:3);
   n=segment(s,4:6);
   r=segment(s,7);
   P1=segment(s,8:10);
   P2=segment(loop(mod(i,loopsize)+1),8:10);
   angle=segment(s,11);

   if r_sphere==Rp
       rj=Rp;
   else
       rj=Rj(s);
   end
   
   flag=0;
   if r<Rp && r_sphere>Rp
       flag=1;
   end
   if r_sphere>Rp %%%% to mesh the spherical patch on the SES
        c=c_sphere+(c-c_sphere)*(r_sphere-Rp)/r_sphere;
        r=r*(r_sphere-Rp)/r_sphere;
        P1=c_sphere+(P1-c_sphere)*(r_sphere-Rp)/r_sphere;
   end
   %%%if r_sphere>Rp && r<Rp, then the cusp case, we need to modify
   %%%the arc_division, flag=1!!!
   
   [Ae,Nae,P,Np]=arc_division(c,r,P1,P2,angle,n,Ae,Nae,P,Np,r_sphere,flag,d,Rp);
     
end


if Nae==0
    return;
end

Ae(Nae,2)=Np0+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Divide the Arc into Several Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ae,Nae,P,Np]=arc_division(c,r,P1,P2,angle,n,Ae,Nae,P,Np,r_sphere,flag,d,Rp)
%%c is the center, r is the arc, P1 and P2 are the starting and ending points, angle is the angle from P1 to P2
%%the clockwise direction, n always points to the intersection point k1
global rj;

theta0=pi/3;%% define the maximum angle change
N_division=max(real(floor(angle/theta0)+1),real(floor(r*angle/d)+1));
if flag==1
    N_arc=real(floor(angle/theta0)+1);
    N_division=N_arc*real(floor(r*angle/d/N_arc)+1);
end


if angle<2*pi && r_sphere>Rp && r_sphere<rj+Rp
    r_new=r*(rj*r_sphere)/((r_sphere-Rp)*(rj+Rp));
    N_division=max(real(floor(angle/theta0)+1),real(floor(r_new*angle/d)+1));%%modify the division number of the arc
end


u=(P1-c)/r;
v=[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];

if N_division > 1000
    
end

points=zeros(N_division,3);
edges=zeros(N_division,2);
for j=0:N_division-1
    angle_j=j/N_division*angle;
    P_j=r*cos(angle_j)*u+r*sin(angle_j)*v+c;
    points(j+1,:)=P_j;
    
    %plot3(P_j(1),P_j(2),P_j(3),'b.')
    
    edges(j+1,:)=[Np+j+1,Np+j+2];
end

%P2 = r*cos(angle)*u+r*sin(angle)*v+c;
if N_division==1 && norm(P2-P1)<10^-10 && angle < 0.1
   return; 
end

P=[P;points];
Np=Np+N_division;

Ae=[Ae;edges];
Nae=Nae+N_division;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Divide the loop into several edgeds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ae,Nae,P,Np]=circle_division(c_sphere,r_sphere,circle,P,Np,d,Rp)
%%%[c,n,r,spoint,epoint,angle]=segment,n outwards the sphere,angle is clockwise 

Ae=[];
Nae=0;

c=circle(1:3);
n=circle(4:6);
r=circle(7);

flag=0;
if r<Rp && r_sphere>Rp
   flag=1;
end

if r_sphere>Rp %%%% to mesh the spherical patch on the SES
    c=c_sphere+(c-c_sphere)*(r_sphere-Rp)/r_sphere;
    r=r*(r_sphere-Rp)/r_sphere;
end

[v1,v2]=orthogonalvectors(n);%%n pointing outhside
P1=c+r*v1;

Np0=Np;
[Ae,Nae,P,Np]=arc_division(c,r,P1,P1,2*pi,n,Ae,Nae,P,Np,r_sphere,flag,d,Rp);


%modify the active edge set
Ae(Nae,2)=Np0+1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%visu a spherical patch%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visuspherical1(c_sphere,T,Nt,P,Np)
% c
% c   viz the Concave SES spherical patch
% c
global Para;

T=real(T);
P=real(P);

% global Vertices;
% global Triangle;
% sv=size(Vertices,1);
% st=size(Triangle,1);
% Vertices=[Vertices(1:sv,:);P(1:Np,:)];
% Triangle=[Triangle(1:st,:);T(1:Nt,:)+sv];

X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];

patch(X,Y,Z,'b','LineStyle','none','DiffuseStrength',0.8,'FaceAlpha',Para.viz_SES(2),'FaceLighting','gouraud')%,'LineStyle','none','BackFaceLighting','lit');%
hold on;
end

function visuspherical2(c_sphere,r_sphere,T,Nt,P,Np)
% c
% c   viz the Convex SES spherical patch
% c
global Para;

global num_triangles;
num_triangles = num_triangles+Nt;


T=real(T);
P=real(P);

% global Vertices;
% global Triangle;
% sv=size(Vertices,1);
% st=size(Triangle,1);
% Vertices=[Vertices(1:sv,:);P(1:Np,:)];
% Triangle=[Triangle(1:st,:);T(1:Nt,:)+sv];

X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];

patch(X,Y,Z,'r','LineStyle','none','DiffuseStrength',0.8,'FaceAlpha',Para.viz_SES(2),'FaceLighting','gouraud');%
hold on;

end

function visuspherical3(c_sphere,r_sphere,T,Nt,P,Np)
% c
% c   viz the SAS spherical patches
% c
global Para;
Rp = Para.radius_probe;

P=ones(Np,1)*c_sphere+(P-ones(Np,1)*c_sphere)*((r_sphere)/(r_sphere-Rp));
X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];
patch(X,Y,Z,'r','LineStyle','none','DiffuseStrength',0.8,'FaceAlpha',Para.viz_SAS(2));%,'LineStyle','none'
hold on;

area=0;
for i=1:Nt
    d1=norm(P(T(i,1),:)-P(T(i,2),:));
    d2=norm(P(T(i,2),:)-P(T(i,3),:));
    d3=norm(P(T(i,3),:)-P(T(i,1),:));
    p=(d1+d2+d3)/2;
    add=sqrt(p*(p-d1)*(p-d2)*(p-d3));
    area=area+add;
end

global P_sas;
global T_sas;
sv=size(P_sas,1);
st=size(T_sas,1);
P_sas=[P_sas(1:sv,:);P];
T_sas=[T_sas(1:st,:);T+sv];
global ind_sas;
ind_sas = [ind_sas,st+Nt];

end
