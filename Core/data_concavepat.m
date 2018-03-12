function data_concavepat(i,hight_set,K,Kn,I,Iijk,C,probe,direction)
%c
%c compute the data structure of a concave spherical patch, including segments and circles
%c
epsilon=10^(-10);

i0=hight_set(i);%% i0-th intersection point
atoms = Iijk(i0,:);
x=I(i0,:);
ci=C(Iijk(i0,1),:);
cj=C(Iijk(i0,2),:);
ck=C(Iijk(i0,3),:);
ni=(ci-x)/norm(ci-x);
nj=(cj-x)/norm(cj-x);
nk=(ck-x)/norm(ck-x);
nij=[ni(2)*nj(3)-ni(3)*nj(2),ni(3)*nj(1)-ni(1)*nj(3),ni(1)*nj(2)-ni(2)*nj(1)];
njk=[nj(2)*nk(3)-nj(3)*nk(2),nj(3)*nk(1)-nj(1)*nk(3),nj(1)*nk(2)-nj(2)*nk(1)];
nki=[nk(2)*ni(3)-nk(3)*ni(2),nk(3)*ni(1)-nk(1)*ni(3),nk(1)*ni(2)-nk(2)*ni(1)];
nij=nij/norm(nij);
njk=njk/norm(njk);
nki=nki/norm(nki);
xi=x+probe*ni;%% xi,xj,xk denote the three touching vertices
xj=x+probe*nj;
xk=x+probe*nk;

%%%%to compute the I_boundary
circle_ij=zeros(2*Kn,5);
ncircle_ij=0;
circle_jk=zeros(2*Kn,5);
ncircle_jk=0;
circle_ki=zeros(2*Kn,5);
ncircle_ki=0;

%%%%%%%%% modify K and Kn by removing some unuseful nearby probe
Kn0 = 0;
K0 = zeros(Kn,1);

remove = zeros(Kn,1);
Num_ij = 0;
Num_jk = 0;
Num_ki = 0;
Kij = 0; % Kij records the useful probe intersection cicj
Kjk = 0;
Kki = 0;
for j = 1:Kn
    j0 = hight_set(K(j));
    x0 = I(j0,:);
    if any(Iijk(i0,1)-Iijk(j0,:)==0) == 1 && any(Iijk(i0,2)-Iijk(j0,:)==0) == 1
         remove(j) = 1;
         Num_ij = Num_ij+1;
         if Num_ij == 1
             Kij = j;
             u = -ni-(-ni*(cj-ci)'/norm(cj-ci))*(cj-ci)/norm(cj-ci);
             v = (x0-ci)-((x0-ci)*(cj-ci)'/norm(cj-ci))*(cj-ci)/norm(cj-ci);
             alphaij = alpha(u,v,cj-ci);
             if direction(i0,1) == -1
                 alphaij = 2*pi - alphaij;
             end
         else
             u = -ni-(-ni*(cj-ci)'/norm(cj-ci))*(cj-ci)/norm(cj-ci);
             v = (x0-ci)-((x0-ci)*(cj-ci)'/norm(cj-ci))*(cj-ci)/norm(cj-ci);
             alpha0 = alpha(u,v,cj-ci);
             if direction(i0,1) == -1
                 alpha0 = 2*pi - alpha0;
             end
             
             if alpha0 > alphaij
                Kij = j;
                alphaij = alpha0;
             end
         end
    end
    
    if any(Iijk(i0,2)-Iijk(j0,:)==0) == 1 && any(Iijk(i0,3)-Iijk(j0,:)==0) == 1
         remove(j) = 1;
         Num_jk = Num_jk+1;
         if Num_jk == 1
             Kjk = j;
             u = -nj-(-nj*(ck-cj)'/norm(ck-cj))*(ck-cj)/norm(ck-cj);
             v = (x0-cj)-((x0-cj)*(ck-cj)'/norm(ck-cj))*(ck-cj)/norm(ck-cj);
             alphajk = alpha(u,v,ck-cj);
             if direction(i0,1) == -1
                 alphajk = 2*pi - alphajk;
             end
         else
             u = -nj-(-nj*(ck-cj)'/norm(ck-cj))*(ck-cj)/norm(ck-cj);
             v = (x0-cj)-((x0-cj)*(ck-cj)'/norm(ck-cj))*(ck-cj)/norm(ck-cj);
             alpha0 = alpha(u,v,ck-cj);
             if direction(i0,1) == -1
                 alpha0 = 2*pi - alpha0;
             end
             if alpha0 > alphajk
                Kjk = j;
                alphajk = alpha0;
             end
         end
    end
    
    if any(Iijk(i0,3)-Iijk(j0,:)==0) == 1 && any(Iijk(i0,1)-Iijk(j0,:)==0) == 1
         remove(j) = 1;
         Num_ki = Num_ki+1;
         if Num_ki == 1
             Kki = j;
             u = -nk-(-nk*(ci-ck)'/norm(ci-ck))*(ci-ck)/norm(ci-ck);
             v = (x0-ck)-((x0-ck)*(ci-ck)'/norm(ci-ck))*(ci-ck)/norm(ci-ck);
             alphaki = alpha(u,v,ci-ck);
             if direction(i0,1) == -1
                 alphaki = 2*pi - alphaki;
             end
         else
             u = -nk-(-nk*(ci-ck)'/norm(ci-ck))*(ci-ck)/norm(ci-ck);
             v = (x0-ck)-((x0-ck)*(ci-ck)'/norm(ci-ck))*(ci-ck)/norm(ci-ck);
             alpha0 = alpha(u,v,ci-ck);
             if direction(i0,1) == -1
                 alpha0 = 2*pi - alpha0;
             end
             if alpha0 > alphaki
                Kki = j;
                alphaki = alpha0;
             end
         end
    end
end

if Kij>0
    remove(Kij) = 0;
end
if Kjk>0
    remove(Kjk) = 0;
end
if Kki>0
    remove(Kki) = 0;
end

for j = 1:Kn
    if remove(j)==0
        Kn0 = Kn0+1;
        K0(Kn0) = K(j);
    end
end
K0 = K0(1:Kn0);
K = K0;
Kn = Kn0;



%%% to compute out all the circles on the three corresponding
%%% planes:x-xi-xj,x-xj-xk,x-xk-xi
for j=1:Kn
    if K(j)<1
        
    end
    j0=hight_set(K(j));%%j0-th intersection point
    x0=I(j0,:);%%center of the probe
    
    tij=(x-x0)*nij';
    if abs(tij)<probe
        cij=x0+tij*nij;%%center of the intersection circle corresponding to ij-plane
        rij=sqrt(probe^2-tij^2);
        ncircle_ij=ncircle_ij+1;
        circle_ij(ncircle_ij,:)=[cij,rij,j];
    end
   
    tjk=(x-x0)*njk';
    if abs(tjk)<probe
        cjk=x0+tjk*njk;%%center of the intersection circle corresponding to jk-plane
        rjk=sqrt(probe^2-tjk^2);
        ncircle_jk=ncircle_jk+1;
        circle_jk(ncircle_jk,:)=[cjk,rjk,j];
    end
    
    tki=(x-x0)*nki';
    if abs(tki)<probe
        cki=x0+tki*nki;%%center of the intersection circle corresponding to ki-plane
        rki=sqrt(probe^2-tki^2);
        ncircle_ki=ncircle_ki+1;
        circle_ki(ncircle_ki,:)=[cki,rki,j];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Compute the intersection point on the boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boundary_ij=zeros(2*Kn,4);%%Kn is the number of corresponding probes
boundary_jk=zeros(2*Kn,4);
boundary_ki=zeros(2*Kn,4);
N_ij=2;N_jk=2;N_ki=2;
boundary_ij(1:2,:)=[xi,0;xj,0];
boundary_jk(1:2,:)=[xj,0;xk,0];
boundary_ki(1:2,:)=[xk,0;xi,0];


for k=1:ncircle_ij
    cij=circle_ij(k,1:3);
    rij=circle_ij(k,4);
    j=circle_ij(k,5);
    if norm(cij-x)<probe+rij && norm(cij-x)>probe-rij
        [x1,x2]=intersection_circle(x,probe,cij,rij,nij);
        if test1(x1,x,xi,xj,xk,k,circle_ij,ncircle_ij) && min(sum(abs(boundary_ij(1:N_ij,1:3)-ones(N_ij,1)*x1),2))>epsilon
            N_ij=N_ij+1;
            boundary_ij(N_ij,:)=[x1,j];
            
        end
        if test1(x2,x,xi,xj,xk,k,circle_ij,ncircle_ij) && min(sum(abs(boundary_ij(1:N_ij,1:3)-ones(N_ij,1)*x2),2))>epsilon
            N_ij=N_ij+1;
            boundary_ij(N_ij,:)=[x2,j];
        end
    end
end

for k = 1:ncircle_jk
    cjk = circle_jk(k,1:3);
    rjk = circle_jk(k,4);
    j = circle_jk(k,5);
    if norm(cjk-x) < probe+rjk && norm(cjk-x) > probe-rjk
        [x1,x2] = intersection_circle(x,probe,cjk,rjk,njk);
        if test1(x1,x,xj,xk,xi,k,circle_jk,ncircle_jk) && min(sum(abs(boundary_jk(1:N_jk,1:3)-ones(N_jk,1)*x1),2))>epsilon
            N_jk = N_jk+1;
            boundary_jk(N_jk,:) = [x1,j];
        end
        if test1(x2,x,xj,xk,xi,k,circle_jk,ncircle_jk) && min(sum(abs(boundary_jk(1:N_jk,1:3)-ones(N_jk,1)*x2),2))>epsilon
            N_jk = N_jk+1;
            boundary_jk(N_jk,:) = [x2,j];
        end
    end
end

for k=1:ncircle_ki
    cki=circle_ki(k,1:3);
    rki=circle_ki(k,4);
    j=circle_ki(k,5);
    if norm(cki-x)<probe+rki && norm(cki-x)>probe-rki
        [x1,x2]=intersection_circle(x,probe,cki,rki,nki);
        if test1(x1,x,xk,xi,xj,k,circle_ki,ncircle_ki) && min(sum(abs(boundary_ki(1:N_ki,1:3)-ones(N_ki,1)*x1),2))>epsilon
            N_ki=N_ki+1;
            boundary_ki(N_ki,:)=[x1,j];
        end
        if test1(x2,x,xk,xi,xj,k,circle_ki,ncircle_ki) && min(sum(abs(boundary_ki(1:N_ki,1:3)-ones(N_ki,1)*x2),2))>epsilon
            N_ki=N_ki+1;
            boundary_ki(N_ki,:)=[x2,j];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Compute I_interior recording all the intersection points interior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_interior=zeros(Kn*(Kn-1),5);
N_interior=0;

I_probe0 = zeros(N_ij+N_jk+N_ki,3);
I_probe0(1:N_ij,:)=boundary_ij(1:N_ij,1:3);
I_probe0(N_ij+1:N_ij+N_jk,:)=boundary_jk(1:N_jk,1:3);
I_probe0(N_ij+N_jk+1:N_ij+N_jk+N_ki,:)=boundary_ki(1:N_ki,1:3);

for j=1:Kn
    j0=hight_set(K(j));%%j0-th intersection point
    x0=I(j0,:);%%center of the probe
    for k=j+1:Kn
        k0=hight_set(K(k));
        y0=I(k0,:);
        %%%x,x0,y0
        [x1,x2,true]=intersection(x,x0,y0,probe);
        if true==1
           if test2(x1,x,xi,xj,xk) %%%%test if x1 is interior                    
                   true_test2 = 1; %% test if x1 is contained by another probe
                   for t = 1:Kn
                       if t~=j && t~= k
                           t0=hight_set(K(t));
                           z0=I(t0,:);
                           if norm(x1-z0)<probe
                               true_test2 = 0;
                               break;
                           end
                       end
                   end
                   if true_test2 == 1
                       N_interior=N_interior+1;
                       I_interior(N_interior,:)=[x1,j,k];%%%the intersection point on the probe might be eaten by another, but we record it
                   end
                %end
           end
           
           if test2(x2,x,xi,xj,xk) %%%test if x1 is interior
               true_test2 = 1;
               for t = 1:Kn
                   if t~=j && t~= k
                       t0=hight_set(K(t));
                       z0=I(t0,:);
                       if norm(x2-z0)<probe
                           true_test2 = 0;
                           break;
                       end
                   end
               end
               if true_test2 == 1
                   N_interior=N_interior+1;
                   I_interior(N_interior,:)=[x2,j,k];
               end
           end
           
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Compute the segments on the boudaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
boundary_ij=sort_boundary(x,boundary_ij,N_ij,nij);%%% associated with x,probe,nij
boundary_jk=sort_boundary(x,boundary_jk,N_jk,njk);
boundary_ki=sort_boundary(x,boundary_ki,N_ki,nki);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Compute the interior segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_boundary=N_ij+N_jk+N_ki;
N_probe=N_ij+N_jk+N_ki+N_interior;
I_probe=zeros(N_probe,3);
I_probe(1:N_ij,:)=boundary_ij(1:N_ij,1:3);
I_probe(N_ij+1:N_ij+N_jk,:)=boundary_jk(1:N_jk,1:3);
I_probe(N_ij+N_jk+1:N_ij+N_jk+N_ki,:)=boundary_ki(1:N_ki,1:3);
I_probe(N_boundary+1:N_probe,:)=I_interior(1:N_interior,1:3);

atom=zeros(Kn,10);%%at most 10 intersection points,it records all the I_probe on the K(j)-th probe
Natom=zeros(Kn,1);

for j=1:Kn
   for k=1:N_probe
       if k<=N_ij
           k0=k;
           if boundary_ij(k0,4)==j
               Natom(j)=Natom(j)+1;
               atom(j,Natom(j))=k;
           end
       elseif k<=N_ij+N_jk
           k0=k-N_ij;
           if boundary_jk(k0,4)==j
               Natom(j)=Natom(j)+1;
               atom(j,Natom(j))=k;
           end
       elseif k<=N_boundary
           k0=k-N_ij-N_jk;
           if boundary_ki(k0,4)==j
               Natom(j)=Natom(j)+1;
               atom(j,Natom(j))=k;
           end
       else
           k0=k-N_boundary;
           if I_interior(k0,4)==j
               Natom(j)=Natom(j)+1;
               atom(j,Natom(j))=k;
           end
           if I_interior(k0,5)==j
               Natom(j)=Natom(j)+1;
               atom(j,Natom(j))=k;
           end
       end
   end
end


%%%%consider the interior circles
               
global arg_eSAS;
global V_cSES;
global V_eSES;
circle_interior=zeros(Kn,7);
N_circle=0;
global index_I_circle;
index_I_circle = zeros(Kn,1);

for j=1:Kn
    if Natom(j)==0
        x1=0.5*(I(hight_set(K(j)),:)+x);
        if test2(x1,x,xi,xj,xk) %% there is no intersection points,which means that test
            r1=sqrt(probe^2-norm(x-x1)^2);
            n1=direction(i0,1)*(x1-x)/norm(x1-x); %%% clockwise rotation
            true1=1;
            for k=1:Kn
                dist1 = norm(x1-I(hight_set(K(k)),:));
                
                if j~=k && norm(x1-I(hight_set(K(k)),:))<probe
                    true1=0;
                    break;
                end
            end
            if true1==1
                N_circle=N_circle+1;
                circle_interior(N_circle,:)=[x1,n1,r1];%%%[c,n,r]
                
                index_I_circle(N_circle,1) = hight_set(K(j));
 
                if arg_eSAS==1
                    V_eSES=V_eSES-pi*r1^2*norm(x1-x)/3;
                else
                    V_cSES=V_cSES-pi*r1^2*norm(x1-x)/3;
                end
            end
        end
    end
end

circle_interior = circle_interior(1:N_circle,:);
index_I_circle = index_I_circle(1:N_circle,1);


%%%%%consider the interior segments

segment=zeros(N_probe,10);%%%segment records all the interior and boundary segments 
for j=1:2:(N_boundary-1)
    %%%%[c,n,r,angle,k1,k2],k1,k2 are index of I_probe
    if j<=N_ij
        uj=I_probe(j,:)-x;
        vj=I_probe(j+1,:)-x;
        aj=acos(uj*vj'/(probe^2));
        
        if j<N_ij-1
            j1=j+1;
        else
            j1=j+2;
        end
        segment((j+1)/2,:)=[x,nij,probe,aj,j,j1];%%%right-hand rule
    elseif j<=N_ij+N_jk
        uj=I_probe(j,:)-x;
        vj=I_probe(j+1,:)-x;
        aj=acos(uj*vj'/(probe^2));
        
        if j<N_ij+N_jk-1
            j1=j+1;
        else
            j1=j+2;
        end
        segment((j+1)/2,:)=[x,njk,probe,aj,j,j1];%%%right-hand rule
    else
        uj=I_probe(j,:)-x;
        vj=I_probe(j+1,:)-x;
        aj=acos(uj*vj'/(probe^2));
        
        if j<N_boundary-1
            j1=j+1;
        else
            j1=1;
        end
        segment((j+1)/2,:)=[x,nki,probe,aj,j,j1];%%%right-hand rule
    end
end

N_segment=N_boundary/2;

%%%%%%%%%% interior
global index_I_segment;
index_I_segment=zeros(N_probe,1);

for j=1:Kn
    if Natom(j)>0

        x1=0.5*(I(hight_set(K(j)),:)+x);
        r1=sqrt(probe^2-norm(x-x1)^2);
        n1=direction(i0,1)*(I(hight_set(K(j)),:)-x)/norm(I(hight_set(K(j)),:)-x); %%% clockwise rotation,if direct = 1, pointing x1 ; else, pointing x
        
        [A,K1,K2,N_new]=sort_segment(j,atom(j,1:Natom(j)),Natom(j),x1,n1,r1,I,hight_set,K,Kn,I_probe,probe,x,xi,xj,xk);
       
        N_segment=N_segment+N_new/2;
        if floor(N_segment)~= N_segment
            disp('Fail to create an interior segment on concave patch!')
        end
        
        segment(N_segment-floor(N_new/2)+1:N_segment,:)=[ones(floor(N_new/2),1)*[x1,n1,r1],A,K1,K2];%%%%[c,n,r,angle,k1,k2],k1,k2 are index of I_probe
        index_I_segment(N_segment-floor(N_new/2)+1:N_segment,1) = hight_set(K(j))*ones(floor(N_new/2),1);
    end
end

segment = segment(1:N_segment,:);
index_I_segment = index_I_segment(1:N_segment,1);

global Para Figs;

if N_boundary>6 || N_segment+N_circle>3
    if arg_eSAS 
        if Para.arg_viz == 1 && Para.viz_sing == 1
            % viz singular arcs
            figure(Figs.ext)
            visu_circlesegment(circle_interior,N_circle,segment,N_segment,I_probe,N_boundary);
            
        end
    end
   
end

global NB;
NB = [N_ij/2,N_jk/2,N_ki/2];
%visu_circlesegment(circle_interior,N_circle,segment,N_segment,I_probe);

construct_concavepat(i0,atoms,I,I_probe,N_probe,segment,N_segment,circle_interior,N_circle,direction(i0,1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%visualize the interior arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visu_circlesegment(circle_interior,N_circle,segment,N_segment,I_probe,N_boundary)

% visu singlar arcs
ind_start = N_boundary/2+1;

division=80;
for i= ind_start:N_segment
    v1=I_probe(segment(i,9),:)-segment(i,1:3);
    v1=v1/norm(v1);
    n=segment(i,4:6);
    v2=[v1(2)*n(3)-v1(3)*n(2),v1(3)*n(1)-v1(1)*n(3),v1(1)*n(2)-v1(2)*n(1)];
    v2=-v2;%%% always clockwise

    alpha0=segment(i,8);
    theta=linspace(0,alpha0,division)';   
    P=ones(division,1)*segment(i,1:3)+segment(i,7)*cos(theta)*v1+segment(i,7)*sin(theta)*v2; 
    P=real(P);
    plot3(P(1:division,1),P(1:division,2),P(1:division,3),'g','LineWidth',2)%'color',cc(1+floor(power(angle(i)/(2*pi),1)*100),:));
    hold on;
    
    %str=strcat(num2str(i));
    %text(P(division/2,1),P(division/2,2),P(division/2,3),str);
    
end

for i=1:N_circle
    
    Cc=circle_interior(i,1:3);
    Cn=circle_interior(i,4:6);
    Cr=circle_interior(i,7);
    
    [v1,v2]=orthogonalvectors(Cn);
    alpha=2*pi;
    theta=linspace(0,alpha,division)';   
    P=ones(division,1)*Cc+Cr*cos(theta)*v1+Cr*sin(theta)*v2;
    plot3(P(:,1),P(:,2),P(:,3),'g','LineWidth',2);

    hold on;
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%sort the interior intersection points to get arcs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,K1,K2,N_new]=sort_segment(j,atom,N,c,n,r,I,hight_set,K,Kn,I_probe,probe,x,xi,xj,xk)
area_left=pi*r^2;
epsilon=10^-10;
point=I_probe(atom',:);%%atoms records all the index of intersetion points
alpha1=zeros(1,N);
x0=point(1,1:3);%%start point
for k=2:N
    alpha1(k)=alpha(x0-c,point(k,1:3)-c,n);
end
[alpha2,order]=sort(alpha1);

A=zeros(floor(N/2),1);
K1=zeros(floor(N/2),1);
K2=zeros(floor(N/2),1);
N_A=0;
K0=[];
for k=1:N
    if norm(point(order(k),:)-point(order(mod(k,N)+1),:))>epsilon
        true=1;
        x_middle=0.5*(point(order(k),:)+point(order(mod(k,N)+1),:));
        if norm(x_middle-c)~=0
            x_middle=c+(x_middle-c)/norm(x_middle-c)*r;
        else
            v1=point(order(k),:)-c;
            u=-[v1(2)*n(3)-v1(3)*n(2),v1(3)*n(1)-v1(1)*n(3),v1(1)*n(2)-v1(2)*n(1)];
            x_middle=c+r*u/norm(u);
        end
        
        if alpha(point(order(k),:)-c,x_middle-c,n)>pi
            x_middle=2*c-x_middle;
        end
        
        for i=1:Kn
            center=I(hight_set(K(i)),:);
            if i~=j && norm(x_middle-center)<probe 
                true=0;
                break;
            end
        end
        if true==1 && test2(x_middle,x,xi,xj,xk)
            N_A=N_A+1;
            theta=alpha2(mod(k,N)+1)-alpha2(k);
            if theta<0
                A(N_A)=theta+2*pi;
            else
                A(N_A)=theta;
            end
            K1(N_A)=atom(order(k));
            K2(N_A)=atom(order(mod(k,N)+1));
            
            K0=[K0;k];
            
        end
        

    end
    
end

N_new=2*N_A;
A=A(1:N_A,:);
K1=K1(1:N_A,:);
K2=K2(1:N_A,:);

for i=1:N_A
    k=K0(i);
    k0=K0(mod(i-2,N_A)+1);
    k1=mod(k0,N)+1;
    angle_arc=alpha2(k)-alpha2(k1);
    if angle_arc<pi
        area_arc=angle_arc*r^2/2-r*sin(angle_arc/2)*r*cos(angle_arc/2);
    else
        area_arc=angle_arc*r^2/2-r*sin(angle_arc/2)*r*cos(angle_arc/2);
    end
    area_left=area_left-area_arc;
end

volume_delete=area_left*norm(c-x)/3;
global arg_eSAS;
global V_cSES;
global V_eSES;
if arg_eSAS ==1
    V_eSES=V_eSES-volume_delete;
else
    V_cSES=V_cSES-volume_delete;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%sort the intersection points to get arcs on the boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boundary=sort_boundary(c,boundary,N,n)
alpha1=zeros(1,N);
x0=boundary(1,1:3);
for k=2:N
    alpha1(k)=alpha(x0-c,boundary(k,1:3)-c,n);
end
[alpha2,order]=sort(alpha1);
boundary=boundary(order,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%test if the point x1 is inside the tretrahedron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function true=test2(x1,x,xi,xj,xk)
true=0;
D0=sign(det([xi-x;xj-x;xk-x]));
D1=sign(det([x1-x;xj-x;xk-x]));
D2=sign(det([xi-x;x1-x;xk-x]));
D3=sign(det([xi-x;xj-x;x1-x]));
if D0==D1 && D0==D2 && D0==D3
   true=1;
   return;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%test if x1 is on the boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function true=test1(x1,x,xi,xj,xk,k,circle,ncircle)
%%test if x1 is on the boundary (one arc)
epsilon=10^(-10);
true=1;
for i=1:ncircle
   if i~=k
       if norm(x1-circle(i,1:3))<=circle(i,4)-epsilon
           true=0;
           return;
       end
   end
end

if true==1
   D0=sign(det([xi-x;xj-x;xk-x]));
   D1=sign(det([x1-x;xj-x;xk-x]));
   D2=sign(det([xi-x;x1-x;xk-x]));
   if D0~=D1 || D0~=D2
       true=0;
       return;
   end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%compute the intersection points of three probes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,x2,true]=intersection(c1,c2,c3,r)
A=0.5*(c1+c2);
B=0.5*(c2+c3);
cij=c1-c2;
ckj=c3-c2;
cross=[cij(2)*ckj(3)-cij(3)*ckj(2),cij(3)*ckj(1)-cij(1)*ckj(3),cij(1)*ckj(2)-cij(2)*ckj(1)];
n=cross/norm(cross);
u=(c3-c2)-(c3-c2)*(c1-c2)'*(c1-c2)/(norm(c1-c2)^2);
t=(B-A)*(c3-c2)'/(u*(c3-c2)');
X1=A+t*u;
c=-(norm(X1-c1))^2+r^2;
if c>0
   x1=X1-sqrt(c)*n;
   x2=X1+sqrt(c)*n;
   true=1;
else
   x1=[0,0,0];
   x2=[0,0,0];
   true=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%compute the intersection point of two circles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1,x2]=intersection_circle(c1,r1,c2,r2,n)%%%x1, x2 are the two intersection of two circles
d=norm(c1-c2);
t=(r1^2-r2^2+d^2)/(2*d);
v=(c2-c1)/norm(c2-c1);
O=c1+t*v;

u=[n(2)*v(3)-n(3)*v(2),n(3)*v(1)-n(1)*v(3),n(1)*v(2)-n(2)*v(1)];
u=u/norm(u);
x1=O+sqrt(r1^2-t^2)*u;
x2=O-sqrt(r1^2-t^2)*u;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%compute the angle between two points on the same circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha=alpha(u,v,n)%the angle between two neighbor edges e and f, countorclockwise direction
t=sign(det([u;v;n]));%clockwise
if t>0
    alpha=acos(u*v'/(norm(u)*norm(v)));
else
    alpha=2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end