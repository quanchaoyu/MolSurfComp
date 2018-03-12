%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% mesh a spherical patch %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Nt,Ae,Nae,P,Np] = advancing_front_approach(c_sphere,r_sphere,N,T,Nt,Ae,Nae,P,Np,d,tolerance)
%%%%initially, T=[],Nt=0;

global Para;

global active;
global nactive;

Rp = Para.radius_probe;

theta=0.25;
alpha1=5/3*pi;%the angle condition parameters
h=d*sqrt(3)/2;

for k=1:N 
     %visu(P,T,Nt,'r')
     
     if Nt>floor(N/4*3)
         %disp('The mesh of the spherical patch is possibly overlapped!');
         %visu(P,Np,T,Nt)
     end
     
     if (mod(Nt,10)==1)
         
     end
     
    if Nae < 3
        Ae(1:Nae,:)=zeros(Nae,2);%empty
        Nae=0;
        break;
    end

    if Nae==3 %% the end conditon
        %number of triangles
        T=[T;Ae(1,1),Ae(1,2),Ae(2,2)];%add a triangle
        Nt=Nt+1;
        
        Nae=0;
        Ae(1:3,:)=zeros(3,2);%empty
        break;
    end
    
    if Nae>=4 || Nt==1 
        
        n1=normal_sphere(c_sphere,P(Ae(1,1),:));
        n2=normal_sphere(c_sphere,P(Ae(1,2),:));
        angle1=angle_sphere(Ae(Nae,:),Ae(1,:),n1,P);%%the left angle
        angle2=angle_sphere(Ae(1,:),Ae(2,:),n2,P);%%the right angle
        
        ne=ne_sphere(Ae(1,:),n1,n2,P);%%the tangent vector outwards the meshed region
        p=testpoint_sphere(Ae(1,:),ne,P,h);%%p is the test point
        x=map_sphere(c_sphere,r_sphere,p,Rp);%% x is on the sphere
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%dist and edgedist condition (including neighbors)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        index_np=0;%% index_np records the number of the point in P
        index_nactive=0;%% index_nactive records the number of activefront
        PA1 = P(Ae(1,1),:);
        PA2 = P(Ae(1,2),:);
        Dist1 = sqrt(sum((ones(Nae,1)*P(Ae(1,1),:)-P(Ae(:,1),:)).*(ones(Nae,1)*P(Ae(1,1),:)-P(Ae(:,1),:)),2));
        Dist2 = sqrt(sum((ones(Nae,1)*P(Ae(1,2),:)-P(Ae(:,1),:)).*(ones(Nae,1)*P(Ae(1,2),:)-P(Ae(:,1),:)),2));
        
        for m=3:Nae %%check if there exists a nonneighbor point on the loop close to the new point
%            dist1=sqrt(sum((PA1-P(Ae(m,1),:)).*(PA1-P(Ae(m,1),:))));
%            dist2=sqrt(sum((PA2-P(Ae(m,1),:)).*(PA2-P(Ae(m,1),:))));
%            dist1 = ;
%            dist2 = Dist2(m);
%            dist1=norm(P(Ae(1,1),:)-P(Ae(m,1),:));
%            dist2=norm(P(Ae(1,2),:)-P(Ae(m,1),:));
            if (Dist1(m)<tolerance || Dist2(m)<tolerance || Dist1(m)+Dist2(m)<tolerance*theta+sqrt((PA1(1)-PA2(1))^2+(PA1(2)-PA2(2))^2+(PA1(3)-PA2(3))^2))
                if ne*(P(Ae(m,1),:)- PA1)'>0 && (Ae(1,1)~=Ae(m,1)&&Ae(1,2)~=Ae(m,1))
                    flag=1;
                    if m==3 && angle2<5*pi/4
                        flag=0;
                    elseif m==Nae && angle1<5*pi/4
                        flag=0;
                    elseif m>3 && m<Nae 
                        angle1_m=angle_vectors(P(Ae(Nae,1),:)-P(Ae(1,1),:),P(Ae(m,1),:)-P(Ae(1,1),:),n1);
                        angle2_m=angle_vectors(P(Ae(m,1),:)-P(Ae(2,1),:),P(Ae(2,2),:)-P(Ae(2,1),:),n2);
                        if angle1_m<angle1 || angle2_m<angle2
                           flag=0; 
                        end
                    end

                    if flag==1
                        if index_np==0
                            index_np=m;
                            edgedist=Dist1(m)+Dist2(m);%%the distance from the point to the edge
                        elseif Dist1(m)+Dist2(m)<edgedist
                            index_np=m;
                            edgedist=Dist1(m)+Dist2(m);
                        elseif Ae(m,1)==Ae(index_np,1)
                            n_m=(P(Ae(m,1),:)-c_sphere)/norm(P(Ae(m,1),:)-c_sphere);
                            a1=angle_vectors(P(Ae(1,1),:)-P(Ae(m,1),:),P(Ae(m,2),:)-P(Ae(m,1),:),n_m);
                            a2=angle_vectors(P(Ae(1,1),:)-P(Ae(index_np,1),:),P(Ae(index_np,2),:)-P(Ae(index_np,1),:),n_m);
                            if m==Nae || a1>a2
                                index_np=m;
                            end
                        end
                    end
                end
            end
        end
        
        if nactive>0
            %pause
        end
        for i=1:nactive %check if there exists a nonneighbor point on other loop close to the new point
            if active(i).meshed==0
                Ae0=active(i).Ae;
                Nae0=active(i).Nae;
                for m=1:Nae0
                    dist1=norm(P(Ae(1,1),:)-P(Ae0(m,1),:));
                    dist2=norm(P(Ae(1,2),:)-P(Ae0(m,1),:));
                    if (dist1<tolerance || dist2<tolerance || dist1+dist2<tolerance*theta+norm(P(Ae(1,1),:)-P(Ae(1,2),:)))&& ne*(P(Ae0(m,1),:)-P(Ae(1,1),:))'>0
                        if index_np==0
                            index_np=m;
                            edgedist=dist1+dist2;%%the distance from the point to the first active edge
                            index_nactive=i;
                        elseif dist1+dist2<edgedist
                            index_np=m;
                            edgedist=dist1+dist2;
                            index_nactive=i;
                        end
                    end
                end
            end
         end


        if index_np>1
            if index_nactive==0
                if index_np>3&&index_np<Nae
                    [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor1_sphere(index_np,c_sphere,r_sphere,N,T,Nt,Ae,Nae,P,Np,d,tolerance);
                    break;
                elseif index_np==Nae
                    [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(1,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                elseif index_np==3
                    [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(2,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                end
                continue;
            elseif index_nactive>0
                [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor2_sphere(index_np,c_sphere,r_sphere,N,index_nactive,T,Nt,Ae,Nae,P,Np,d,tolerance);
                break;
            end
        end
               
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Test x: general dist condition with resp. to x
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        index_np=0;%% index_np records the number of the point in P
        index_nactive=0;%% index_nactive records the number of activefront

        PA1 = P(Ae(1,1),:);
        Dist = sqrt(sum((ones(Nae,1)*x-P(Ae(:,1),:)).*(ones(Nae,1)*x-P(Ae(:,1),:)),2)); % the distance from x to each front point
        
        for m=3:Nae %%check if there exists a nonneighbor point on the loop close to the new point
            %dist=norm(x-P(Ae(m,1),:));

            if (Dist(m)<tolerance || Dist(m-1)+Dist(m)<tolerance*theta+norm(P(Ae(m-1,1),:)-P(Ae(m-1,2),:))|| Dist(m)+Dist(mod(m,Nae)+1)<tolerance*theta+norm(P(Ae(m,1),:)-P(Ae(m,2),:))) && ne*(P(Ae(m,1),:)-PA1)'>0 &&(Ae(1,1)~=Ae(m,1)&&Ae(1,2)~=Ae(m,1))
                dist1=norm(P(Ae(1,1),:)-P(Ae(m,1),:));
                dist2=norm(P(Ae(1,2),:)-P(Ae(m,1),:));
                
                flag=1;
                if m==3 && angle2<5*pi/4
                    flag=0;
                elseif m==Nae && angle1<5*pi/4
                    flag=0;
                elseif m>3 && m<Nae 
                    angle1_m=angle_vectors(P(Ae(Nae,1),:)-P(Ae(1,1),:),P(Ae(m,1),:)-P(Ae(1,1),:),n1);
                    angle2_m=angle_vectors(P(Ae(m,1),:)-P(Ae(2,1),:),P(Ae(2,2),:)-P(Ae(2,1),:),n2);
                    if angle1_m<angle1 || angle2_m<angle2
                       flag=0; 
                    end
                end
                
                if flag==1
                    if index_np==0
                        index_np=m;
                        edgedist=dist1+dist2;%%the distance from the point to the edge
                    elseif dist1+dist2<edgedist
                        index_np=m;
                        edgedist=dist1+dist2;
                    elseif Ae(m,1)==Ae(index_np,1)
                        n_m=(P(Ae(m,1),:)-c_sphere)/norm(P(Ae(m,1),:)-c_sphere);
                        a1=angle_vectors(P(Ae(1,1),:)-P(Ae(m,1),:),P(Ae(m,2),:)-P(Ae(m,1),:),n_m);
                        a2=angle_vectors(P(Ae(1,1),:)-P(Ae(index_np,1),:),P(Ae(index_np,2),:)-P(Ae(index_np,1),:),n_m);
                        if m==Nae || a1>a2
                            index_np=m;
                        end
                    end
                end
                
            end
        end

        for i=1:nactive %check if there exists a nonneighbor point on other loop close to the new point
            if active(i).meshed==0
                Ae0=active(i).Ae;
                Nae0=active(i).Nae;
                for m=1:Nae0
                    dist=norm(x-P(Ae0(m,1),:));

                    if dist<tolerance && ne*(P(Ae0(m,1),:)-P(Ae(1,1),:))'>0
                        dist1=norm(P(Ae(1,1),:)-P(Ae0(m,1),:));
                        dist2=norm(P(Ae(1,2),:)-P(Ae0(m,1),:));

                        if index_np==0
                            index_np=m;
                            edgedist=dist1+dist2;%%the distance from the point to the edge
                            index_nactive=i;
                        elseif dist1+dist2<edgedist
                            index_np=m;
                            edgedist=dist1+dist2;
                            index_nactive=i;
                        end
                    end
                end
            end
         end


        if index_np>0
            if index_nactive==0
                if index_np>3&&index_np<Nae
                    [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor1_sphere(index_np,c_sphere,r_sphere,N,T,Nt,Ae,Nae,P,Np,d,tolerance);
                    break;
                elseif index_np==Nae
                    [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(1,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                elseif index_np==3
                    [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(2,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                end
                continue;
            elseif index_nactive>0
                [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor2_sphere(index_np,c_sphere,r_sphere,N,index_nactive,T,Nt,Ae,Nae,P,Np,d,tolerance);
                break;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%angle condition
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if angle1>alpha1&&angle2>alpha1 
            % in this case, we chose the one with the shortest edge distance
            
            edgedist_Nae=norm(P(Ae(Nae,1),:)-P(Ae(1,1),:))+norm(P(Ae(Nae,1),:)-P(Ae(1,2),:));
            edgedist_2=norm(P(Ae(2,1),:)-P(Ae(1,1),:))+norm(P(Ae(2,1),:)-P(Ae(1,2),:));
            if edgedist_Nae<edgedist_2
                [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(1,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
            else
                [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(2,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
            end
            continue;
            
        elseif angle1>alpha1
            angle1_2=angle_vectors(P(Ae(Nae,1),:)-P(Ae(2,1),:),P(Ae(2,2),:)-P(Ae(2,1),:),n2);
            if angle1_2>angle2
                [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(1,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                continue;
            end
        elseif angle2>alpha1
            angle1_1=angle_vectors(P(Ae(Nae,1),:)-P(Ae(1,1),:),P(Ae(2,2),:)-P(Ae(1,1),:),n1);
            if angle1_1>angle1
                [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(2,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance);
                continue;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Add a new point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                
        [T,Nt,Ae,Nae,P,Np]=addnewpoint_sphere(x,T,Nt,Ae,Nae,P,Np);
                
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Collapse the neighbor edge (left and right)%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Nt,Ae,Nae,P,Np]=collapse_neighbor_sphere(index_case,c_sphere,r_sphere,T,Nt,Ae,Nae,P,Np,tolerance)
%%%index_case decides we connect the left or the right edge
global Para;
Rp = Para.radius_probe;
gamma=1.8;%%a parameter which controles the length of edges
if index_case==1%%%Connect the left edge
    if norm(P(Ae(Nae,1),:)-P(Ae(1,2),:))<gamma*tolerance
        T=[T;Ae(Nae,1),Ae(1,2),Ae(1,1)];
        Nt=Nt+1;
        Ae=[Ae(2:Nae-1,:);Ae(Nae,1),Ae(1,2)];
        Nae=Nae-1;
    
    else
        p=0.5*(P(Ae(Nae,1),:)+P(Ae(1,2),:));
        x=map_sphere(c_sphere,r_sphere,p,Rp);
        
        if Nae==4 && norm(x-P(Ae(3,1),:))<10^-10
            T=[T;Ae(Nae,1),Ae(3,1),Ae(1,1);Ae(3,1),Ae(1,2),Ae(1,1)];
            Nt=Nt+2;
            Ae=[];
            Nae=0;
        else
            P=[P;x];
            Np=Np+1;
            T=[T;Ae(Nae,1),Np,Ae(1,1);Np,Ae(1,2),Ae(1,1)];
            Nt=Nt+2;
            Ae=[Ae(2:Nae-1,:);Ae(Nae,1),Np;Np,Ae(1,2)];
        end
        
    end
        
else%%%Connect the right edge
    if norm(P(Ae(1,1),:)-P(Ae(2,2),:))<gamma*tolerance
        T=[T;Ae(1,1),Ae(2,2),Ae(1,2)];
        Nt=Nt+1;
        Ae=[Ae(3:Nae,:);Ae(1,1),Ae(2,2)];
        Nae=Nae-1;
    
    else
        p=0.5*(P(Ae(1,1),:)+P(Ae(2,2),:));
        x=map_sphere(c_sphere,r_sphere,p,Rp);
        
        if Nae==4 && norm(x-P(Ae(Nae,1),:))<10^-10
            T=[T;Ae(1,1),Ae(Nae,1),Ae(2,1);Ae(Nae,1),Ae(2,2),Ae(2,1)];
            Nt=Nt+2;
            Ae=[];
            Nae=0;
        else
            P=[P;x];
            Np=Np+1;
            T=[T;Ae(1,1),Np,Ae(2,1);Np,Ae(2,2),Ae(2,1)];
            Nt=Nt+2;
            Ae=[Ae(3:Nae,:);Ae(1,1),Np;Np,Ae(2,2)];
        end
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%Collapse a nonneighbor point %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor1_sphere(index_np,c_sphere,r_sphere,N,T,Nt,Ae,Nae,P,Np,d,tolerance)
global Para;
Rp = Para.radius_probe;
gamma=1.8;

m=index_np;
if norm(P(Ae(m,1),:)-P(Ae(1,1),:))<=gamma*tolerance && norm(P(Ae(m,1),:)-P(Ae(1,2),:))<=gamma*tolerance
    T=[T;Ae(m,1),Ae(1,2),Ae(1,1)];
    Nt=Nt+1;
    Ae1=[Ae(2:m-1,:);Ae(m,1),Ae(1,2)];
    Nae1=m-1;
    Ae2=[Ae(m:Nae,:);Ae(1,1),Ae(m,1)];
    Nae2=Nae-m+2;
elseif norm(P(Ae(m,1),:)-P(Ae(1,1),:))>gamma*tolerance && norm(P(Ae(m,1),:)-P(Ae(1,2),:))>gamma*tolerance
    
    p=0.5*(P(Ae(1,1),:)+P(Ae(m,1),:));
    x=map_sphere(c_sphere,r_sphere,p,Rp);
    P=[P;x];
    Np=Np+1;
    p=0.5*(P(Ae(2,1),:)+P(Ae(m,1),:));
    x=map_sphere(c_sphere,r_sphere,p,Rp);
    P=[P;x];
    Np=Np+1;
        
    if norm(P(Ae(m,1),:)-P(Ae(1,1),:))>norm(P(Ae(m,1),:)-P(Ae(1,2),:))
        T=[T;Ae(1,1),Np-1,Ae(2,1);Np-1,Ae(m,1),Np;Np,Ae(2,1),Np-1];
        Nt=Nt+3;
    else
        T=[T;Ae(1,1),Np-1,Np;Np-1,Ae(m,1),Np;Np,Ae(2,1),Ae(1,1)];
        Nt=Nt+3;
    end
    
    Ae1=[Ae(2:m-1,:);Ae(m,1),Np;Np,Ae(1,2)];
    Nae1=m;

    Ae2=[Ae(m:Nae,:);Ae(1,1),Np-1;Np-1,Ae(m,1)];
    Nae2=Nae-m+3;
    
elseif norm(P(Ae(m,1),:)-P(Ae(1,1),:))>gamma*tolerance
    p=0.5*(P(Ae(1,1),:)+P(Ae(m,1),:));
    x=map_sphere(c_sphere,r_sphere,p,Rp);
    P=[P;x];
    Np=Np+1;
    
    T=[T;Ae(1,1),Np,Ae(2,1);Np,Ae(m,1),Ae(1,2)];
    Nt=Nt+2;
    
    Ae1=[Ae(2:m-1,:);Ae(m,1),Ae(1,2)];
    Nae1=m-1;
    Ae2=[Ae(m:Nae,:);Ae(1,1),Np;Np,Ae(m,1)];
    Nae2=Nae-m+3;
    
elseif norm(P(Ae(m,1),:)-P(Ae(1,2),:))>gamma*tolerance
    p=0.5*(P(Ae(2,1),:)+P(Ae(m,1),:));
    x=map_sphere(c_sphere,r_sphere,p,Rp);
    P=[P;x];
    Np=Np+1;
    
    T=[T;Ae(m,1),Np,Ae(1,1);Np,Ae(2,1),Ae(1,1)];
    Nt=Nt+2;
    
    Ae1=[Ae(2:m-1,:);Ae(m,1),Np;Np,Ae(1,2)];
    Nae1=m;
    Ae2=[Ae(m:Nae,:);Ae(1,1),Ae(m,1)];
    Nae2=Nae-m+2;
end

[T,Nt,Ae,Nae,P,Np] = advancing_front_approach(c_sphere,r_sphere,N,T,Nt,Ae1,Nae1,P,Np,d,tolerance);
[T,Nt,Ae,Nae,P,Np] = advancing_front_approach(c_sphere,r_sphere,N,T,Nt,Ae2,Nae2,P,Np,d,tolerance);
end

function [T,Nt,Ae,Nae,P,Np]=collapse_nonneighbor2_sphere(index_np,c_sphere,r_sphere,N,index_nactive,T,Nt,Ae,Nae,P,Np,d,tolerance)
global active;
active(index_nactive).meshed=1;
Ae0=active(index_nactive).Ae;
Nae0=active(index_nactive).Nae;

m=index_np;
T=[T;Ae0(m,1),Ae(1,2),Ae(1,1)];
Nt=Nt+1;

Ae1=[Ae(2:Nae,:);Ae(1,1),Ae0(m,1);Ae0(m:Nae0,:);Ae0(1:m-1,:);Ae0(m,1),Ae(1,2)];
Nae1=Nae+Nae0+1;

[T,Nt,Ae,Nae,P,Np] = advancing_front_approach(c_sphere,r_sphere,N,T,Nt,Ae1,Nae1,P,Np,d,tolerance);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Add a new point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,Nt,Ae,Nae,P,Np]=addnewpoint_sphere(x,T,Nt,Ae,Nae,P,Np)
P=[P;x];
Np=Np+1;
T=[T;Ae(1,1),Np,Ae(1,2)];
Nt=Nt+1;
Ae=[Ae(2:Nae,:);Ae(1,1),Np;Np,Ae(1,2)];
Nae=Nae+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Compute the angle between two neighbor edges%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha=angle_sphere(e,f,n,P)
%the angle between two neighbor edges e and f, countorclockwise direction
u=P(e(1),:)-P(e(2),:);
v=P(f(2),:)-P(f(1),:);
u = u-(u*n')*n;% project the vector u to the tangent vector of sphere
v = v-(v*n')*n;
t=sign(det([u;v;n]));
if t<0
    alpha = acos(u*v'/(sqrt(u(1)^2+u(2)^2+u(3)^2)*sqrt(v(1)^2+v(2)^2+v(3)^2)));
else
    alpha = 2*pi-acos(u*v'/(sqrt(u(1)^2+u(2)^2+u(3)^2)*sqrt(v(1)^2+v(2)^2+v(3)^2)));
end
end

function alpha = angle_vectors(u,v,n)
%the angle between two vectors, countorclockwise direction
u = u-(u*n')*n;% project the vector u to the tangent vector of sphere
v = v-(v*n')*n;
t = sign(det([u;v;n]));
if t < 0
    alpha = acos(u*v'/(norm(u)*norm(v)));
else
    alpha = 2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Compute the normal vector outwards the sphere%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normal=normal_sphere(c_sphere,x)
normal=(x-c_sphere)/sqrt((x(1)-c_sphere(1))^2+(x(2)-c_sphere(2))^2+(x(3)-c_sphere(3))^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Add a test point%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P=testpoint_sphere(e,ne,Pses,h)%add a new points
P=0.5*(Pses(e(1),:)+Pses(e(2),:))+h*ne;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Compute ne  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ne=ne_sphere(e,n1,n2,Pses)%%ne is the tangent vetor outwards the edge
u=Pses(e(2),:)-Pses(e(1),:);
v=0.5*n1+0.5*n2;
u=u-(u*v')*v; % project u to the correponding tangent vector to the sphere
ne=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3),u(1)*v(2)-u(2)*v(1)];
ne=ne/sqrt(ne(1)^2+ne(2)^2+ne(3)^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Map any point p to the sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x_sphere=map_sphere(c_sphere,r_sphere,p,Rp)%%maps any point p to x_sphere on the sphere
if r_sphere>Rp %%%% to mesh the spherical patch on the SES
    r=r_sphere-Rp;
else
    r=r_sphere;
end
u=(p-c_sphere)/norm(p-c_sphere);
x_sphere=c_sphere+r*u;
end

