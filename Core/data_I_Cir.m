function [high_I,hightvalue] = data_I_Cir

%c
%c   1) Compute the intersection points, loops (circle) and spherical patches on the SAS
%c   2) Input the computed data for meshing concave patches, convex patches and also for meshing toroidal patches using global variables
%c
global Para Geom Inter;
M = Geom.M;
C = Geom.centers;
R = Geom.R;
Rp = Para.radius_probe;
    
global DataI DataCir;

for i = 1:M
    Inter.M_int(i,1:Inter.num_int(i)) = sort(Inter.M_int(i,1:Inter.num_int(i)));
end

kmax = size(Inter.M_int,2); %  the max intersection number
inter = Inter.M_int; %  intersection matrix
Row = Inter.num_int; %  intersection number for each SAS-ball

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute all intersection points on the SAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = zeros(floor(M*kmax*(kmax-1)/3),3); % I records the coordinates of all intersection points on the SAS
high_I = ones(floor(M*kmax*(kmax-1)/3),1); % the i-th intersection point x with three atomic centers ci,cj,ck, if the high from x to plan (ci,cj,ck) is smallar than probe, high_I(i) = 1; otherwise, high_I(i) =0 
hightvalue = zeros(size(high_I));

Iijk = zeros(floor(M*kmax*(kmax-1)/3),3); % Iijk(i,:) gives the three indices (i,j,k) of SAS-balls forming the i-th intersection point
Ii = zeros(M,kmax*(kmax-1)); % Ii(i,:) gives indices of intersection points on the i-th SAS-ball
In = zeros(M,1); % In(i) gives the number of intersection points on the i-th SAS-ball
direction = zeros(floor(M*kmax*(kmax-1)/3),3); % if ci,cj,ck,x satisfies the right-hand rule, direction = 1; otherwise, direction = -1
s = 0;

circle = [];
I_circle = zeros(M,kmax,2*(kmax-1));% I_circle(i,row,:) records the indices of intersection point on the intersection circle of the i-th SAS-ball and the inter(row)-th SAS-ball
I_circle_num = zeros(M,kmax); % it records the number of intersection points on an intersection circle.

count = 0;
count_nonfree = 0;

for i = 1:M
    for row1 = 1:Row(i)
        j = inter(i,row1);
        if j > i
            
            circletest = 1;
            
            A = circlecenter(C(i,:),C(j,:),R(i)+Rp,R(j)+Rp); % compute the center of an intersection circle

            nij = (C(j,:)-C(i,:))/norm(C(j,:)-C(i,:));
            rij = sqrt((R(i)+Rp)^2-norm(C(i,:)-A)^2);
            
            for row = 1:Row(i) % Test if a circle is completedly covered by a sphere 
                k = inter(i,row);
                if k ~= j 
                    a = A-C(k,:);
                    normal = C(j,:)-C(i,:);
                    dist_a = sqrt(a(1)^2+a(2)^2+a(3)^2);
                    dist_normal = sqrt(normal(1)^2+normal(2)^2+normal(3)^2);
                    ctheta = (a(1)*normal(1)+a(2)*normal(2)+a(3)*normal(3))/(dist_a*dist_normal);%cos theta;
                    if (sqrt(1-ctheta^2)*dist_a+rij)^2+(ctheta*dist_a)^2-(R(k)+Rp)^2 < 0
                        circletest = 0; % In this case, the cicle is entirely covered by a sphere Sk
                        
                        
                        break;
                        
                    end
                end
            end
            
            if circletest == 1
                for row2 = 1:Row(i)
                    k = inter(i,row2);
                    
                    if k < j && circletest == 1% Test a circle is a whole circle
                        a = A-C(k,:);
                        normal = C(j,:)-C(i,:);
                        dist_a = sqrt(a(1)^2+a(2)^2+a(3)^2);
                        dist_normal = sqrt(normal(1)^2+normal(2)^2+normal(3)^2);
                        ctheta = (a(1)*normal(1)+a(2)*normal(2)+a(3)*normal(3))/(dist_a*dist_normal);%cos theta;
                        if (-sqrt(1-ctheta^2)*dist_a+rij)^2+(ctheta*dist_a)^2-(R(k)+Rp)^2 < 0
                            circletest = 0;
                        end
                    end
                    
                    
                    if k > j && any(inter(j,1:Row(j)) == k)
                        
                        row3 = find(inter(j,1:Row(j)) == k);
                        
                        B = circlecenter(C(k,:),C(j,:),R(k)+Rp,R(j)+Rp);
                        cij = C(i,:)-C(j,:);
                        ckj = C(k,:)-C(j,:);
                        cross = [cij(2)*ckj(3)-cij(3)*ckj(2),cij(3)*ckj(1)-cij(1)*ckj(3),cij(1)*ckj(2)-cij(2)*ckj(1)];
                        n = cross/norm(cross);
                        u = (C(k,:)-C(j,:))-(C(k,:)-C(j,:))*(C(i,:)-C(j,:))'*(C(i,:)-C(j,:))/((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2+(C(i,3)-C(j,3))^2);
                        t = (B-A)*(C(k,:)-C(j,:))'/(u*(C(k,:)-C(j,:))');
                        X1 = A+t*u;
                        c = -(norm(X1-C(i,:)))^2+(R(i)+Rp)^2;
                        
                        
                        if c > 0
                            circletest = 0;
                        
                            s1 = -sqrt(c);s2 = sqrt(c);% the hight of the tetrahedron (x,ci,cj,ck)
                            Y1 = X1+s1*n;
                            Y2 = X1+s2*n;
                            
                            true_hight = 1;
                            if s2 < Rp
                                true_hight = 0;
                            end

                            true1 = 1;
                            true2 = 1;
                            for row = 1:Row(i)
                                l = inter(i,row);
                                dist1 = sqrt((C(l,1)-Y1(1))^2+(C(l,2)-Y1(2))^2+(C(l,3)-Y1(3))^2);
                                
                                if dist1 < R(l)+Rp
                                    if l ~= j && l ~= k
                                        true1 = 0;
                                        break;
                                    end
                                end
                            end
                            for row = 1:Row(i)
                                l = inter(i,row);
                                dist2 = sqrt((C(l,1)-Y2(1))^2+(C(l,2)-Y2(2))^2+(C(l,3)-Y2(3))^2);
                                if dist2 < R(l)+Rp
                                    if l ~= j && l ~= k
                                        true2 = 0;
                                        break;
                                    end
                                end
                            end
                            

                            if true1==1

                                s=s+1;
                                I(s,:)=Y1;
                                Iijk(s,:)=[i,j,k];
                                direction(s,:)=[1,1,-1];

                                In(i) = In(i)+1;
                                Ii(i,In(i)) = s;
                                In(j) = In(j)+1;
                                Ii(j,In(j)) = s;
                                In(k) = In(k)+1;
                                Ii(k,In(k)) = s;
                                
                                I_circle_num(i,row1) = I_circle_num(i,row1)+1;
                                I_circle(i,row1,I_circle_num(i,row1)) = s;

                                I_circle_num(j,row3) = I_circle_num(j,row3)+1;
                                I_circle(j,row3,I_circle_num(j,row3)) = s;
                                
                                I_circle_num(i,row2) = I_circle_num(i,row2)+1;
                                I_circle(i,row2,I_circle_num(i,row2)) = s;
                                
                                hightvalue(s) = s2;
                                if true_hight == 0
                                    high_I(s) = 0;% it means that the hight of the tretrohedron (x,ci,cj,ck) is small so that singularity might occur
                                end
                            end
                            if true2 == 1

                                s = s+1;
                                I(s,:) = Y2;
                                Iijk(s,:) = [i,j,k];
                                direction(s,:) = [-1,-1,1];

                                In(i) = In(i)+1;
                                Ii(i,In(i)) = s;
                                In(j) = In(j)+1;
                                Ii(j,In(j)) = s;
                                In(k) = In(k)+1;
                                Ii(k,In(k)) = s;
                                
                                I_circle_num(i,row1) = I_circle_num(i,row1)+1;
                                I_circle(i,row1,I_circle_num(i,row1)) = s;

                                I_circle_num(j,row3) = I_circle_num(j,row3)+1;
                                I_circle(j,row3,I_circle_num(j,row3)) = s;
                                
                                I_circle_num(i,row2) = I_circle_num(i,row2)+1;
                                I_circle(i,row2,I_circle_num(i,row2)) = s;
                                
                                hightvalue(s) = s2;
                                if true_hight == 0
                                    high_I(s) = 0;% it means that the hight of the tretrohedron (x,ci,cj,ck) is small so that singularity might occur
                                end
                            end

                        end
                    end
                end
            end
            
            if circletest == 1
                circle = [circle;i,j,A,nij,rij];
            else
                count_nonfree = count_nonfree+1;
            end
            
        end
    end
end
I = I(1:s,:);
Iijk = Iijk(1:s,:);

Ii = Ii(:,1:max(In)); %   In(i): total number of intersection points on i-th SAS-sphere
direction = direction(1:s,:);

DataI.I = I;
DataI.nI = s;
DataI.direction = direction;
DataI.Iijk = Iijk;
DataI.Ii = Ii;
DataI.In = In;
DataI.I_circle = I_circle;
DataI.I_circle_num = I_circle_num;

high_I = high_I(1:s,1);
hightvalue = hightvalue(1:s,1);
DataI.high_I = high_I;
DataI.hightvalue = hightvalue;


cmax = 10; % max number of circles on the i-th sphere
ncircle = size(circle,1);
circleindex = zeros(M,cmax);% circleindex(i,:) -- indices of circles on the i-th SAS-sphere  
ncircleindex = zeros(M,1);% ncircleindex(i) -- number of circles on the i-th SAS-sphere

for k=1:ncircle
    i=circle(k,1);
    j=circle(k,2);
    ncircleindex(i,1)=ncircleindex(i,1)+1;
    circleindex(i,ncircleindex(i,1))=k;
    ncircleindex(j,1)=ncircleindex(j,1)+1;
    circleindex(j,ncircleindex(j,1))=k;
end

DataCir.circle = circle;
DataCir.ncircle = ncircle;
DataCir.circleindex = circleindex;
DataCir.ncircleindex = ncircleindex;

end


function O=circlecenter(c1,c2,r1,r2)
%c
%c  Compute the center of an intersection circle
%c
d=sqrt((c1(1)-c2(1))^2+(c1(2)-c2(2))^2+(c1(3)-c2(3))^2);
t=(r1^2-r2^2+d^2)/(2*d);
O=c1+(c2-c1)*t/d;
end