function true = interiorloop_dist(point,ci,center,ri,i,R,probe,loop,loopsize,I,segment,ncrasegment)
% c
% c  check if a point on the sphere is in the interior of a loop, used for constructing a convex patch  
% c

true=0;

nearest=0;
for k=1:loopsize
    if segment(loop(k),1)==i
        j=segment(loop(k),2);
    else
        j=segment(loop(k),1);
    end
    
    v1=(center(j,:)-ci)/norm(center(j,:)-ci)*ri;
    v2=I(segment(loop(k),3),:)-ci;
    v3=point-ci;
    theta=real(acos(v3*v1'/(ri^2))-acos(v2*v1'/(ri^2)));
    
    if k==1 || theta<theta0
       theta0=theta;
       nearest=k;
    end
end


% if theta0 < -eps % the point is covered by a neighboring SAS-ball
%     true = 0;
%     return;
% end

K = zeros(10,1);
nK = 0;
if segment(loop(nearest),1)==i
    j_nearest = segment(loop(nearest),2);
else
    j_nearest = segment(loop(nearest),1);
end

for k = 1:loopsize
    if segment(loop(k),1)==i
        j=segment(loop(k),2);
    else
        j=segment(loop(k),1);
    end
    if j == j_nearest
        nK = nK+1;
        K(nK) = k;
    end
end

j = j_nearest;

v1=(center(j,:)-ci)/norm(center(j,:)-ci);%*ri;
v3=point-ci;

v=v3-(v1*v3')*v1;%/(ri^2);
    
for k0 = 1:nK
    k = K(k0);
    
    direct=segment(loop(k),5);
    n=ncrasegment(loop(k),1:3);
    A=ncrasegment(loop(k),4:6);
    alpha1=real(alpha(direct,I(segment(loop(k),3),:)-A,v,n));
    alpha2=real(alpha(direct,I(segment(loop(k),3),:)-A,I(segment(loop(k),4),:)-A,n));
    
    if segment(loop(k),1)==i
        j=segment(loop(k),2);
    else
        j=segment(loop(k),1);
    end
    if alpha1 < alpha2 && norm(center(j,:)-point)>=R(j)+probe
        true = 1;
        return;
    end
end

end

function alpha=alpha(direct,u,v,n)%the angle between two neighbor edges e and f, countorclockwise direction
t=sign(det([u;v;n]));%direct=1 means clockwise
if direct*t>0
    alpha=acos(u*v'/(norm(u)*norm(v)));
else
    alpha=2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end