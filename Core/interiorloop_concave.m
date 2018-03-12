function true=interiorloop_concave(point,x,probe,loop,loopsize,I_probe,segment,direct)
% c
% c   check if a point on the sphere is in the interior of a loop, used for constructing a concave patch     
% c   [c,n,r,angle,k1,k2],k1,k2 are index of I_probe
% c

true=0;

nearest=0;
for k=1:loopsize
    v1=segment(loop(k),4:6)*probe;
    if direct==-1
        v1=-v1;
    end
    v2=I_probe(segment(loop(k),9),:)-x;
    v3=point-x;
    theta=acos(v3*v1'/(probe^2))-acos(v2*v1'/(probe^2));
    if k==1 || theta<theta0
       theta0=theta;
       nearest=k;
    end
end

if theta0 < 0 % the point is covered by a neighboring SAS-ball
    true = 0;
    return;
end

K = zeros(10,1);
nK = 0;

for k = 1:loopsize
    if  norm(segment(loop(k),1:3)-segment(loop(nearest),1:3))== 0 && abs(segment(loop(k),7)-segment(loop(nearest),7))== 0
        nK = nK+1;
        K(nK) = k;
    end
end

if nK == 0
    disp('error in interiorloop_concave');
end

for k0 = 1:nK
    k = K(k0);

    v1=segment(loop(k),4:6)*probe;
    if direct==-1
        v1=-v1;
    end
    v3=point-x;

    v=v3-(v1*v3')*v1/(probe^2);

    n=segment(loop(k),4:6);
    A=segment(loop(k),1:3);
    alpha1=real(alpha(I_probe(segment(loop(k),9),:)-A,v,n));
    alpha2=real(alpha(I_probe(segment(loop(k),9),:)-A,I_probe(segment(loop(k),10),:)-A,n));

    if alpha1<alpha2
        true=1;
        return;
    end
    
end

end

function alpha=alpha(u,v,n)%the angle between two neighbor edges e and f, countorclockwise direction
t=sign(det([u;v;n]));%direct=1 means clockwise
if t>0
    alpha=acos(u*v'/(norm(u)*norm(v)));
else
    alpha=2*pi-acos(u*v'/(norm(u)*norm(v)));
end
end