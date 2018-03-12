%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%compute the area of a spherical patch%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=area_spherical(c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Gauss-Bonnet Theorem%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%c_sphere denotes the sphere center, r_sphere denotes the sphere radius
%%%loops,loopsize,nloops denote the composition of this spherical patches
%%%[c,n,r,spoint,angle]=segment,n outwards the sphere,angle is clockwise 
%%%segment denote all the cirulcar arc on the boundary of this spherical patch
%%%the direction of an edge is defined to satisfy the right hand rule of n
%%%circle=[c,n,r],pointing outside the sphere
%%%there is only one patch with patchsize loops
Ke_sum=0;
Euler=2;%%Euler Characteristic
Kv_sum=0;

%patchesize=1;

for i=1:patchesize
    k=patches(i);%%the k-th loop or the -k-th circle
    if k>0
        for j=1:loopsize(k)
            k0=loops(k,j);
            phi=segment(k0,11);
            r=segment(k0,7);
            c=segment(k0,1:3);
            n=segment(k0,4:6); %% n points outwards from c_sphere to the center of the other sphere
            spoint=segment(k0,8:10);
            if (c-c_sphere)*n'>0
                Ke_sum=Ke_sum-sqrt(r_sphere^2-r^2)/r_sphere*phi;
            else
                Ke_sum=Ke_sum+sqrt(r_sphere^2-r^2)/r_sphere*phi;
            end
            
            %%%%%%%%%%%
            c1=segment(loops(k,mod(j-2,loopsize(k))+1),1:3);
            n1=segment(loops(k,mod(j-2,loopsize(k))+1),4:6);
            v1=spoint-c1;
            v2=spoint-c;
            
            u1=[n1(2)*v1(3)-n1(3)*v1(2),n1(3)*v1(1)-n1(1)*v1(3),n1(1)*v1(2)-n1(2)*v1(1)];
            u2=[n(2)*v2(3)-n(3)*v2(2),n(3)*v2(1)-n(1)*v2(3),n(1)*v2(2)-n(2)*v2(1)];
            
            alpha=acos(u1*u2'/(norm(u1)*norm(u2)));
            Kv_sum=Kv_sum+alpha;
            %%%%%%%%%%%
        end
        
        Euler=Euler-1;
    else
        k=-k;
        r=circle(k,7);
        c=circle(k,1:3);
        n=circle(k,4:6);
        if (c-c_sphere)*n'>0
            Ke_sum=Ke_sum-sqrt(r_sphere^2-r^2)/r_sphere*2*pi;
        else
            Ke_sum=Ke_sum+sqrt(r_sphere^2-r^2)/r_sphere*2*pi;
        end
        
        Euler=Euler-1;
    end
end

A=r_sphere^2*(2*pi*Euler-Ke_sum-Kv_sum);%%%the exact area of the volume

%A
end