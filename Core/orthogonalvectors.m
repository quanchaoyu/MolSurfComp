function [vector1,vector2] = orthogonalvectors(n)
%c
%c return the authorgonal vectors associated with the normal vector n
%c Remark: vector1,vector2,n obey right-hand rule; -n and n returns the
%c         same vector1
%c
a=n(1);b=n(2);c=n(3); %c solve a*x+b*y+c*z=0

if a~=0 && b~=0
    x=-sign(a)*b;y=sign(a)*a;z=0;
elseif b~=0 && c~=0
    x=0;y=-sign(b)*c;z=sign(b)*b;
elseif c~=0 && a~=0
    x=sign(c)*c;y=0;z=-sign(c)*a;
elseif a~=0
    x=0;y=1;z=0;
elseif b~=0
    x=0;y=0;z=1;
elseif c~=0
    x=1;y=0;z=0;
end

vector1=[x,y,z]/norm([x,y,z]);
vector2=[b*z-c*y,c*x-a*z,a*y-b*x];
vector2=vector2/norm(vector2);

end