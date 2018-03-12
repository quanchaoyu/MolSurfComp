function NV = compute_NV(T,P,c_sphere,arg_NV)

% compute the normal vector for each triangle on spherical patches

N = size(T,1);
NV = zeros(N,3);
for i = 1:N
    a = P(T(i,1),:);
    b = P(T(i,2),:);
    c = P(T(i,3),:);
    
    V1 = b-a;
    V2 = c-a;
    
    nv = cross(V1,V2);
    nv = nv/norm(nv);
    
    s = sign(nv*(a-c_sphere)');
    
    if s ~= arg_NV
        nv = -nv;
    end
    
    NV(i,:) = nv;
end

end