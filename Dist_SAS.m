function dist = Dist_SAS(point)
%
%   input: point is a 3D coordinate
%   output: dist is the distance between this point and SAS
%
%InitCode;
DataGlob;

%   find a ball containing point, otherwise, the point is outside SAS
index_ball = 0;
for i = 1:M
    if norm(point-C(i,:)) < R(i)+Rp
        index_ball = i;
        break;
    end
end

%   if point is outside SAS
if index_ball == 0
    %disp('Point outside the SAS:')
    for i = 1:M
        dist0 = norm(point-C(i,:))-R(i)-Rp;
        if i == 1
            dist = dist0;
        elseif dist0 < dist
            dist = dist0;
        end
    end
    %disp(['    distance = ',num2str(dist)])
    return;
end

%   if point is inside SAS, compute 
%       dist: distance
%       Xp: closest point
%       grad_fsas: gradient of fsas
%       laplace_fsas: laplace of fsas

%disp('Point inside the SAS:')
[dist,Xp,fses,~,grad_fsas,laplace_fsas] = fsas(point,index_ball,...
    C,R,M,Rp,I,inter,Row,Iijk,Ii,In,ext_I,seg,ncrasegment,...
    circle,circleindex,ncircleindex,ext_circle,...
    ext_patch,patches,patchesize,patches_index,loops,loopsize,loops_index);

% disp(['    distance = ',num2str(dist)])
% disp(['    gradient of distance function = ',num2str(grad_fsas)])
% disp(['    Laplacian of distance function = ',num2str(laplace_fsas)])
% disp(['    closest point = ',num2str(Xp)])

end