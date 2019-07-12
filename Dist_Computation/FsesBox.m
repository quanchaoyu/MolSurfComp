function [Fses,Fsas,X,Y,Z] = FsesBox(Rp,np,filename,SliceInfo)


MolSurfComp_WithInput(Rp,filename)
DataGlob;

disp(['Number of atoms: ',num2str(length(R))])

% find box containing all SAS-spheres
Rmax = max(R);
Mx = max(Geom.centers(:,1)+Rmax+Rp);
mx = min(Geom.centers(:,1)-Rmax-Rp);
My = max(Geom.centers(:,2)+Rmax+Rp);
my = min(Geom.centers(:,2)-Rmax-Rp);
Mz = max(Geom.centers(:,3)+Rmax+Rp);
mz = min(Geom.centers(:,3)-Rmax-Rp);
xl = mx:(Mx-mx)/(np-1):Mx;
yl = my:(My-my)/(np-1):My;
zl = mz:(Mz-mz)/(np-1):Mz;
if(nargin==4)
    switch SliceInfo{1,1}
        case char('x')
            [X,Y,Z] = meshgrid(SliceInfo{2,1},yl,zl);
        case char('y')
            [X,Y,Z] = meshgrid(xl,SliceInfo{2,1},zl);
        case char('z')
            [X,Y,Z] = meshgrid(xl,yl,SliceInfo{2,1});
        case char('xAv')
            [X,Y,Z] = meshgrid((Mx+mx)/2,yl,zl);
        case char('yAv')
            [X,Y,Z] = meshgrid(xl,(My+my)/2,zl);
        case char('zAv')
            [X,Y,Z] = meshgrid(xl,yl,(Mz+mz)/2);
    end            
else
    [X,Y,Z] = meshgrid(xl,yl,zl);
end
S = size(X);
Xv = reshape(X,prod(S),1);
Yv = reshape(Y,prod(S),1);
Zv = reshape(Z,prod(S),1);

    
% allocate memory
Fsesv = zeros(size(Xv));
Fsasv = zeros(size(Xv));

% compute the SAS-distance at each point of the grid
for i = 1:length(Xv)
    x = [Xv(i,1),Yv(i,1),Zv(i,1)];
    % evaluate distance function
    fsasx = Dist_SAS(x);
    % change sign of SAS-distance
    Fsasv(i,1) = -fsasx;
    % compute SES-distance
    Fsesv(i,1) = -fsasx - Rp;
end

% reshape 
if(size(S,2)==2)
    Fses = reshape(Fsesv,S(1),S(2));
    Fsas = reshape(Fsasv,S(1),S(2));
else
    Fses = reshape(Fsesv,S(1),S(2),S(3));
    Fsas = reshape(Fsasv,S(1),S(2),S(3));
end