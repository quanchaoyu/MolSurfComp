function [Fses,Fsas,X,Y,Z] = FsesBox(Rp,np,filename)

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
[X,Y,Z] = meshgrid(xl,yl,zl);
Xv = reshape(X,np^3,1);
Yv = reshape(Y,np^3,1);
Zv = reshape(Z,np^3,1);

% allocate memory
Fsesv = zeros(size(Xv));
Fsasv = zeros(size(Xv));

% compute the SAS-distance at each point
for i = 1:length(Xv)
    x = [Xv(i,1),Yv(i,1),Zv(i,1)];
    fsasx = Dist_SAS(x);
    Fsasv(i,1) = -fsasx;
    Fsesv(i,1) = -fsasx - Rp;
end

% reshape 
Fses = reshape(Fsesv,np,np,np);
Fsas = reshape(Fsasv,np,np,np);
