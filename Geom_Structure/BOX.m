function smallboxes=BOX(Geom,rmax)
% c
% c   create a large box which enclose all centers of atoms, composed of
% c   small boxes
% c
xmin = Geom.centers(1,1);
xmax = xmin;
ymin = Geom.centers(1,2);
ymax = ymin;
zmin = Geom.centers(1,3);
zmax = zmin;

for i = 2:Geom.M
    if Geom.centers(i,1) < xmin
        xmin = Geom.centers(i,1);
    elseif Geom.centers(i,1) > xmax
        xmax = Geom.centers(i,1);
    end
    if Geom.centers(i,2) < ymin
        ymin = Geom.centers(i,2);
    elseif Geom.centers(i,2) > ymax
        ymax = Geom.centers(i,2);
    end
    if Geom.centers(i,3) < zmin
        zmin = Geom.centers(i,3);
    elseif Geom.centers(i,3) > zmax
        zmax = Geom.centers(i,3);
    end
end

xmin = floor(xmin/rmax)*rmax;
xmax = (floor(xmax/rmax)+1)*rmax;
ymin = floor(ymin/rmax)*rmax;
ymax = (floor(ymax/rmax)+1)*rmax;
zmin = floor(zmin/rmax)*rmax;
zmax = (floor(zmax/rmax)+1)*rmax;

smallboxes = [floor(xmin/rmax),floor(xmax/rmax)+1;floor(ymin/rmax),floor(ymax/rmax)+1;floor(zmin/rmax),floor(zmax/rmax)+1];
end