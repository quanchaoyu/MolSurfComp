function VisualizeSlice

%names = ["1B17","1CRN","1ETN","1GZI","1YJO","3DIK","101M","3DIK"];
%names =
%["carbo","glutaredoxin","hiv-1-gp41","l-plectasin","ubch5b","vancomycin"];
file.name = '3spheres';
file.format = '.xyzr';

np = 200;
Rp = 1;

SliceInfo{1,1} = 'zAv';
SliceInfo{2,1} = -8.5;

[Fses,Fsas,X,Y,Z] = FsesBox(Rp,np,file,SliceInfo);

figure(1)
clf
hold on
surf(X,Y,Z,Fses,'EdgeColor','none')
[~,h] = contour(X,Y,Fses,[0 0],'b','LineWidth',3);
h.ContourZLevel = Z(1,1);
view(0,90)
colorbar
title('Distance to SES')
movegui('northwest')

figure(2)
clf
hold on
surf(X,Y,Z,Fsas,'EdgeColor','none')
[~,h] = contour(X,Y,Fsas,[0 0],'r','LineWidth',3);
h.ContourZLevel = Z(1,1);
view(0,90)
colorbar
title('Distance to SAS')
movegui('northwest')
