function VisualizeSlice

%names = ["1B17","1CRN","1ETN","1GZI","1YJO","3DIK","101M","3DIK"];
filename = 'benzene';
np = 100;
Rp = 0;

SliceInfo{1,1} = 'z';
SliceInfo{2,1} = -8.560666666666666;

[Fses,Fsas,X,Y,Z] = FsesBox(Rp,np,filename,SliceInfo);

figure(1)
clf
surf(X,Y,Z,Fses)
view(0,90)
colorbar
title('Distance to SES')
movegui('northwest')

figure(2)
clf
surf(X,Y,Z,Fsas)
view(0,90)
colorbar
title('Distance to SAS')
movegui('northwest')
