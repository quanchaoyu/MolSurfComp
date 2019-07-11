function VisualizeSlice

%names = ["1B17","1CRN","1ETN","1GZI","1YJO","3DIK","101M","3DIK"];
filename = 'benzene';
np = 40;
Rp = 0;

[Fses,Fsas,X,Y,Z] = FsesBox(Rp,np,filename);

figure(1)
clf
ind = floor(np/2);
surf(squeeze(X(:,:,ind)),squeeze(Y(:,:,ind)),squeeze(Z(:,:,ind)),squeeze(Fses(:,:,ind)))
%surf(squeeze(X(:,:,ind)),squeeze(Y(:,:,ind)),squeeze(Z(:,:,ind)),max(0,squeeze(Fses(:,:,ind))))
view(0,90)
colorbar
title('Distance to SES')
movegui('northwest')

figure(2)
clf
ind = floor(np/2);
surf(squeeze(X(:,:,ind)),squeeze(Y(:,:,ind)),squeeze(Z(:,:,ind)),squeeze(Fsas(:,:,ind)))
%surf(squeeze(X(:,:,ind)),squeeze(Y(:,:,ind)),squeeze(Z(:,:,ind)),max(0,squeeze(Fsas(:,:,ind))))
view(0,90)
colorbar
title('Distance to SAS')
movegui('northwest')
