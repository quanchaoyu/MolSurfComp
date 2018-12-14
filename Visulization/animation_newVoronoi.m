function animation_newVoronoi
opengl hardware

fig = openfig('new_Voro3.fig'); %('CirSeg_SES.fig'); 

nfigure = 0;
filename0 = 'new_Voronoi'; %'CirSeg_SES'; 

n = 100;
for i = 0:n
    filename = [filename0,'.gif'];
    nfigure = nfigure+1;
    %direction = [0 0 1];
    %rotate3d(fig,direction,10)
    view(150+360/n*i,20);

    axis([-12 12 -12 12 -12 12])
    if i == 0
        zoom(2)
    end
    drawnow

    set(gcf, 'units','normalized')
    %movegui(gcf,'center')
    frame = getframe(gcf);
    im0 = frame2im(frame);
    im = imresize(im0, [1000, 1000], 'bicubic');
    [imind,cm] = rgb2ind(im,256);
    if nfigure == 1;
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
    
    file = [filename0,'-',num2str(i+1)];
    print(file,'-dpng')
end

end