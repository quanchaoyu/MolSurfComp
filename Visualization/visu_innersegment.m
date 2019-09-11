function visu_innersegment
DataGlob;

division=40; % each segment is divided to 100 
n=ncrasegment(:,1:3);
center=ncrasegment(:,4:6);
r=ncrasegment(:,7);
angle=ncrasegment(:,8);

hold on;
for i=1:nsegment
    v1=I(seg(i,3),:)-center(i,:);
    v1=v1/norm(v1);
    
    v2=[v1(2)*n(i,3)-v1(3)*n(i,2),v1(3)*n(i,1)-v1(1)*n(i,3),v1(1)*n(i,2)-v1(2)*n(i,1)];
    if seg(i,5)==1
        v2=-v2;
    end
    alpha=angle(i);
    theta=linspace(0,alpha,division)';   
    P=ones(division,1)*center(i,:)+r(i)*cos(theta)*v1+r(i)*sin(theta)*v2; 
    P=real(P);
    if ext_segment(i)==1
        %plot3(P(1:division,1),P(1:division,2),P(1:division,3),'y')%'color',[1,0.5,0])%'color',cc(1+floor(power(angle(i)/(2*pi),1)*100),:));
    else
        plot3(P(:,1),P(:,2),P(:,3),'r-','LineWidth',1);
        plot3(P([1,division],1),P([1,division],2),P([1,division],3),'r.','MarkerSize',8);
    end
    
    %plot3(P([1,division],1),P([1,division],2),P([1,division],3),'b.','MarkerSize',20);

end

for k=1:ncircle
    Cc=circle(k,3:5);
    Cn=circle(k,6:8);
    Cr=circle(k,9);
    
    [v1,v2]=orthogonalvectors(Cn);
    alpha=2*pi;
    theta=linspace(0,alpha,division)';   
    P=ones(division,1)*Cc+Cr*cos(theta)*v1+Cr*sin(theta)*v2;
    if ext_circle(k)==1
        %plot3(P(:,1),P(:,2),P(:,3),'color',[1,0.5,0]);
    else
        plot3(P(:,1),P(:,2),P(:,3),'m');
    end
end

end
