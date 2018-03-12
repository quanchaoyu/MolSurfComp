function visu(P,T,Nt,c) 
%%P is the set of vertices, T is the set of triangles, c is the color

T=real(T);
P=real(P);

X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];

patch(real(X),real(Y),real(Z),c,'LineStyle','none','DiffuseStrength',0.8,'FaceAlpha',1,'FaceLighting','gouraud')%,'BackFaceLighting','unlit');
hold on;
end