%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%%%%%%% Visualize toroidal patches%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visutorpat(P,Np,T,Nt,normal)
global Para;
al = Para.viz_SES(2);

T=real(T);
P=real(P);
global num_triangles;
num_triangles = num_triangles+Nt;

global Vertices;
global Triangle;
global NormalVects;
if Para.out_STL
    sv=size(Vertices,1);
    st=size(Triangle,1);
    Vertices=[Vertices(1:sv,:);P(1:Np,:),zeros(Np,1)];
    Triangle=[Triangle(1:st,:);T(1:Nt,:)+sv,ones(Nt,1)*508];
    T = T(1:Nt,:);
    NV = compute_NV_toroidal(T,P,normal);
    NormalVects = [NormalVects(1:st,:);NV];
end

X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];

patch(real(X),real(Y),real(Z),'y','LineStyle','none','FaceAlpha',al);
hold on;

end