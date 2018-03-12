function visusphere(c,r,alphaface)
%figure('Position',[700 100 600 600])
global nsphere;
global T_sphere;
global Nt_sphere;
global P_sphere;
global Np_sphere;

nsphere = 0;

if nsphere == 0
    readsphere;
end
nsphere=nsphere+1;
%{
rid=[];
for i=1:Nt
    for j=i+1:Nt
       if (T(i,1)==T(j,2)&&T(i,2)==T(j,3)&&T(i,3)==T(j,1))||(T(i,1)==T(j,3)&&T(i,2)==T(j,1)&&T(i,3)==T(j,2))
           rid=[rid,j];
       end
    end
end
rid
T(rid,:)=[];
Nt=Nt-length(rid);
%}

P=P_sphere.*(r/1.2);
P=P+ones(Np_sphere,1)*c;
T=T_sphere;
Nt=Nt_sphere;
Np = Np_sphere;

P=real(P);
%line([P(E(1:Ne,1),1)';P(E(1:Ne,2),1)'],[P(E(1:Ne,1),2)';P(E(1:Ne,2),2)'],[P(E(1:Ne,1),3)';P(E(1:Ne,2),3)'],'LineStyle','-');%,'Marker','.'
X=[P(T(1:Nt,1),1)';P(T(1:Nt,2),1)';P(T(1:Nt,3),1)'];
Y=[P(T(1:Nt,1),2)';P(T(1:Nt,2),2)';P(T(1:Nt,3),2)'];
Z=[P(T(1:Nt,1),3)';P(T(1:Nt,2),3)';P(T(1:Nt,3),3)'];
% for i = 1:Nt_sphere
%     patch(X(:,i),Y(:,i),Z(:,i),'y','FaceAlpha',0.1);%'DiffuseStrength',0.8,
%     refined_triangle(c,r,X(:,i),Y(:,i),Z(:,i),1);
%     hold on;
% end


patch(X,Y,Z,'y','LineStyle','none','FaceAlpha',alphaface);%'DiffuseStrength',0.8,
hold on;

end

function readsphere
global T_sphere;
global Nt_sphere;
global P_sphere;
global Np_sphere;

T_sphere=textread('T_sphere.txt');
Nt_sphere=textread('Nt_sphere.txt');
P_sphere=textread('Pses_sphere.txt');
Np_sphere=textread('Np_sphere.txt');

end
