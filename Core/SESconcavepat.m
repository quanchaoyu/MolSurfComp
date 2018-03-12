function SESconcavepat(s,I,Iijk,Ii,In,C,R,Rp,ext_I,hight,direction,intersection)
%c
%c  Compute Segments and Circles of a Concave Patch 
%c
global arg_eSAS;

global DataAV;

global V_cSAS;
global V_eSAS;
global V_cSES;
global V_eSES;

%c
%c  Compute the Volume of cSAS and eSAS
%c
V_cSAS=0;
V_eSAS=0;

for i=1:s
    
    ci=C(Iijk(i,1),:);
    cj=C(Iijk(i,2),:);
    ck=C(Iijk(i,3),:);
    x=I(i,:);
    direct=direction(i,1);
    ni=(ci-x)/norm(ci-x);
    nj=(cj-x)/norm(cj-x);
    nk=(ck-x)/norm(ck-x);
    nij=[ni(2)*nj(3)-ni(3)*nj(2),ni(3)*nj(1)-ni(1)*nj(3),ni(1)*nj(2)-ni(2)*nj(1)];
    njk=[nj(2)*nk(3)-nj(3)*nk(2),nj(3)*nk(1)-nj(1)*nk(3),nj(1)*nk(2)-nj(2)*nk(1)];
    nki=[nk(2)*ni(3)-nk(3)*ni(2),nk(3)*ni(1)-nk(1)*ni(3),nk(1)*ni(2)-nk(2)*ni(1)];
    nij=direct*nij/norm(nij);
    njk=direct*njk/norm(njk);
    nki=direct*nki/norm(nki);
    
    z1=(x(3)+ci(3)+cj(3))/3;
    n1=nij(3);
    d1=norm(x-ci);
    d2=norm(x-cj);
    d3=norm(ci-cj);
    p=(d1+d2+d3)/2;
    area1=sqrt(p*(p-d1)*(p-d2)*(p-d3));
    volume1=area1*n1*z1;

    z2=(x(3)+cj(3)+ck(3))/3;
    n2=njk(3);
    d1=norm(x-cj);
    d2=norm(x-ck);
    d3=norm(cj-ck);
    p=(d1+d2+d3)/2;
    area2=sqrt(p*(p-d1)*(p-d2)*(p-d3));
    volume2=area2*n2*z2;

    z3=(x(3)+ck(3)+ci(3))/3;
    n3=nki(3);
    d1=norm(x-ck);
    d2=norm(x-ci);
    d3=norm(ck-ci);
    p=(d1+d2+d3)/2;
    area3=sqrt(p*(p-d1)*(p-d2)*(p-d3));
    volume3=area3*n3*z3;

    if ext_I(i)==1
        V_cSAS=V_cSAS+volume1+volume2+volume3;
        V_eSAS=V_eSAS+volume1+volume2+volume3;

    else
        V_cSAS=V_cSAS+volume1+volume2+volume3;
    end
end
DataAV.Vcsas=DataAV.Vcsas+V_cSAS;
DataAV.Vesas=DataAV.Vesas+V_eSAS;

V_cSES=V_cSAS;
V_eSES=V_eSAS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eSAS case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hight_set=zeros(s-sum(hight),1);
nhight=0;

arg_eSAS=1;

global index_externalconcavepatch;
index_externalconcavepatch = 0;

inverse_hight = zeros(s,1);
for i=1:s % s is the number of intersection points
   if hight(i)==0 && ext_I(i)==1
       %  ext_I(i) == 1 means the i-th intersection point is on the eSAS
       nhight=nhight+1;
       hight_set(nhight,1)=i;
       inverse_hight(i) = nhight;
   end
end

probemax = 200;
neighbor_I=zeros(nhight,probemax);
N_neighbor=zeros(nhight,1);
for i = 1:nhight
    x1 = I(hight_set(i),:);
    i1 = Iijk(hight_set(i),1);
    for j = 1:intersection.num_int(i1)
        j1 = intersection.M_int(i1,j);
        for k = 1:In(j1)
            k1 = Ii(j1,k);
            x2 = I(k1,:);
            if k1~=hight_set(i) && norm(x1-x2)<2*Rp && any(inverse_hight(k1)==neighbor_I(i,:)) == 0
                N_neighbor(i) = N_neighbor(i)+1;
                neighbor_I(i,N_neighbor(i)) = inverse_hight(k1);
            end
        end
    end
end

for i=1:nhight
   if N_neighbor(i)>0 
      % N_neighbor(i)>0 means that hight_set(i) is a really singular concave patch
      %if ext_I
      data_concavepat(i,hight_set,neighbor_I(i,1:N_neighbor(i)),N_neighbor(i),I,Iijk,C,Rp,direction);
   else
      % in this case, the concave SES patch violates the hight condition, but is still triangle-shaped
      data_concavepat(i,hight_set,neighbor_I(i,1:N_neighbor(i)),N_neighbor(i),I,Iijk,C,Rp,direction);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cSAS case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global index_completeconcavepatch;
index_completeconcavepatch = 0;

hight_set=zeros(s-sum(hight),1);
nhight=0;
arg_eSAS=0;

inverse_hight = zeros(s,1);
for i=1:s % s is the number of intersection points
   if hight(i)==0 
       nhight=nhight+1;
       hight_set(nhight,1)=i;
       inverse_hight(i) = nhight;
   end
end

probemax = 200;
neighbor_I=zeros(nhight,probemax);
N_neighbor=zeros(nhight,1);
for i = 1:nhight
    x1 = I(hight_set(i),:);
    i1 = Iijk(hight_set(i),1);
    for j = 1:intersection.num_int(i1)
        j1 = intersection.M_int(i1,j);
        for k = 1:In(j1)
            k1 = Ii(j1,k);
            x2 = I(k1,:);
            if k1~=hight_set(i) && norm(x1-x2)<2*Rp && any(inverse_hight(k1)==neighbor_I(i,:)) == 0
                N_neighbor(i) = N_neighbor(i)+1;
                neighbor_I(i,N_neighbor(i)) = inverse_hight(k1);
            end
        end
    end
end

for i=1:nhight
   if N_neighbor(i)>0 
      % N_neighbor(i)>0 means that hight_set(i) is a really singular concave patch
      data_concavepat(i,hight_set,neighbor_I(i,1:N_neighbor(i)),N_neighbor(i),I,Iijk,C,Rp,direction);
   else
      % in this case, the concave SES patch violate the hight condition, but is still triangle-shaped
      data_concavepat(i,hight_set,neighbor_I(i,1:N_neighbor(i)),N_neighbor(i),I,Iijk,C,Rp,direction);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Concave Triangles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Para Figs;

arg_eSAS=1;

for i=1:s
    if hight(i)==1 
        segment0=zeros(3,11);
        
        % segment:[i,j,p1,p2,direct],ncrasegment->[c,n,r,spoint,angle]
        c=I(i,:);
        
        ni=(C(Iijk(i,1),:)-c)/norm(C(Iijk(i,1),:)-c);
        nj=(C(Iijk(i,2),:)-c)/norm(C(Iijk(i,2),:)-c);
        nk=(C(Iijk(i,3),:)-c)/norm(C(Iijk(i,3),:)-c);
        nij=[ni(2)*nj(3)-ni(3)*nj(2),ni(3)*nj(1)-ni(1)*nj(3),ni(1)*nj(2)-ni(2)*nj(1)];
        njk=[nj(2)*nk(3)-nj(3)*nk(2),nj(3)*nk(1)-nj(1)*nk(3),nj(1)*nk(2)-nj(2)*nk(1)];
        nki=[nk(2)*ni(3)-nk(3)*ni(2),nk(3)*ni(1)-nk(1)*ni(3),nk(1)*ni(2)-nk(2)*ni(1)];
        nij=nij/norm(nij);
        njk=njk/norm(njk);
        nki=nki/norm(nki);
        spoint_i=c+(C(Iijk(i,1),:)-c)*Rp/(R(Iijk(i,1))+Rp);
        spoint_j=c+(C(Iijk(i,2),:)-c)*Rp/(R(Iijk(i,2))+Rp);
        spoint_k=c+(C(Iijk(i,3),:)-c)*Rp/(R(Iijk(i,3))+Rp);
        angle_i=acos((spoint_i-c)*(spoint_j-c)'/Rp^2);
        angle_j=acos((spoint_j-c)*(spoint_k-c)'/Rp^2);
        angle_k=acos((spoint_k-c)*(spoint_i-c)'/Rp^2);
        loops=[1,2,3];
        if direction(i,1)==-1
            nij=-nij;njk=-njk;nki=-nki;
            si=spoint_i;
            spoint_i=spoint_j;
            spoint_j=spoint_k;
            spoint_k=si;
            loops=[1,3,2];
        end
        segment0(1,:)=[c,nij,Rp,spoint_i,angle_i];
        segment0(2,:)=[c,njk,Rp,spoint_j,angle_j];
        segment0(3,:)=[c,nki,Rp,spoint_k,angle_k];
        
        Area=area_spherical(I(i,:),Rp,loops,3,segment0,[],1,1);
        if ext_I(i)==1
            index_externalconcavepatch = index_externalconcavepatch+1;
            index_concavepatch = index_externalconcavepatch;
            atoms = Iijk(i,:);
            NB = [1,1,1];
            index_I = i;
            
            if Para.out_MolStrc
                output_SES_concavepatches(atoms,NB,zeros(3,1),[],index_I,index_concavepatch,1,I(i,:),Rp,loops,3,segment0,[],1,1);
            end
            
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext
                    figure(Figs.ext)
                    hold on;
                end
                mesh_sphpat(I(i,:),Rp,loops,3,segment0,[],1,1);%%% Just below: it is added the area of concave patch on the eSES
            end
            CompAreaVol_concave(I(i,:),Rp,loops,3,segment0,[],1,1); %% add!
            
            DataAV.Acses=DataAV.Acses+Area;
            V_cSES=V_cSES-Area*Rp/3;
        else
            
            if Para.arg_meshing(1) 
                if Para.arg_viz && Para.arg_int
                    figure(Figs.int)
                    hold on;
                end
                mesh_sphpat(I(i,:),Rp,loops,3,segment0,[],1,1);
            end
            
            DataAV.Acses=DataAV.Acses+Area;
            V_cSES=V_cSES-Area*Rp/3;

        end
        
        index_completeconcavepatch = index_completeconcavepatch+1;
        index_concavepatch = index_completeconcavepatch;
        atoms = Iijk(i,:);
        NB = [1,1,1];
        index_I = i;
        
        if Para.out_MolStrc
            output_SES_concavepatches(atoms,NB,zeros(3,1),[],index_I,index_concavepatch,0,I(i,:),Rp,loops,3,segment0,[],1,1);
        end
    end
end

DataAV.Vcses=DataAV.Vcses+V_cSES;
DataAV.Veses=DataAV.Veses+V_eSES;
end