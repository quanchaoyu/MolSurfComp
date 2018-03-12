function data_SEStorpat
% c
% c   Mesh toroidal patches with or without singularities
% c
% c   segment -- [i,j,k1,k2,direct], where (i,j) are the indices of two
% c              corresponding SAS-balls, (k1,k2) are indices the starting and ending
% c              points, direct gives the direction of the segment
% c   ncrasegment -- [normal vector,center,radius,angle]
% c   circle -- [i,j,center,normal,radius]
% c
global  Geom Para;
global DataI DataSeg DataCir;
I = DataI.I;
segment = DataSeg.segment;
nsegment = DataSeg.nsegment;
ncrasegment = DataSeg.ncrasegment;
circle = DataCir.circle;
ncircle = DataCir.ncircle;

Rp = Para.radius_probe;
C = Geom.centers;
R = Geom.R;

global DataAV;

global Ext Figs;
ext_segment = Ext.segment;
ext_circle = Ext.circle;

global arg_eSAS;
if arg_eSAS == 0
    ext_segment = ones(size(ext_segment));
    ext_circle = ones(size(ext_circle));
end

% if Rp == 0
%    return; 
% end


global index_toroidalpatch; % count the number of toroidal patches
index_toroidalpatch = 0;

%% Arc case

for i=1:nsegment
    
    r=ncrasegment(i,7);%  the radius of the i-th arc
    A=ncrasegment(i,4:6);
    n=ncrasegment(i,1:3);%  from ci to cj
    P_k1=I(segment(i,3),:);
    P_k2=I(segment(i,4),:);
    direct=segment(i,5);
    angle=ncrasegment(i,8);
    ri=R(segment(i,1));
    rj=R(segment(i,2));
    ci=C(segment(i,1),:);
    cj=C(segment(i,2),:);
    

    %%
    if r < Rp && (ci-A)*(cj-A)'<0
        A1=A-sqrt(Rp^2-r^2)*n;
        theta1=acos(r/(ri+Rp))-acos(r/Rp);
        c1=ci*Rp/(Rp+ri)+A*ri/(Rp+ri);
        r1=ri/(Rp+ri)*r;%%c1,r1 denote the center and the radius of the arc on the VdW sphere
        
        if ext_segment(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_cusp(c1,r1,A1,direct,n,angle,r,A,P_k1,P_k2,theta1,Rp,i);
            end
        else
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
                mesh_cusp(c1,r1,A1,direct,n,angle,r,A,P_k1,P_k2,theta1,Rp,i);
            end
        end
        
        A2=A+sqrt(Rp^2-r^2)*n;
        theta2=acos(r/(rj+Rp))-acos(r/Rp);
        
        c1=cj*Rp/(Rp+rj)+A*rj/(Rp+rj);
        r1=rj/(Rp+rj)*r;
        
        theta=acos(r/Rp);
        As=angle*Rp*(r*(theta1+theta2)-Rp*(sin(theta1+theta)+sin(theta2+theta)-2*sin(theta)));
        V_delete=angle*Rp^2*(r*(theta1+theta2)/2-Rp/3*(sin(theta1+theta)+sin(theta2+theta)-2*sin(theta)))+angle*r^2*sqrt(Rp^2-r^2)/3;
        Vs_sas=angle/2*r^2*norm(ci-cj)/3;
        Vs_ses=Vs_sas-V_delete;
        
        if ext_segment(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_cusp(c1,r1,A2,direct,n,angle,r,A,P_k1,P_k2,theta2,Rp,i);
            end
            DataAV.Acses=DataAV.Acses+As;
            DataAV.Aeses=DataAV.Aeses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            DataAV.Vesas=DataAV.Vesas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            DataAV.Veses=DataAV.Veses+Vs_ses;
        else
            DataAV.Acses=DataAV.Acses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
            	mesh_cusp(c1,r1,A2,direct,n,angle,r,A,P_k1,P_k2,theta2,Rp,i);
            end
        end
        
        %%%%% output the toroidal patch corresponding to the SAS-segment
        if Para.out_MolStrc
            sing = 1; 
            if direct == 1
                P1 = P_k1;
                P2 = P_k2;
                index_I1 = segment(i,3);
                index_I2 = segment(i,4);
            else
                P1 = P_k2;
                P2 = P_k1;
                index_I1 = segment(i,4);
                index_I2 = segment(i,3);
            end

            SASsegment = [A,n,r,P1,P2,angle]; % [c,n,r,s,e,angle], n points from ci to cj

            c_segment1 = ci+ri/(ri+Rp)*(A-ci);
            c_segment2 = cj+rj/(rj+Rp)*(A-cj);
            n_segment1 = -n; % the interior is on the right-hand along the segment
            n_segment2 = n;
            r_segment1 = r*ri/(ri+Rp);
            r_segment2 = r*rj/(rj+Rp);
            P1_segment1 = ci+ri/(ri+Rp)*(P2-ci);%%
            P2_segment1 = ci+ri/(ri+Rp)*(P1-ci);
            P1_segment2 = cj+rj/(rj+Rp)*(P1-cj);
            P2_segment2 = cj+rj/(rj+Rp)*(P2-cj);

            s1 = [c_segment1,n_segment1,r_segment1,P1_segment1,P2_segment1,angle,1]; % [c,n,r,s,e,angle,right-hand]
            s2 = [c_segment2,n_segment2,r_segment2,P1_segment2,P2_segment2,angle,1];

            ni1 = ci-P1;nj1 = cj-P1;
            n1_twosides = [ni1(2)*nj1(3)-ni1(3)*nj1(2),ni1(3)*nj1(1)-ni1(1)*nj1(3),ni1(1)*nj1(2)-ni1(2)*nj1(1)];
            n1_twosides = n1_twosides/norm(n1_twosides);

            ni2 = ci-P2;nj2 = cj-P2;
            n2_twosides = -[ni2(2)*nj2(3)-ni2(3)*nj2(2),ni2(3)*nj2(1)-ni2(1)*nj2(3),ni2(1)*nj2(2)-ni2(2)*nj2(1)];
            n2_twosides = n2_twosides/norm(n2_twosides);

            segment1_twosides = [P1,n1_twosides,Rp,P2_segment1,A1,theta1,1;P2,n2_twosides,Rp,A1,P1_segment1,theta1,1];
            segment2_twosides = [P2,n2_twosides,Rp,P2_segment2,A2,theta2,1;P1,n1_twosides,Rp,A2,P1_segment2,theta2,1];

            segment1 = [s1;segment1_twosides];
            segment2 = [s2;segment2_twosides];

            segment1_label = [1,segment(i,1),i,0,0;2,index_I1,segment(i,1),segment(i,2),1;2,index_I2,segment(i,1),segment(i,2),1]; % [type,index_atom,index_segment], [2,index_segment,atom_i,atom_j,1 or -1]
            segment2_label = [1,segment(i,2),i,0,0;2,index_I1,segment(i,1),segment(i,2),-1;2,index_I2,segment(i,1),segment(i,2),-1]; 

            arg_eSES = ext_segment(i);
            arg_SAScircle = 0;

            theta3 = theta1+theta2;

            IJ = segment(i,1:2);
            output_SES_toroidalpatches1(sing,A1,A2,[theta1,theta2,theta3],SASsegment,segment1,segment2,segment1_label,segment2_label,[],[],arg_eSES,arg_SAScircle,IJ);
        end
    else 
        theta1=acos(r/(ri+Rp));
        theta2=acos(r/(rj+Rp));
        
        if norm(ci-cj)>max(sqrt((ri+Rp)^2-r^2),sqrt((rj+Rp)^2-r^2))
            theta=theta1+theta2;
            As=angle*Rp*(r*(theta)-Rp*(sin(theta1)+sin(theta2)));
            V_delete=angle*Rp^2*(r*(theta)/2-Rp/3*(sin(theta1)+sin(theta2)));
        else
            theta=abs(theta1-theta2);
            As=angle*Rp*(r*(theta)-Rp*(abs(sin(theta1)-sin(theta2))));
            V_delete=angle*Rp^2*(r*(theta)/2-Rp/3*(abs(sin(theta1)-sin(theta2))));
        end
        Vs_sas=angle/2*r^2*norm(ci-cj)/3;
        Vs_ses=Vs_sas-V_delete;

        if ext_segment(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_toroide(ci,cj,ri,rj,direct,n,angle,r,A,P_k1,P_k2,theta1,theta2,Rp,i);
            end
            DataAV.Acses=DataAV.Acses+As;
            DataAV.Aeses=DataAV.Aeses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            DataAV.Vesas=DataAV.Vesas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            DataAV.Veses=DataAV.Veses+Vs_ses;
        else
            DataAV.Acses=DataAV.Acses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
                mesh_toroide(ci,cj,ri,rj,direct,n,angle,r,A,P_k1,P_k2,theta1,theta2,Rp,i);
            end

        end
        
        %%%%% output the toroidal patch corresponding to the SAS-segment
        if Para.out_MolStrc
            sing = 0; 
            if direct == 1
                P1 = P_k1;
                P2 = P_k2;
                index_I1 = segment(i,3);
                index_I2 = segment(i,4);
            else
                P1 = P_k2;
                P2 = P_k1;
                index_I1 = segment(i,4);
                index_I2 = segment(i,3);
            end

            SASsegment = [A,n,r,P1,P2,angle]; % [c,n,r,s,e,angle], n points from ci to cj

            c_segment1 = ci+ri/(ri+Rp)*(A-ci);
            c_segment2 = cj+rj/(rj+Rp)*(A-cj);
            n_segment1 = -n; % the interior is on the right-hand along the segment
            n_segment2 = n;
            r_segment1 = r*ri/(ri+Rp);
            r_segment2 = r*rj/(rj+Rp);
            P1_segment1 = ci+ri/(ri+Rp)*(P2-ci);
            P2_segment1 = ci+ri/(ri+Rp)*(P1-ci);
            P1_segment2 = cj+rj/(rj+Rp)*(P1-cj);
            P2_segment2 = cj+rj/(rj+Rp)*(P2-cj);

            s1 = [c_segment1,n_segment1,r_segment1,P1_segment1,P2_segment1,angle,1]; % [c,n,r,s,e,angle,right-hand]
            s2 = [c_segment2,n_segment2,r_segment2,P1_segment2,P2_segment2,angle,1];

            ni1 = ci-P1;nj1 = cj-P1;
            n1_twosides = [ni1(2)*nj1(3)-ni1(3)*nj1(2),ni1(3)*nj1(1)-ni1(1)*nj1(3),ni1(1)*nj1(2)-ni1(2)*nj1(1)];
            n1_twosides = n1_twosides/norm(n1_twosides);

            ni2 = ci-P2;nj2 = cj-P2;
            n2_twosides = -[ni2(2)*nj2(3)-ni2(3)*nj2(2),ni2(3)*nj2(1)-ni2(1)*nj2(3),ni2(1)*nj2(2)-ni2(2)*nj2(1)];
            n2_twosides = n2_twosides/norm(n2_twosides);

            if norm(ci-cj)>max(sqrt((ri+Rp)^2-r^2),sqrt((rj+Rp)^2-r^2))
                theta3=theta1+theta2;
            else
                theta3=abs(theta1-theta2);
            end

            segment1_side = [P1,n1_twosides,Rp,P2_segment1,P1_segment2,theta3,1];
            segment2_side = [P2,n2_twosides,Rp,P2_segment2,P1_segment1,theta3,1];

            segment0 = [s1;segment1_side;s2;segment2_side];

            segment0_label = [1,segment(i,1),i,0,0;2,index_I1,segment(i,1),segment(i,2),0; 1,segment(i,2),i,0,0;2,index_I2,segment(i,1),segment(i,2),0]; % [type,index_atom,index_segment], [2,index_segment,atom_i,atom_j,1 or -1]

            arg_eSES = ext_segment(i);
            arg_SAScircle = 0;

            IJ = segment(i,1:2);
            output_SES_toroidalpatches1(sing,[],[],[theta1,theta2,theta3],SASsegment,[],[],[],[],segment0,segment0_label,arg_eSES,arg_SAScircle,IJ);
        end
    end
    
    
    
end

%% Circle case
for i=1:ncircle
    r=circle(i,9);%% the radius of the i-th circle
    A=circle(i,3:5);
    n=circle(i,6:8);
    direct=1;
    [v1,v2]=orthogonalvectors(n);
    P_k1=A+r*v1;

    ri=R(circle(i,1));
    rj=R(circle(i,2));
    ci=C(circle(i,1),:);
    cj=C(circle(i,2),:);
        
    %%
    if r<Rp && (ci-A)*(cj-A)'<0
        A1=A+sqrt(Rp^2-r^2)*n;
        A2=A-sqrt(Rp^2-r^2)*n;
        
        theta1=acos(r/(ri+Rp))-acos(r/Rp);
        theta2=acos(r/(rj+Rp))-acos(r/Rp);

        c1=ci*Rp/(Rp+ri)+A*ri/(Rp+ri);
        r1=ri/(Rp+ri)*r;
        if ext_circle(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_cusp(c1,r1,A1,direct,n,2*pi,r,A,P_k1,P_k1,theta1,Rp,-i);
            end
        else
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
                mesh_cusp(c1,r1,A1,direct,n,2*pi,r,A,P_k1,P_k1,theta1,Rp,-i);
            end
        end
        c1=cj*Rp/(Rp+rj)+A*rj/(Rp+rj);
        r1=rj/(Rp+rj)*r;
        
        theta=acos(r/Rp);
        As=2*pi*Rp*(r*(theta1+theta2)-Rp*(sin(theta1+theta)+sin(theta2+theta)-2*sin(theta)));
        V_delete=2*pi*Rp^2*(r*(theta1+theta2)/2-Rp/3*(sin(theta1+theta)+sin(theta2+theta)-2*sin(theta)))+2*pi*r^2*sqrt(Rp^2-r^2)/3;
        Vs_sas=pi*r^2*norm(ci-cj)/3;
        Vs_ses=Vs_sas-V_delete;
        
        if ext_circle(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_cusp(c1,r1,A2,direct,n,2*pi,r,A,P_k1,P_k1,theta2,Rp,-i);
            end
            DataAV.Acses=DataAV.Acses+As;
            DataAV.Aeses=DataAV.Aeses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            DataAV.Vesas=DataAV.Vesas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            DataAV.Veses=DataAV.Veses+Vs_ses;
        else
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
                mesh_cusp(c1,r1,A2,direct,n,2*pi,r,A,P_k1,P_k1,theta2,Rp,-i);
            end
            
            DataAV.Acses=DataAV.Acses+As;
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
        end
        
        %%%%% output the toroidal patch corresponding to the SAS circle
        if Para.out_MolStrc
            sing = 1;

            SAScircle = [A,n,r];

            c_circle1 = ci+ri/(ri+Rp)*(A-ci);
            c_circle2 = cj+rj/(rj+Rp)*(A-cj);
            n_circle1 = -n; % the interior in on the right-hand along the circle
            n_circle2 = n;
            r_circle1 = r*ri/(ri+Rp);
            r_circle2 = r*rj/(rj+Rp);
            circle1 = [c_circle1,n_circle1,r_circle1,1]; % [c,n,r,right-hand]
            circle2 = [c_circle2,n_circle2,r_circle2,1];
            circle1_label = [1,circle(i,1),-i]; %[type,index_atom,-index_circle]
            circle2_label = [1,circle(i,2),-i]; 

            arg_eSES = ext_circle(i);
            arg_SAScircle = 1;

            theta3 = theta1+theta2;
            IJ = circle(i,1:2);
            output_SES_toroidalpatches2(sing,A1,A2,[theta1,theta2,theta3],SAScircle,circle1,circle2,circle1_label,circle2_label,arg_eSES,arg_SAScircle,IJ);
        end
    else 
        
        theta1=acos(r/(ri+Rp));
        theta2=acos(r/(rj+Rp));
        
        if norm(ci-cj)>max(sqrt((ri+Rp)^2-r^2),sqrt((rj+Rp)^2-r^2))
            theta=theta1+theta2;
            As=2*pi*Rp*(r*(theta)-Rp*(sin(theta1)+sin(theta2)));
            V_delete=2*pi*Rp^2*(r*(theta)/2-Rp/3*(sin(theta1)+sin(theta2)));
        else
            theta=abs(theta1-theta2);
            As=2*pi*Rp*(r*(theta)-Rp*(abs(sin(theta1)-sin(theta2))));
            V_delete=2*pi*Rp^2*(r*(theta)/2-Rp/3*(abs(sin(theta1)-sin(theta2))));
        end
        Vs_sas=pi*r^2*norm(ci-cj)/3;
        Vs_ses=Vs_sas-V_delete;
                
        if ext_circle(i)==1
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_ext && Para.viz_SES(1)
                    figure(Figs.ext)
                    hold on;
                end
                mesh_toroide(ci,cj,ri,rj,direct,n,2*pi,r,A,P_k1,P_k1,theta1,theta2,Rp,-i);
            end
            DataAV.Acses=DataAV.Acses+As;
            DataAV.Aeses=DataAV.Aeses+As;
            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;
            DataAV.Vesas=DataAV.Vesas+Vs_sas;
            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
            DataAV.Veses=DataAV.Veses+Vs_ses;
        else
            if Para.arg_meshing(1)
                if Para.arg_viz && Para.arg_int && Para.viz_SES(1)
                    figure(Figs.int)
                    hold on;
                end
                mesh_toroide(ci,cj,ri,rj,direct,n,2*pi,r,A,P_k1,P_k1,theta1,theta2,Rp,-i);
            end
            DataAV.Acses=DataAV.Acses+As;            
            DataAV.Vcsas=DataAV.Vcsas+Vs_sas;            
            DataAV.Vcses=DataAV.Vcses+Vs_ses;
        end
        
        %%%%% output the toroidal patch corresponding to the SAS circle
        if Para.out_MolStrc
            sing = 0;
            A1 = [];
            A2 = [];

            SAScircle = [A,n,r];

            c_circle1 = ci+ri/(ri+Rp)*(A-ci);
            c_circle2 = cj+rj/(rj+Rp)*(A-cj);
            n_circle1 = -n; % the interior in on the right-hand along the circle
            n_circle2 = n;
            r_circle1 = r*ri/(ri+Rp);
            r_circle2 = r*rj/(rj+Rp);
            circle1 = [c_circle1,n_circle1,r_circle1,1]; % [c,n,r,right-hand]
            circle2 = [c_circle2,n_circle2,r_circle2,1];
            circle1_label = [1,circle(i,1),-i]; %[type,index_atom,-index_circle]
            circle2_label = [1,circle(i,2),-i]; 

            arg_eSES = ext_circle(i);
            arg_SAScircle = 1;

            if norm(ci-cj)>max(sqrt((ri+Rp)^2-r^2),sqrt((rj+Rp)^2-r^2))
                theta3=theta1+theta2;
            else
                theta3=abs(theta1-theta2);
            end

            IJ = circle(i,1:2);
            output_SES_toroidalpatches2(sing,A1,A2,[theta1,theta2,theta3],SAScircle,circle1,circle2,circle1_label,circle2_label,arg_eSES,arg_SAScircle,IJ);
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mesh a thin toroide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh_toroide(ci,cj,ri,rj,direct,n,angle,r,A,P_k1,P_k2,theta1,theta2,Rp,k)
global Para;
d = Para.arg_meshing(2);
if norm(ci-cj)>max(sqrt((ri+Rp)^2-r^2),sqrt((rj+Rp)^2-r^2))
    flag=1;
    theta=theta1+theta2;
else
    flag=0;
    theta=abs(theta1-theta2);
end

theta0=pi/3;
% N_probe=real(floor(Rp*(theta)/d)+1);
N_probe=real(max(floor(Rp*(theta)/d)+1,floor(theta/theta0)+1));%% to gurantee the length of each edge is greater than tolerance
N_arc=real(max(floor(angle/theta0)+1,floor(r*angle*max(ri,rj)/((Rp+max(ri,rj))*d))+1));%% length constraint and (angle constraint, impossible)

N_point=(N_probe+1)*(N_arc+1);%% the number of points
P=zeros(N_point,3);%% record all the coordinates of points
Np=0;

points1=zeros(N_probe+1,3);
points2=zeros(N_probe+1,3);

for i=0:N_probe
   for j=0:N_arc
       Np=Np+1;
       P(Np,:)=coord_toroide(i,j,N_probe,N_arc,P_k1,r,A,Rp,theta1,theta2,angle,direct,n,ci,flag);
   end
   points1(i+1,:)=P(Np-N_arc,:);
   points2(i+1,:)=P(Np,:);
end

N_triangle=2*N_probe*N_arc;
T=zeros(N_triangle,3);
Nt=0;
for i=1:N_probe
   for j=1:N_arc
       i1=i-1;j1=j;
       i2=i;j2=j-1;
       i3=i;j3=j;

       Nt=Nt+1;
       [k1,k2,k3]=index_P_toroide(i1,j1,i2,j2,i3,j3,N_arc);
       T(Nt,:)=[k1,k2,k3];

       i1=i-1;j1=j;
       i2=i;j2=j-1;
       i3=i-1;j3=j-1;

       Nt=Nt+1;
       [k1,k2,k3]=index_P_toroide(i1,j1,i2,j2,i3,j3,N_arc);
       T(Nt,:)=[k1,k2,k3];       
   end
end

normal = zeros(Nt,3);
for i = 1:Nt
    p = (P(T(i,1),:)+P(T(i,2),:)+P(T(i,3),:))/3;
    nt = (p-A)-((p-A)*n')*n;
    xp = A+r*nt/norm(nt);
    normal(i,:) = (xp-p)/norm(xp-p);
end

if Para.arg_viz && Para.viz_SES(1)
    visutorpat(P,Np,T,Nt,normal);
end

global Points_SES;
Points_SES = [Points_SES;P];

end

function [k1,k2,k3]=index_P_toroide(i1,j1,i2,j2,i3,j3,N_arc)
k1=i1*(N_arc+1)+j1+1;
k2=i2*(N_arc+1)+j2+1;
k3=i3*(N_arc+1)+j3+1;
end
%%return the coordinate of the (i,j) point on the toroide
function x=coord_toroide(i,j,N_probe,N_arc,P_k1,r,A,Rp,theta1,theta2,angle,direct,n,ci,flag)

if flag==1
    theta_i=i*(theta1+theta2)/N_probe;
    d=abs(norm(ci-A)-tan(theta1-theta_i)*r);
else
    theta_i=i*(theta1-theta2)/N_probe;
    d=abs(norm(ci-A)-tan(theta1-theta_i)*r);
end
Ai=ci+d*n;

u=(P_k1-A)/r;
v=direct*[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];
angle_j=j/N_arc*angle;
P_j=r*cos(angle_j)*u+r*sin(angle_j)*v+A;

%%%%%%
x=P_j+Rp*(Ai-P_j)/norm(Ai-P_j);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mesh a cusp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh_cusp(c1,r1,A1,direct,n,angle,r,A,P_k1,P_k2,theta1,Rp,k)
global Para;
d = Para.arg_meshing(2);

N_probe=real(floor(Rp*theta1/d)+1);% corresponds to the arc on the probe at fixed point, to gurantee the length of each edge is greater than tolerance
theta0=pi/3;
N_arc=real(floor(angle/theta0)+1); % corresponds to the circular arc
%N_arc=1;

N_point=1;%%%the number of points

t=zeros(N_probe+1,1);
for i=1:N_probe
   theta_i=acos(r/Rp)+theta1/N_probe*i;
   L_i=r/cos(theta_i);
   r_i=(L_i-Rp)/L_i*r;
   angle_division=angle/N_arc;
   t(i+1)=floor(r_i*angle_division/d)+1;
   N_point=N_point+t(i+1)*N_arc+1;
end
P=zeros(real(N_point),3);%% record all the coordinates of points
Np=0;

points1=zeros(N_probe+1,3);
points2=zeros(N_probe+1,3);

for i=0:N_probe
   for j=0:(t(i+1)*N_arc)
       Np=Np+1;
       P(Np,:)=coord_cusp(i,j,N_probe,N_arc,P_k1,A1,r,A,Rp,theta1,angle,direct,n,t);
       if j==0
           points1(i+1,:)=P(Np,:);
       end
       if j==t(i+1)*N_arc
           points2(i+1,:)=P(Np,:);
       end

   end
end
points0=P((N_point-t(N_probe+1)*N_arc):N_point,:);

N_triangle=N_probe^2*N_arc;
T=zeros(N_triangle,3);
Nt=0;
for i=1:N_probe
   for j=1:(t(i+1)*N_arc)
       if i<N_probe
           if t(i+1)==t(i)+1
               i1=i-1;j1=j-floor((j-1)/t(i+1))-1;
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];
           else
               i1=i-1;j1=j;
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];
           end
           
           if t(i+1)==t(i+2)-1
               i1=i+1;j1=j+floor((j-1)/t(i+1));
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];    
           else
               i1=i+1;j1=j-1;
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];
           end
           
       else
           
           if t(i+1)==t(i)+1
               i1=i-1;j1=j-floor((j-1)/t(i+1))-1;
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];
           else
               i1=i-1;j1=j;
               i2=i;j2=j-1;
               i3=i;j3=j;

               Nt=Nt+1;
               [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t);
               T(Nt,:)=[k1,k2,k3];
           end
           
       end
   end
end

normal = zeros(Nt,3);
for i = 1:Nt
    p = (P(T(i,1),:)+P(T(i,2),:)+P(T(i,3),:))/3;
    nt = (p-A)-((p-A)*n'/norm(n))*n/norm(n);
    xp = A+r*nt/norm(nt);
    
    n0 = (xp-p)/norm(xp-p);
    
    n_1 = P(T(i,2),:)-P(T(i,1),:);
    n_2 = P(T(i,3),:)-P(T(i,1),:);
    n = [n_1(2)*n_2(3)-n_1(3)*n_2(2),n_1(3)*n_2(1)-n_1(1)*n_2(3),n_1(1)*n_2(2)-n_1(2)*n_2(1)];
    n = n/norm(n);
    if n*n0'<0
        n = -n;
    end
    
    normal(i,:) = (xp-p)/norm(xp-p);
end

if Para.arg_viz && Para.viz_SES(1)
    visutorpat(P,Np,T,Nt,normal);
end

end

function [k1,k2,k3]=index_P(i1,j1,i2,j2,i3,j3,N_arc,t)%%k1,k2,k3 denote the three vertex of the triangle
k1=0;k2=0;k3=0;

for i=0:(i1-1)
    k1=k1+t(i+1)*N_arc+1;
end
k1=k1+j1+1;

for i=0:(i2-1)
    k2=k2+t(i+1)*N_arc+1;
end
k2=k2+j2+1;

for i=0:(i3-1)
    k3=k3+t(i+1)*N_arc+1;
end
k3=k3+j3+1;

end

%% return the coordinate of the (i,j) point on the cusp
function x=coord_cusp(i,j,N_probe,N_arc,P_k1,A1,r,A,Rp,theta1,angle,direct,n,t)

if i==0
    x=A1;
    return;
end

theta=acos(r/Rp);
radius_i=(r/cos(theta+theta1*i/N_probe)-Rp)/(r/cos(theta+theta1*i/N_probe))*r;
Ai=A+sqrt((r/cos(theta+theta1*i/N_probe))^2-r^2)*(A1-A)/norm(A1-A);

division=t(i+1)*N_arc;

u=(P_k1-A)/r;
v=direct*[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];
angle_j=j/division*angle;
P_j=r*cos(angle_j)*u+r*sin(angle_j)*v+A;

x=radius_i/r*P_j+(r-radius_i)/r*Ai;
end

function P_mid = midpt(P_k1,A,r,n,direct,angle_mid)

u=(P_k1-A)/r;
v=direct*[n(2)*u(3)-n(3)*u(2),n(3)*u(1)-n(1)*u(3),n(1)*u(2)-n(2)*u(1)];
P_mid=r*cos(angle_mid)*u+r*sin(angle_mid)*v+A;

end
