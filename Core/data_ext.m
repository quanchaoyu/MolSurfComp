function data_ext

DataGlob; % call Para, Geom, Inter, data of I Cir Seg Loop Pat, AreaVolume

global Ext;
%% neighborship between patches
global neighbor;
global Nneighbor;
neighbor=zeros(npatches,max(Row)); %  the j-th row records all the neighbor patches to the j-th spherical patch 
Nneighbor=zeros(npatches,1); %  the j-th row records all the neighbor patches number to the j-th spherical patch

% for i=1:M
%    if patches_index(i,1)>0 %  i denotes the i-th sphere
%        for j=patches_index(i,1):patches_index(i,2) 
%             count = 0;
%             for k=1:patchesize(j,1) %  k denotes the k-th loop on the boundary of the j-th spherical patch
%                if patches(j,k)>0 %  patches(j,k)>0 denotes the index of the loop on the i-th sphere
% 
%                    ln=loops_index(i,1)+patches(j,k)-1; %  the loop number
%                    count = count+loopsize(ln);
% 
%                else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere
%                     count = count+1;
%                end
%             end
%             if Nneighbor(j) - count
%                 Nneighbor(j) - count
%             end
%        end
%    end
% end


for i=1:M
   if patches_index(i,1)>0 %  i denotes the i-th sphere
       for j=patches_index(i,1):patches_index(i,2) %% j denotes the j-th spherical patch
           for k=1:patchesize(j,1) %  k denotes the k-th loop on the boundary of the j-th spherical patch
               if patches(j,k)>0 %  patches(j,k)>0 denotes the index of the loop on the i-th sphere
                   
                   ln=loops_index(i,1)+patches(j,k)-1; %  the loop number
                   for s0=1:loopsize(ln)
                       sn=loops(ln,s0); %  the segment number
                       j0=neighbor_patch(sn,seg(sn,1:2),i,patches,patchesize,patches_index,loops,loopsize,loops_index,circleindex);
                       Nneighbor(j)=Nneighbor(j)+1;
                       neighbor(j,Nneighbor(j))=j0;
                   end
                   
               else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere
                   
                   cn=-circleindex(i,-patches(j,k)); %  the circle number
                   j0=neighbor_patch(cn,circle(-cn,1:2),i,patches,patchesize,patches_index,loops,loopsize,loops_index,circleindex);
                   Nneighbor(j)=Nneighbor(j)+1;
                   neighbor(j,Nneighbor(j))=j0;
               end
           end
       end
   end
end

%% ext_patch, ext_segment, ext_I

%c  obtain an initial patch outside
[~,i_left]=min(C(1:M,1)-R); %  the leftmost sphere indice
t=1;
if patches_index(i_left,1) == 0
    [~,i_left]=min(C(1:M,2)-R);
end

for j=patches_index(i_left,1):patches_index(i_left,2)
    
    if j==patches_index(i_left,1)
        
        j_left=j;
        
        if j == 0
            display('Error: there exists isolated SAS-ball!');
        end
        
        for k=1:patchesize(j,1)
            if patches(j,k)>0 %  patches(j,k)>0 denotes the index of the loop on the i-th sphere
                   
               ln=loops_index(i_left,1)+patches(j,k)-1; %  the loop number
               for s0=1:loopsize(ln)
                   sn=loops(ln,s0); %  the segment number
                   if t==1
                       x_left=I(seg(sn,3),1);
                       t=0;
                   elseif x_left>I(seg(sn,3),1)
                       x_left=I(seg(sn,3),1);
                   end
               end

           else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere

               cn=-circleindex(i_left,-patches(j,k)); %  the circle number
               if t==1
                   x_left=circle(-cn,3);
                   t=0;
               elseif x_left>circle(-cn,3)
                   x_left=circle(-cn,3);
               end
            end
        end
        
    else
        
        for k=1:patchesize(j,1)
            if patches(j,k)>0 %  patches(j,k)>0 denotes the index of the loop on the i-th sphere
                   
               ln=loops_index(i_left,1)+patches(j,k)-1; %% the loop number
               for s0=1:loopsize(ln)
                   sn=loops(ln,s0); %  the segment number
                   if x_left>I(seg(sn,3),1)
                       j_left=j;
                       x_left=I(seg(sn,3),1);
                   end
               end

           else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere

               cn=-circleindex(i_left,-patches(j,k)); %  the circle number
               if x_left>circle(-cn,3)
                   j_left=j;
                   x_left=circle(-cn,3);
               end
            end
        end
        
    end
end

ext_patch=zeros(npatches,1); %  ext_patch(i) denote if the i-th spherical patch is on the eSAS
ext_patch(j_left)=1;

ext_patchset=zeros(npatches,1); %  ext_patchset denotes the set of spherical patches on the eSAS
ext_patchset(1)=j_left;
nset=1;

ext_segment=zeros(nsegment,1); %  ext_segment(i) denotes if the i-th segment is on the eSAS
ext_circle=zeros(ncircle,1);

for i=1:npatches
    if i<=nset
        j=ext_patchset(i);
        
        for k=1:Nneighbor(j)
            n=neighbor(j,k); %  n-th patch is a neighbor to the j-th patch
            if ext_patch(n)==0 %  in this case, the n-th patch is a new neighbor to the j-th patch
                nset=nset+1;
                ext_patchset(nset)=n;
                ext_patch(n)=1;
            end
        end
    else
        break;
    end
end


for i=1:nset
   j=ext_patchset(i); %  the j-th spherical patch
   for k=1:patchesize(j,1)
       
        if patches(j,k)>0 %  patches(j,k)>0 denotes the index of loop on the i-th sphere

           ln=loops_index(patch_atom(j),1)+patches(j,k)-1; %  the loop number
           for s0=1:loopsize(ln)
               sn=loops(ln,s0); %  the segment number
               if ext_segment(sn)==0
                   ext_segment(sn)=1;
               end
           end

       else %  patches(j,k)<0 denotes the index of the circle on the i-th sphere

           cn=-circleindex(patch_atom(j),-patches(j,k)); %% the circle number
           if ext_circle(-cn)==0
               ext_circle(-cn)=1;
           end
        end
    end
end

ext_I=zeros(s,1); %  ext_I denotes if the i-th intersection point is on the eSAS
for i=1:nsegment
   if ext_segment(i)==1 
       j=seg(i,3);
       if ext_I(j)==0
           ext_I(j)=1;
       end
       j=seg(i,4);
       if ext_I(j)==0
           ext_I(j)=1;
       end
   end
end

Ext.I = ext_I;
Ext.circle = ext_circle;
Ext.segment = ext_segment;
Ext.patch = ext_patch;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return the neighbor patch of the j-th patch w.r.t the n-th segment (or circle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function j = neighbor_patch(n,s,i,patches,patchesize,patches_index,loops,loopsize,loops_index,circleindex)
if s(1) == i %  Assume that [i,j] are the two corresponding spheres
    i0 = s(2);
elseif s(2) == i
    i0 = s(1);
else
    printf('error');
end

for j = patches_index(i0,1):patches_index(i0,2)
    for k = 1:patchesize(j,1)
        if patches(j,k)>0 && n>0 %%% patches(j,k) denote the k-th loop on the j-th spherical patch
            ln = loops_index(i0,1)+patches(j,k)-1; %% t lheoop number
            for t = 1:loopsize(ln)
                if loops(ln,t) == n
                    return;%% return j
                end
            end
        elseif patches(j,k)<0 && n<0 
            cn = -circleindex(i0,-patches(j,k)); %% the circle number
            if cn == n
                return;
            end
        end
    end
end

if j == patches_index(i0,2)
    
end

end

