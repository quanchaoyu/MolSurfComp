%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% data structure of SES spherical patches %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_SESsphpat

DataGlob; % call Para, Geom, Inter, data of I Cir Seg Loop Pat, AreaVolume

%% convex spherical patches
for i=1:npatches
   if ext_patch(i)==1 %  in this case, the spherical patch is on the CSAS
       i0=patch_atom(i); %  i0 is the atom where the spherical patch is
       
       %c  mesh spherical SAS-patches
       %mesh_sphpat(C(i0,:),R(i0)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j));
       
       DataAV.Aesas=DataAV.Aesas+Area_sphpat(i);
       DataAV.Aeses=DataAV.Aeses+(R(i0)/(R(i0)+Rp))^2*Area_sphpat(i);
       
       DataAV.Vesas=DataAV.Vesas+Area_sphpat(i)*(R(i0)+Rp)/3;
       DataAV.Veses=DataAV.Veses+(R(i0)/(R(i0)+Rp))^2*Area_sphpat(i)*R(i0)/3;
   end
end

% visulize the convex patches
global arg_eSAS;
arg_eSAS = 1;

index_convexpatch = 0;
index_SASpatch = 0;

global Para Figs;

global Rj; %
Rj=zeros(nsegment,1);
segment0=zeros(nsegment,11);
for i=1:M
    if nsatom(i)>0
        nloops=loops_index(i,1)-1;
        nloops_i=loops_index(i,2)-loops_index(i,1)+1;
        loops_i=loops(nloops+1:nloops+nloops_i,:);
        loopsize_i=loopsize(nloops+1:nloops+nloops_i,1);
    elseif ncircleindex(i)>0
        loops_i=[];
        loopsize_i=[];
        nloops_i=0;
    end
    
    if nsatom(i)>0||ncircleindex(i)>0
  
        numpatches=patches_index(i,1)-1;
        npatches_i=patches_index(i,2)-patches_index(i,1)+1;
        patches_i=patches(numpatches+1:numpatches+npatches_i,:);
        patchesize_i=patchesize(numpatches+1:numpatches+npatches_i,1);

        %   modify loops_i,segment,circle so that mesh_sphpat can be used
        [loops_i0,segment0,circle0] = mod_seg_loop_cir(i,nloops_i,loops_i,loopsize_i,circleindex,ncircleindex,satom,nsatom,C,R,seg,nsegment,ncrasegment,segment0);
        
        for j=1:npatches_i
                            
            index_convexpatch = index_convexpatch+1;
            extSES = ext_patch(patches_index(i,1)+j-1);
            
            if Para.out_MolStrc
                output_SES_convexpatches(i,index_convexpatch,extSES,C(i,:),R(i)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j),circleindex);
            end
            
            index_SASpatch = index_SASpatch+1;
            extSAS = extSES;
            
            if Para.out_MolStrc
                output_SAS_patches(i,index_SASpatch,extSAS,C(i,:),R(i)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j),circleindex);
            end
            
            % mesh and visulize each convex patch
            if ext_patch(patches_index(i,1)+j-1)==1
                if Para.arg_meshing(1) 
                    if Para.arg_viz && Para.arg_ext
                        figure(Figs.ext)
                        hold on;
                        mesh_sphpat(C(i,:),R(i)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j));
                    end
                end
                
            else
                if Para.arg_meshing(1) 
                    if Para.arg_viz && Para.arg_int
                        figure(Figs.int)
                        hold on;
                        mesh_sphpat(C(i,:),R(i)+Rp,loops_i0,loopsize_i,segment0,circle0,patches_i(j,:),patchesize_i(j));
                    end
                end
                
            end
            
        end
    end
end


%%  construct, mesh, visulize concave patches
SESconcavepat(s,I,Iijk,Ii,In,C,R,Rp,ext_I,high_I,direction,Inter);

end
