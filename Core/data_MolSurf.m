%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_MolSurf
%c
%c  compute the construction of SAS, including intersection points, circular arcs (circles) and spherical patches
%c
global Para;
global Figs;

if(Para.arg_viz) && (Para.arg_ext)
    Figs.ext = figure('Position',[0 600 600 500],'Name','External SAS & SES');
    FigSetting;
end

if (Para.arg_viz) && (Para.arg_int)
    Figs.int = figure('Position',[0 0 600 500],'Name','Internal SES');
    FigSetting;
end

%% SAS structure
data_I_Cir; % compute intersection points and intersection circles
data_Seg_Pat; % compute circular segments and spherical patches
data_ext; % compute data of exterior SAS

if(Para.arg_viz) && (Para.arg_ext)
    figure(Figs.ext);
    visu_extsegment;
    %FigSetting;
end

%% SES structure
data_SESsphpat; % data of SES spherical patches
data_SEStorpat; % data of SES toroidal patches

end




