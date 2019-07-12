%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MolSurfComp
% c
% c
% c   Molecular Surface Computation: SES-singularitied are completely
% c   characterized!
% c   
% c   Copyright: Chaoyu QUAN & Benjamin STAMM, March 2018.
% c   Written by Chaoyu QUAN.
% c   Email: quanchaoyu@gmail.com (C. QUAN) or bsberkeley@gmail.com (B. STAMM)
% c
% c   References: 
% c   [1] C. Quan, B. Stamm, Mathematical analysis and calculation of molecular surfaces.
% c   Journal of Computational Physics, 2016, 322: 760-782.
% c
% c   [2] C. Quan, B. Stamm, Meshing molecular surfaces based on analytical implicit
% c   representation. Journal of Molecular Graphics and Modelling, 2017, 71: 200-210.
% c
% c
% c
% c   Geom: data structure of molecule
% c         Geom.M -- number of atoms
% c         Geom.R -- Van der Waals radii
% c         Geom.centers -- center coordinates
% c   
% c   Inter: for each SAS-ball, record its intersected SAS-balls
% c          
% c   FUNCTION data_MolSurf: compute the data structures of SAS components & 
% c   SES patches, then mesh molecular surfaces.
% c



%% Initialization
clearvars -global % clear all
close all

time_init = cputime; % initial time

global Para Geom Inter Figs; % parameters, geometrical structure, intersection info, figures

InitCode; % add pathes

%% Read input file
read_inputs; % read inputs

if strcmp(Para.format, '.pdb')
    Geom = read_PDB(Para.filename); % read a pdb-file, UFF VDW-radii
elseif strcmp(Para.format,'.xyzr') 
    Geom = read_xyzr(Para.filename); % 4 colums: center (x,y,z) & radius r.
elseif strcmp(Para.format,'.xyz')
    Geom = read_xyz(Para.filename);
    %Geom.R = 1.1*Geom.R;
else
    disp('It is not a readable file!');
end

if Para.radius_probe == 0
    disp('You are computing a VDW surface.');
end

%% Compute intersectionship
Inter = interstructure(Geom,Para.radius_probe); % compute the SAS balls intersected with a SAS ball

%% Compute data structure of mol surfs, mesh the surface
global Vertices Triangle NormalVects;
if Para.out_STL
    Vertices = [];
    Triangle = [];
    NormalVects = [];
end

data_MolSurf;

%% Fill inner holes
Para.arg_fillholes = false;
Para.viz_holes = false;
if Para.arg_fillholes
    fillhole;
end

%% Display molecular areas and volumes
if Para.out_AreaVol
    display_area_volume; % display a table of molecular areas and volumes
end

%% Viz points, circles & segments on the SAS
if Para.viz_CirSeg(1) || Para.viz_CirSeg(2) 
    Figs.CirSeg = figure('Position',[700 600 600 500],'Name','Intersection Points, Circles & Segments on the SAS');
    if Para.viz_CirSeg(1)
        visu_extsegment;
    end
    if Para.viz_CirSeg(2) 
        visu_innersegment;
    end
    FigSetting;
end

%% Output STL file
if (Para.out_STL)
    writeSTL(Para.filename);
end

%% Display the total run time
RunTime_total = cputime-time_init;
disp(['Total Run Time: ',num2str(RunTime_total),' seconds'])

end

