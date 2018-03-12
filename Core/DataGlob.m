global Para Geom Inter;
M = Geom.M;
C = Geom.centers;
R = Geom.R;
Rp = Para.radius_probe;

inter = Inter.M_int; %  intersection matrix
Row = Inter.num_int; %  intersection number for each SAS-ball

global DataI DataSeg DataCir DataLoop DataPat;

I = DataI.I;
s = DataI.nI;
Ii = DataI.Ii;
Iijk = DataI.Iijk;
In = DataI.In;
direction = DataI.direction;
I_circle = DataI.I_circle;
I_circle_num = DataI.I_circle_num;
high_I = DataI.high_I;
hightvalue = DataI.hightvalue;

circle = DataCir.circle;
ncircle = DataCir.ncircle;
circleindex = DataCir.circleindex;
ncircleindex = DataCir.ncircleindex;

seg = DataSeg.segment;
nsegment = DataSeg.nsegment;
ncrasegment = DataSeg.ncrasegment;
satom = DataSeg.satom;
nsatom = DataSeg.nsatom;

loops = DataLoop.loops;
nloops = DataLoop.nloops;
loopsize = DataLoop.loopsize;
loops_index = DataLoop.loops_index;

patches = DataPat.patches;
npatches = DataPat.npatches;
patchesize = DataPat.patchesize;
patches_index = DataPat.patches_index;
patch_atom = DataPat.patch_atom; % patch index to atom index

global DataAV;
Area_sphpat = DataAV.Area_sphpat;

global Ext;
if ~isempty(Ext)
    ext_I = Ext.I;
    ext_circle = Ext.circle;
    ext_segment = Ext.segment;
    ext_patch = Ext.patch;
end