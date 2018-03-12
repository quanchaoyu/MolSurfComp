function CompAreaVol_concave(c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize)
global Para;
Rp = Para.radius_probe;
%% Area & Volume computation
global DataAV;
global arg_eSAS;

global V_cSES;
global V_eSES;

if arg_eSAS==1
    Area=area_spherical(c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize);
    DataAV.Aeses=DataAV.Aeses+Area;
    V_eSES=V_eSES-Area*Rp/3;
elseif arg_eSAS==0
    Area=area_spherical(c_sphere,r_sphere,loops,loopsize,segment,circle,patches,patchesize);
    DataAV.Acses=DataAV.Acses+Area;
    V_cSES=V_cSES-Area*Rp/3;
end

end