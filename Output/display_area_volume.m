%%  Display Molecular Areas and Volumes
function display_area_volume
% c
% c   output the areas and volumes of (complete or exterior) SES and SAS
% c

format long;%shortE;
global DataAV;

disp('-------------------------------------');
disp('! Table of Different Areas & Volumes');
disp('-------------------------------------');

L = {'External';'Complete';'Difference'};
Asas = [DataAV.Aesas;DataAV.Acsas;DataAV.Aesas-DataAV.Acsas];
Vsas = [DataAV.Vesas;DataAV.Vcsas;DataAV.Vesas-DataAV.Vcsas];
Ases = [DataAV.Aeses;DataAV.Acses;DataAV.Aeses-DataAV.Acses];
Vses = [DataAV.Veses;DataAV.Vcses;DataAV.Veses-DataAV.Vcses];

AreaVol = table(Asas,Vsas,Ases,Vses,...
    'RowNames',L);
disp(AreaVol);

end
