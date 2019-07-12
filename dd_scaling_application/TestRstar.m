function [F,names,dF,dvdW] = TestRstar

names = ["1B17","1CRN","1ETN","1GZI","1YJO","3DIK","101M","3DIK","carbo","glutaredoxin","hiv-1-gp41","l-plectasin","ubch5b","vancomycin"];
names = "1B17";
format = ".pdb";
% names =
% ["carbo.xyz","glutaredoxin","hiv-1-gp41","l-plectasin","ubch5b","vancomycin"];
np = 10;

for i=1:length(names)
    char(names(i))
    
    file.name = char(names(i));
    file.format = char(format(i));
    
    f = FindRstar(np,file);
    s = size(f);
    F(1:s(1),1:s(2),i) = f;
    dvdW(i) = f(1,1);
    dF(i) = f(end,1);
    F(end,1,i)
    save('tmp_save')
end