function Geom = read_xyzr(filename)
% c
% c   Read xyzr-file and save data to Geom
% c

filename = [filename,'.xyzr'];
Mat = dlmread(filename);
Geom.M = size(Mat,1);
Geom.centers = Mat(:,1:3);
Geom.R = Mat(:,4);

end
