function writeSTL(filename)

global Triangle;
global Vertices;
global NormalVects;

if isempty(Triangle)
    disp('No mesh generated!');
    return;
end

T = Triangle(:,1:3);
x = abs(floor((min(min(Vertices)))))+1;
V = Vertices(:,1:3)+x;
NV = NormalVects;

ID = ['./',filename,'.stl'];

N = size(T,1);

fileID = fopen(ID,'wt');
    fprintf(fileID,'%-16s \n',['solid ',filename]);
    fprintf(fileID,'\n');
    
    for i = 1:N
        v1 = V(T(i,1),:);
        v2 = V(T(i,2),:);
        v3 = V(T(i,3),:);
        
        fprintf(fileID,'%13s','facet normal '); fprintf(fileID,'%d %d %d \n',NV(i,:));
        fprintf(fileID,'%10s \n','outer loop');
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v1);
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v2);
        fprintf(fileID,'%s','vertex '); fprintf(fileID,'%d %d %d \n',v3);
        fprintf(fileID,'%7s \n','endloop');
        fprintf(fileID,'%8s \n','endfacet');
        
        fprintf(fileID,'\n');
        
    end
    
    fprintf(fileID,'%8s \n','endsolid');
fclose(fileID);

%{
for i = 1:N
    v1 = Vertices(T(i,1),1:3);
    v2 = Vertices(T(i,2),1:3);
    v3 = Vertices(T(i,3),1:3);

    line(v1(1)+[0:0.01:1]*NV(i,1),v1(2)+[0:0.01:1]*NV(i,2),v1(3)+[0:0.01:1]*NV(i,3))
    hold on;
end
%}
    
end