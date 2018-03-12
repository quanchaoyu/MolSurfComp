function viz_MSMS(filename)
figure('Color',[1,1,1]);
axis off;
view([40,25]);
light('position',[1,1,1]);
lighting gouraud;
axis equal;

InitCode; %% add pathes

filename1 = [filename,'.face']; %filenames
filename2 = [filename,'.vert']; 

F = textread(filename1);
V = textread(filename2);

P = V(:,1:3);
T = F(:,1:3);
Nt = size(T,1);

T1 = zeros(Nt,3);
T2 = zeros(Nt,3);
T3 = zeros(Nt,3);
Nt1 = 0;
Nt2 = 0;
Nt3 = 0;

for i = 1:Nt
    if F(i,4) == 1
        Nt1 = Nt1+1;
        T1(Nt1,:) = T(i,:);
    elseif F(i,4) == 2
        Nt2 = Nt2+1;
        T2(Nt2,:) = T(i,:);
    elseif F(i,4) == 3
        Nt3 = Nt3+1;
        T3(Nt3,:) = T(i,:);
    end
end
T1 = T1(1:Nt1,:);
T2 = T2(1:Nt2,:);
T3 = T3(1:Nt3,:);

visu(P,T1,Nt1,'r');
visu(P,T2,Nt2,'b');
visu(P,T3,Nt3,'y');
end