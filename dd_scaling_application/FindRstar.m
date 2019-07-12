function f = FindRstar(np,filename)


Rp = 6;
F = FsesBox(Rp,np,filename) - Rp;
DistvdW = max(max(max(F)));
f(1,1) = 0;
f(1,2) = DistvdW;
if(DistvdW<0)
    disp('Problem: vdw-distance is negative')
    return
end

Rp = 3;
cnt = 2;
f(cnt,1) = Rp;
F = FsesBox(Rp,np,filename) - Rp;
f(cnt,2) = max(max(max(F)));

while(f(cnt,2)>=0)
    disp('Initial Rp not large enough, retry')
%     A = [[f(1,1) 1];[f(2,1) 1]];
%     fr = [f(1,2);f(2,2)];
%     c = A\fr;
%     Rpm1 = Rp;
%     Rp = -c(2)/c(1);
%     if(Rp<Rpm1)
    Rp = 2*Rp;
%     end
%    cnt = cnt+1;
    f(cnt,1) = Rp;
    F = FsesBox(Rp,np,filename) - Rp;
    f(cnt,2) = max(max(max(F)));

    f
end
a = f(1,:);
b = f(2,:);
Rpm1 = Rp;

for i=1:10
    f
    cnt = cnt+1;
    A = [[a(1,1) 1];[b(1,1) 1]];
    fr = [a(1,2);b(1,2)];
    c = A\fr;
    Rpm1 = Rp;
    Rp = -c(2)/c(1);
    f(cnt,1) = Rp;
    tic
    F = FsesBox(Rp,np,filename) - Rp;
    f(cnt,2) = max(max(max(F)));
    toc
    if(abs(Rpm1-Rp)/Rp<1e-2)
        break;
    end
    if(f(cnt,2)<0)
        b = [Rp,f(cnt,2)];
    else
        a = [Rp,f(cnt,2)];
    end
end
f


