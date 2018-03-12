function [loops_i0,segment0,circle0] = mod_seg_loop_cir(i,nloops_i,loops_i,loopsize_i,circleindex,ncircleindex,satom,nsatom,C,R,segment,nsegment,ncrasegment)
%c
%c   Modify segment, loops, circles --> segment0, loops0, circles0
%c
global DataI DataCir;
I = DataI.I;
circle = DataCir.circle;

global Rj;

Rj=zeros(nsegment,1);
segment0=zeros(nsegment,11);
loops_i0=loops_i;
circle0=zeros(ncircleindex(i),7);

% segment:[i,j,p1,p2,direct],ncrasegment  -->  segment0:[c,n,r,spoint,angle]
% circle:[i,j,c,n,r]  --> [c,n,r]
% satom(i,:),nsatom(i)
% modify the segment's direction

for k=1:nsatom(i)
    direct_s=segment(satom(i,k),5)*(2*(segment(satom(i,k),1)==i)-1);%%%direct_s=-1 mean the right hand of the segment is interior
    direct_n=2*(segment(satom(i,k),1)==i)-1;%%direct_n=-1 means the normal points inside
    if direct_s==-1 %%%modify the direction of the loop, to ensure the righthand is interior
        segment0(satom(i,k),:)=[ncrasegment(satom(i,k),4:6),direct_n*ncrasegment(satom(i,k),1:3),ncrasegment(satom(i,k),7),I(segment(satom(i,k),4),:),ncrasegment(satom(i,k),8)];
    else
        segment0(satom(i,k),:)=[ncrasegment(satom(i,k),4:6),direct_n*ncrasegment(satom(i,k),1:3),ncrasegment(satom(i,k),7),I(segment(satom(i,k),3),:),ncrasegment(satom(i,k),8)];
    end

    if segment(satom(i,k),1)==i
        Rj(satom(i,k))=R(segment(satom(i,k),2));
    else
        Rj(satom(i,k))=R(segment(satom(i,k),1));
    end
end

%%%modify the loop's direction
for k=1:nloops_i
   loops_i0(k,1:loopsize_i(k))=flipud(loops_i(k,1:loopsize_i(k))')';
end

%%%modify the circle's direction
for k=1:ncircleindex(i)
    if circle(circleindex(i,k),1)==i
        circle0(k,:)=circle(circleindex(i,k),3:9);
    else
        circle0(k,:)=[circle(circleindex(i,k),3:5),-circle(circleindex(i,k),6:8),circle(circleindex(i,k),9)];
    end
end
end