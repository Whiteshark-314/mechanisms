function Three_RRR_trace(active_links,passive_links,end_effector_sideLengths,fixed_coordinates)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);
Cx=sum(Fx)/3;
Cy=sum(Fy)/3;
alpha=0;
v = VideoWriter('3RRR.avi');
v.FrameRate=30;
open(v);
for i=0:1:359
    clf
    cx=(50*cosd(i))+Cx;
    cy=(50*sind(i))+Cy;
    Three_RRR_draw(a,p,e,F,[cx,cy],alpha)
    pause(0.1)
    frame = getframe(gcf);
    for j=1:10
        writeVideo(v,frame);
    end
end
close(v)