function Three_RRR_draw(active_links,passive_links,end_effector_sideLengths,fixed_coordinates,O,alpha)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);
[~,Ex,Ey,thetas_combi] = Three_RRR_ik(a,p,e,F,O,alpha);
fh = figure(1);
fh.WindowState = 'maximized';
for i=1:8
    if i>=5
        j=i+1;
    else
        j=i;
    end
    subplot(3,3,j)
    clf
    for k=1:1
        plot(Fx,Fy,'.r','MarkerSize',20)
        hold on
        plot(Ex,Ey,'.b','MarkerSize',15)
        Px=Fx+(a.*cos(thetas_combi(k,:)))';
        Py=Fy+(a.*sin(thetas_combi(k,:)))';
        plot(Px,Py,'.g','MarkerSize',10)
        plot([Fx,Px]',[Fy,Py]','r', 'LineWidth', 2)
        plot([Px,Ex']',[Py,Ey']','g', 'LineWidth', 2)
        fill(Ex,Ey,'b')
        plot(cx, cy, '.k', 'MarkerSize',25)
        xlim([-350,550])
        ylim([-100,300])
        axis equal
        hold off
    end
    pause(0.1)
    if i==0
        gif('3RRR.gif')
    else
        gif
    end
end
