function Three_RRR_draw(active_links,passive_links,end_effector_sideLengths,fixed_coordinates,O,alpha)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);
[thetas,Ex,Ey] = Three_RRR_ik(a,p,e,F,O,alpha);
unique_list=unique(nchoosek([1,2,1,2,1,2],3),'rows');
thetas_combi=[thetas(1,unique_list(:,1))', thetas(2,unique_list(:,2))', thetas(3,unique_list(:,3))'];
figure(1)
for i=1:8
    if i>=5
        j=i+1;
    else
        j=i;
    end
    subplot(3,3,j)
    plot(Fx,Fy,'.r','MarkerSize',20)
    hold on
    plot(Ex,Ey,'.b','MarkerSize',15)
    Px=Fx+(a.*cos(thetas_combi(i,:)))';
    Py=Fy+(a.*sin(thetas_combi(i,:)))';
    plot(Px,Py,'.g','MarkerSize',10)
    plot([Fx,Px]',[Fy,Py]','r', 'LineWidth', 2)
    plot([Px,Ex']',[Py,Ey']','g', 'LineWidth', 2)
    fill(Ex,Ey,'b')
    plot(O(1), O(2), '.k', 'MarkerSize',25)
    xlim([-100,300])
    ylim([-100,250])
    axis equal
    hold off
end
