function [thetas_1,thetas_2,thetas_3]=Three_RRR_choose(active_links,passive_links,end_effector_sideLengths,fixed_coordinates)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);
Cx=sum(Fx)/3;
Cy=sum(Fy)/3;
alpha=0;
thetas_1=zeros(360,8);
thetas_2=zeros(360,8);
thetas_3=zeros(360,8);
t_1=zeros(2,2);
t_2=zeros(2,2);
t_3=zeros(2,2);
v = VideoWriter('3RRR.avi');
v.FrameRate=30;
open(v);
fh = figure(1);
fh.WindowState = 'maximized';
for i=0:1:359
    cx=(50*cosd(i))+Cx;
    cy=(50*sind(i))+Cy;
    [thetas,Ex,Ey,thetas_combi]=Three_RRR_ik(a,p,e,F,[cx,cy],alpha);
    t_1(2,:)=thetas(1,:);
    t_2(2,:)=thetas(2,:);
    t_3(2,:)=thetas(3,:);
    thetas_1(i+1,:)=thetas_combi(:,1)';
    thetas_2(i+1,:)=thetas_combi(:,2)';
    thetas_3(i+1,:)=thetas_combi(:,3)';
    if i>=1
        t_1(1,:)=thetas_1(i,4:5);
        t_2(1,:)=thetas_2(i,2:3);
        t_3(1,:)=thetas_3(i,1:2);
        if (abs(diff([t_1(1,1),t_1(2,2)]))<0.2)  ||  (abs(diff([t_1(1,2),t_1(2,1)]))<0.2)
            thetas_1(i+1,:)=circshift(thetas_1(i+1,:),4);
        end
        if (abs(diff([t_2(1,1),t_2(2,2)]))<0.2)  ||  (abs(diff([t_2(1,2),t_2(2,1)]))<0.2)
            thetas_2(i+1,:)=circshift(thetas_2(i+1,:),2);
        end
        if (abs(diff([t_3(1,1),t_3(2,2)]))<0.2)  ||  (abs(diff([t_3(1,2),t_3(2,1)]))<0.2)
            thetas_3(i+1,:)=circshift(thetas_3(i+1,:),1);
        end
    else
        t_1(1,:)=thetas(1,:);
        t_2(1,:)=thetas(2,:);
        t_3(1,:)=thetas(3,:);
    end
    thetas_=[thetas_1(i+1,:)',thetas_2(i+1,:)',thetas_3(i+1,:)'];
    clf
    for k=1:8
        if k>=5
            j=k+1;
        else
            j=k;
        end
        subplot(3,3,j)
        plot(Fx,Fy,'.r','MarkerSize',20)
        hold on
        plot(Ex,Ey,'.b','MarkerSize',15)
        Px=Fx+(a.*cos(thetas_(k,:)))';
        Py=Fy+(a.*sin(thetas_(k,:)))';
        plot(Px,Py,'.g','MarkerSize',10)
        plot([Fx,Px]',[Fy,Py]','r', 'LineWidth', 2)
        plot([Px,Ex']',[Py,Ey']','g', 'LineWidth', 2)
        fill(Ex,Ey,'b')
        plot(cx, cy, '.k', 'MarkerSize',25)
        xlim([-100,300])
        ylim([-100,250])
        axis equal
        hold off
    end
    pause(0.1)
    frame = getframe(gcf);
    for j=1:10
        writeVideo(v,frame);
    end
end
close(v)
