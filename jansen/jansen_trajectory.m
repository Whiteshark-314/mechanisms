function [Fx,Fy]=jansen_trajectory(angle_intervals,link_lengths)
    arguments
        angle_intervals (1,1) double
        link_lengths double =[38,41.5,39.3,40.1,55.8,39.4,36.7,65.7,49,50,61.9,7.8,15];
    end
E=zeros(2,ceil(360/angle_intervals)+1);
count=1;
h=figure(1);clf
% filename = 'jansen.gif';
for i=0:angle_intervals:360
    F=jansen_FK(i*pi/180,link_lengths);
    if i==0
        clf
        plot(F(1,1:6),F(2,1:6))
        hold on
        plot(F(1,[1,5]),F(2,[1,5]))
        plot(F(1,[1,6]),F(2,[1,6]))
        plot(F(1,[1,7]),F(2,[1,7]))
        plot(F(1,[4,7]),F(2,[4,7]))
        plot(F(1,[6,8]),F(2,[6,8]))
        plot(F(1,[7:9,7]),F(2,[7:9,7]))
        xlim([-100,100])
        ylim([-100 60])
        axis equal manual
    end
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if i == 0 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
    E(:,count)=F(:,9);
    count=count+1;
end
plot(E(1,:),E(2,:))
Fx=E(1,:);
Fy=E(2,:);
% imwrite(imind,cm,filename,'gif','WriteMode','append');
hold off
