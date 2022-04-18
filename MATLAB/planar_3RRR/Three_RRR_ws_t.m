a=[100 100 100];
p=[75 75 75];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
Fx=F(:,1);
Fy=F(:,2);
iota=acos((e(1)^2-e(2)^2+e(3)^2)/(2*e(1)*e(3)));
ex=[0;e(1);e(3)*cos(iota)];
ey=[0;0;e(3)*sin(iota)];
O=mean([ex,ey]);
%r=[a(1)+p(1)+h(1)/sqrt(3) a(2)+p(2)+h(2)/sqrt(3) a(3)+p(3)+h(3)/sqrt(3)];
ap1=a(1)+p(1);
apex=sqrt(ap1^2-(Fx(2)-Fx(1)-e(1))^2/4)+e(1)*sqrt(3)/2;

[k,l]=circcirc(0,0,e(3),e(1),0,e(2));
e_xy=[0,0;e(1),0;k(1),abs(l(1))];
O_xy=mean(e_xy);
r=[a(1)+p(1)+sqrt(sum((e_xy(1,:)-O_xy).^2)) a(2)+p(2)+sqrt(sum((e_xy(2,:)-O_xy).^2)) a(3)+p(3)+sqrt(sum((e_xy(3,:)-O_xy).^2))];
r2=abs([a(1)-p(1) a(2)-p(2) a(3)-p(3)]);
[k_l,l_l]=circcirc(F(3,1),F(3,2),r(3),F(2,1),F(2,2),r(2));
[k_r,l_r]=circcirc(F(3,1),F(3,2),r(3),F(1,1),F(1,2),r(1));
[k_t,l_t]=circcirc(F(2,1),F(2,2),r(2),F(1,1),F(1,2),r(1));
[min_l,I1]=min(k_l);
[max_r,I2]=max(k_r);
[max_t,I3]=max(l_t);
theta1=linspace(atan2d(l_r(I2)-F(1,2),k_r(I2)-F(1,1)),atan2d(l_t(I3)-F(1,2),k_t(I3)-F(1,1)),1000);
theta2=linspace(atan2d(l_t(I3)-F(2,2),k_t(I3)-F(2,1)),360+atan2d(l_l(I1)-F(2,2),k_l(I1)-F(2,1)),1000);
theta3=linspace(atan2d(l_l(I1)-F(3,2),k_l(I1)-F(3,1)),atan2d(l_r(I2)-F(3,2),k_r(I2)-F(3,1)),1000);
x1=[r(1)*cosd(theta1)+F(1,1), r(2)*cosd(theta2)+F(2,1), r(3)*cosd(theta3)+F(3,1)];
y1=[r(1)*sind(theta1)+F(1,2), r(2)*sind(theta2)+F(2,2), r(3)*sind(theta3)+F(3,2)];

% a=175*sin(acos(2/7))+50/sqrt(3);
% a_=100+100/sqrt(3);
% y=[122:1:160,161:0.1:184.9,185:0.01:192,192.001:0.001:a];
% x=[100+sqrt(-1*a_*(y-a)); 100-sqrt(-1*a_*(y-a))];
% XY=[[x(1,:)';x(2,:)'],[y';y']];
% theta=2*pi/3;
% rmat=[cos(theta) sin(theta);-sin(theta) cos(theta)];
% XY_=[XY(:,1)-100,XY(:,2)-100/sqrt(3)];
% XY_1=XY_*rmat;
% XY_2=XY_1*rmat;
% XY_1(:,1)=XY_1(:,1)+100;
% XY_1(:,2)=XY_1(:,2)+100/sqrt(3);
% XY_2(:,1)=XY_2(:,1)+100;
% XY_2(:,2)=XY_2(:,2)+100/sqrt(3);
% XY=[XY;XY_1;XY_2];
% 
% in=inpolygon(x1,y1,XY(:,1),XY(:,1));
% 
fill(x1,y1,'k')
%fill(XY(:,1),XY(:,2),'w')
axis equal
