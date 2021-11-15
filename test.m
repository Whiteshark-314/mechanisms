centers=[0,0;0,0;200,0;200,0;100,100*sqrt(3);100,100*sqrt(3)];
radii=[175+100/sqrt(3);100/sqrt(3)-25;175+100/sqrt(3);100/sqrt(3)-25;175+100/sqrt(3);100/sqrt(3)-25];
viscircles(centers,radii,'color','k');
hold on
scatter(Ox,Oy,'.b')
% new_c=[100,100/sqrt(3)];
% new_r=roots([1,-(200/sqrt(3)+350),20625+10000/3+35000/sqrt(3)-new_c(2)^2]);
% new_r=175*sin(acos(2/7))-50/sqrt(3);
% viscircles(new_c,new_r,'color','k')
% t=0:0.1:180;
% w=-50;
% a=0.5*(175*sin(acos(2/7))+50/sqrt(3)+100+100/sqrt(3)-w);
% c=0.5*(175*sin(acos(2/7))+50/sqrt(3)-100-100/sqrt(3)-w);
% h=0.5*(175*sin(acos(2/7))+50/sqrt(3)-100-100/sqrt(3)+w);
% k=100;
% b=sqrt(a^2-c^2);
% x=a*cos(t)+h;
% y=b*sin(t)-k;
% theta=pi/2;
% xy=[x',y'];
% rmat=[cos(theta) sin(theta);-sin(theta) cos(theta)];
% fg=xy*rmat;
% scatter(fg(:,1),fg(:,2),'.k')
% axis equal
%Parabola
apex=175*sin(acos(2/7))+50/sqrt(3);
a_=100+100/sqrt(3);
y=[122:1:160,161:0.1:184.9,185:0.01:192,192.001:0.001:apex];
x=[100+sqrt(-1*a_*(y-apex)); 100-sqrt(-1*a_*(y-apex))];
XY=[[x(1,:)';x(2,:)'],[y';y']];
theta=2*pi/3;
rmat=[cos(theta) sin(theta);-sin(theta) cos(theta)];
XY_=[XY(:,1)-100,XY(:,2)-100/sqrt(3)];
XY_1=XY_*rmat;
XY_2=XY_1*rmat;
XY_1(:,1)=XY_1(:,1)+100;
XY_1(:,2)=XY_1(:,2)+100/sqrt(3);
XY_2(:,1)=XY_2(:,1)+100;
XY_2(:,2)=XY_2(:,2)+100/sqrt(3);
scatter(XY(:,1),XY(:,2),'.k')
scatter(XY_1(:,1),XY_1(:,2),'.k')
scatter(XY_2(:,1),XY_2(:,2),'.k')
axis equal
% Catenary
% a=175*sin(acos(2/7))+50/sqrt(3);
% x=-100:1:100;
% y=a+16*(1-cosh(x/40));
% scatter(x+100,y,'.k')
% axis equal