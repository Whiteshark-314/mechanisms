%% Initialize required values
a=[100 100 100];
p=[75 75 75];
e=[75 75 75];
F=[0,0;200,0;100,100*sqrt(3)];
t=pi*[65 130 115]/180;
inc=5;
v=0:inc:360;
n=length(v);
points_ex1=cell(n^3,1);
points_ex2=cell(n^3,1);
points_ex3=cell(n^3,1);
points_ey1=cell(n^3,1);
points_ey2=cell(n^3,1);
points_ey3=cell(n^3,1);
points_Ox=cell(n^3,1);
points_Oy=cell(n^3,1);
condvects=getcondvects_n_k(n,3,v);
D=getcondvects_n_k(n,3,v,'cell');

%% desired workspace envelope
syms x
lhs=cos(x-pi/6+atan((a(1)+p(1))*cos(x)/(F(2,1)-(a(1)+p(1))*sin(x))));
rhs=(e(1)^2+F(2,1)^2-2*(a(1)+p(1))*F(2,1)*cos(x))/(2*e(1)*sqrt((a(1)+p(1))^2+F(2,1)^2-2*(a(1)+p(1))*F(2,1)*cos(x)));
change_over=180*double(vpasolve(lhs-rhs==0,x,0.6))/pi;
L=((a(1)+p(1))^2+(e(1))^2-(a(1)+p(1))^2+F(2,1)^2+(a(1)+p(1))*e(1)/sqrt(3))/(2*F(2,1));
%extreme=2*atand((e(1)+2*sqrt((a(1)+p(1))^2+sqrt(3)*(a(1)+p(1))*e(1)+...
%    e(1)^2-L^2))/(2*(a(1)+p(1))+2*L+sqrt(3)*e(1)));
T=[acosd((F(2,1)-e(1))/(2*(a(1)+p(1)))):-0.1:change_over,change_over:-0.5:30,30]';
g=atand((a(1)+p(1))*sind(T)./(F(2,1)-(a(1)+p(1))*cosd(T)));
f=sqrt((a(1)+p(1))^2+F(2,1)^2-2*(a(1)+p(1))*F(2,1)*cosd(T));
beta=acosd((e(1)^2+f.^2-(a(1)+p(1))^2)./(2*e(1)*f))-g;
Xe1=(a(1)+p(1))*cosd(T)+e(1)*cosd(30+beta)/sqrt(3);
Ye1=(a(1)+p(1))*sind(T)+e(1)*sind(30+beta)/sqrt(3);
Xe2=(a(1)+p(1)+e(1)/sqrt(3))*cosd(T);
Ye2=(a(1)+p(1)+e(1)/sqrt(3))*sind(T);
lXe=length(Xe1);
dist_list=zeros(lXe,1);
I_list=zeros(lXe,1);
for jj=1:lXe
[dist,I]=min(sqrt((Xe1(jj)-Xe2).^2+(Ye1(jj)-Ye2).^2));
dist_list(jj)=dist;I_list(jj)=I;
end
[~,K]=min(dist_list);
X=[Xe2(end:-1:I_list(K)+1);Xe1(K:-1:1)];
Y=[Ye2(end:-1:I_list(K)+1);Ye1(K:-1:1)];
min1=min([X';Y'],[],2);
points_env=[X';Y']-min1;
points_env=[points_env(1,:),-points_env(1,end:-1:1);points_env(2,:),points_env(2,end:-1:1)];
points_env1=points_env+min1;
rotz_120=[cosd(120),-sind(120);sind(120),cosd(120)];
points_env2=rotz_120*points_env;
points_env3=rotz_120'*points_env;
points_envz=[points_env1,points_env2+points_env1(:,end)-points_env2(:,1),points_env3+points_env1(:,1)-points_env3(:,end)];
fill(points_envz(1,:),points_envz(2,:),'k')
hold on
axis equal

%% desired workspace envelope with alpha=0
omega=[30:0.5:acosd((F(2,1)-e(1))/(2*(a(1)+p(1)))),acosd((F(2,1)-e(1))/(2*(a(1)+p(1))))];
points=[(a(1)+p(1))*cosd(omega)+e(1)/2;(a(1)+p(1))*sind(omega)+e(1)/(2*sqrt(3))];
min2=min(points,[],2);
points=[points(1,:);points(2,:)]-min2;
points=[points(1,:),-points(1,end:-1:1);points(2,:),points(2,end:-1:1)];
points1=points+min2;
rotz_120=[cosd(120),-sind(120);sind(120),cosd(120)];
points2=rotz_120*points;
points3=rotz_120'*points;
pointsz=[points1,points2+points1(:,end)-points2(:,1),points3+points1(:,1)-points3(:,end)];
fill(pointsz(1,:),pointsz(2,:),'b')
axis equal
%% desired workspace for different tolerence values
