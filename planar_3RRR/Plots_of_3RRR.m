%% Initialize required values
a=[100 90 100];
p=[100 100 95];
e=[100 100 100];
f=[200 200 200];
A=acosd((e(1)^2+e(3)^2-e(2)^2)/(2*e(1)*e(3)));
B=acosd((e(1)^2+e(2)^2-e(3)^2)/(2*e(1)*e(2)));
C=acosd((e(2)^2+e(3)^2-e(1)^2)/(2*e(2)*e(3)));

iota=acosd((f(1)^2+f(3)^2-f(2)^2)/(2*f(1)*f(3)));
A_1=[1,1,1;-1,1,0;-1,0,1];
A_2=[1,1,1;-1,1,0;-1,0,1];
b_1=[0;f(1);f(3)*cosd(iota)];
b_2=[0;0;f(3)*sind(iota)];
ex=A_1\b_1;
ey=A_2\b_2;

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
plot(points_envz(1,:),points_envz(2,:),'k')
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
%% desired ws envelope with tol and alpha=0
zeta=180-B;
e_points=[0,e(1),e(1)+e(2)*cosd(zeta);0,0,e(2)*sind(zeta)];
e_c=mean(e_points,2);

oMega_A1=[acosd((F(2,1)/2-e_c(1))/(a(1)+p(1))):-0.1:60,60:-0.5:A/2];
oMega_A1=oMega_A1(end:-1:1);
oMega_B1=[180-acosd((F(2,1)/2-(e(1)-e_c(1)))/(a(2)+p(2))):-0.1:60,60:-0.5:A/2];
oMega_B1=oMega_B1(1:length(oMega_A1));
oMega_B1=oMega_B1(end:-1:1);
XA1=(a(1)+p(1))*cosd(oMega_A1)+e_c(1);
XB1=F(2,1)+(a(2)+p(2))*cosd(oMega_B1)-(e(1)-e_c(1));
YA1=(a(1)+p(1))*sind(oMega_A1)+e_c(2);
YB1=(a(2)+p(2))*sind(oMega_B1)+e_c(2);
XAB1=[XA1;XB1];
YAB1=[YA1;YB1];
[YAB_1,I1]=min(YAB1,[],1,'linear');
XAB_1=XAB1(I1);

oMega_B2=[max(oMega_B1):0.1:120,120:0.5:180-B/2];
oMega_A2=[max(oMega_A1):0.1:120,120:0.5:180-B/2];
oMega_A2=oMega_A2(1:length(oMega_B2));
XA2=(a(1)+p(1))*cosd(oMega_A2)+e_c(1);
XB2=F(2,1)+(a(2)+p(2))*cosd(oMega_B2)-(e(1)-e_c(1));
YA2=(a(1)+p(1))*sind(oMega_A2)+e_c(2);
YB2=(a(2)+p(2))*sind(oMega_B2)+e_c(2);
XAB2=[XA2;XB2];
YAB2=[YA2;YB2];
[YAB_2,I2]=min(YAB2,[],1,'linear');
XAB_2=XAB2(I2);

R120=[cosd(120),-sind(120);sind(120),cosd(120)];
R_120=[cosd(-120),-sind(-120);sind(-120),cosd(-120)];

F120=(R120*F')'+[200,0];
e120=R120*e_points+[100;0];
ec120=R120*e_c+[100;0];
oMega_C3=[acosd((F120(1,1)/2-ec120(1))/(a(3)+p(3))):-0.1:60,60:-0.5:C/2];
oMega_C3=oMega_C3(end:-1:1);
oMega_A3=[180-acosd((F120(1,1)/2-(e120(1,1)-ec120(1)))/(a(1)+p(1))):-0.1:60,60:-0.5:C/2];
oMega_A3=oMega_A3(1:length(oMega_C3));
oMega_A3=oMega_A3(end:-1:1);
XA3=F120(1,1)+(a(1)+p(1))*cosd(oMega_A3)-(e120(1,1)-ec120(1));
XC3=(a(3)+p(3))*cosd(oMega_C3)+ec120(1);
YA3=(a(1)+p(1))*sind(oMega_A3)+ec120(2);
YC3=(a(3)+p(3))*sind(oMega_C3)+ec120(2);
XAC3=[XA3;XC3]-200;
YAC3=[YA3;YC3];
[YAC_3,I3]=min(YAC3,[],1,'linear');
XAC_3=XAC3(I3);

oMega_C4=[max(oMega_C3):0.1:120,120:0.5:180-A/2];
oMega_A4=[max(oMega_A3):0.1:120,120:0.5:180-A/2];
oMega_C4=oMega_C4(1:length(oMega_A4));
XA4=F120(1,1)+(a(1)+p(1))*cosd(oMega_A4)-(e120(1,1)-ec120(1));
XC4=(a(3)+p(3))*cosd(oMega_C4)+ec120(1);
YA4=(a(1)+p(1))*sind(oMega_A4)+ec120(2);
YC4=(a(3)+p(3))*sind(oMega_C4)+ec120(2);
XAC4=[XA4;XC4]-200;
YAC4=[YA4;YC4];
[YAC_4,I4]=min(YAC4,[],1,'linear');
XAC_4=XAC4(I4);

F_120=(R_120*(F-[200,0])')';
e_120=R_120*(e_points-[100;0]);
e_c120=R_120*(e_c-[100;0]);
oMega_B5=[acosd((F_120(3,1)/2-e_c120(1))/(a(2)+p(2))):-0.1:60,60:-0.5:B/2];
oMega_B5=oMega_B5(end:-1:1);
oMega_C5=[180-acosd((F_120(3,1)/2-(e_120(1,3)-e_c120(1)))/(a(3)+p(3))):-0.1:60,60:-0.5:B/2];
oMega_C5=oMega_C5(1:length(oMega_B5));
oMega_C5=oMega_C5(end:-1:1);
XC5=F_120(3,1)+(a(3)+p(3))*cosd(oMega_C5)-(e_120(1,3)-e_c120(1));
XB5=(a(2)+p(2))*cosd(oMega_B5)+e_c120(1);
YC5=(a(3)+p(3))*sind(oMega_C5)+e_c120(2);
YB5=(a(2)+p(2))*sind(oMega_B5)+e_c120(2);
XBC5=[XB5;XC5];
YBC5=[YB5;YC5];
[YBC_5,I5]=min(YBC5,[],1,'linear');
XBC_5=XBC5(I5);

oMega_B6=[max(oMega_B5):0.1:120,120:0.5:180-C/2];
oMega_C6=[max(oMega_C5):0.1:120,120:0.5:180-C/2];
oMega_B6=oMega_B6(1:length(oMega_C6));
XC6=F_120(3,1)+(a(3)+p(3))*cosd(oMega_C6)-(e_120(1,3)-e_c120(1));
XB6=(a(2)+p(2))*cosd(oMega_B6)+e_c120(1);
YC6=(a(3)+p(3))*sind(oMega_C6)+e_c120(2);
YB6=(a(2)+p(2))*sind(oMega_B6)+e_c120(2);
XBC6=[XB6;XC6];
YBC6=[YB6;YC6];
[YBC_6,I6]=min(YBC6,[],1,'linear');
XBC_6=XBC6(I6);

XBC=[XBC_5,XBC_6];
YBC=[YBC_5,YBC_6];
BC=R120*[XBC;YBC]+[200;0];
XAC=[XAC_3,XAC_4];
YAC=[YAC_3,YAC_4];
AC=R_120*[XAC;YAC];
XAB=[XAB_1,XAB_2];
YAB=[YAB_1,YAB_2];
AB=[XAB;YAB];
hold on
% plot(AC(1,:),AC(2,:))
% plot(AB(1,:),AB(2,:))
% plot(BC(1,:),BC(2,:))
fill([AB(1,:),BC(1,:),AC(1,:)],[AB(2,:),BC(2,:),AC(2,:)],'c')
axis equal
%% desired workspace for different tolerence values
%% tolerence vs max position error
%% tolerence vs actuator position vs maximum position error
%% tolerence vs orientation angle 1,2,3
%% tolerence vs actuator position vs orientation angle 1,2,3
%% compensation IK approach
