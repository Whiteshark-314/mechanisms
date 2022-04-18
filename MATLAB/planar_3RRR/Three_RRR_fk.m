function [Phi,alpha,Ex,Ey,O]=Three_RRR_fk(active_links,passive_links,end_effector_sideLengths,fixed_coordinates,thetas)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);

Xx_1=(Fx(2)-Fx(1)+a(2)*cos(thetas(2))-a(1)*cos(thetas(1)))/p(1);
Yy_1=(Fy(2)-Fy(1)+a(2)*sin(thetas(2))-a(1)*sin(thetas(1)))/p(1);
pp_1=p(2)/p(1);
ep_1=-e(1)/p(1);

Xx_3=(Fx(2)-Fx(3)+a(2)*cos(thetas(2))-a(3)*cos(thetas(3)))/p(3);
Yy_3=(Fy(2)-Fy(3)+a(2)*sin(thetas(2))-a(3)*sin(thetas(3)))/p(3);
pp_3=p(2)/p(3);
ep_3=e(2)/p(3);

a11=2*Xx_1*pp_1+2*pp_1*ep_1;
a12=0;
a13=2*Xx_1*pp_1-2*pp_1*ep_1;
b11=2*Yy_1*pp_1;
b12=4*pp_1*ep_1;
b13=2*Yy_1*pp_1;
c11=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2-1+2*Xx_1*ep_1;
c12=4*Yy_1*ep_1;
c13=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2-1-2*Xx_1*ep_1;

a31=2*Xx_3*pp_3-pp_3*ep_3;
a32=-2*sqrt(3)*pp_3*ep_3;
a33=2*Xx_3*pp_3+pp_3*ep_3;
b31=2*Yy_3*pp_3+sqrt(3)*pp_3*ep_3;
b32=-2*pp_3*ep_3;
b33=2*Yy_3*pp_3-sqrt(3)*pp_3*ep_3;
c31=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2-1-Xx_3*ep_3+sqrt(3)*Yy_3*ep_3;
c32=-2*sqrt(3)*Xx_3*ep_3-2*Yy_3*ep_3;
c33=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2-1+Xx_3*ep_3-sqrt(3)*Yy_3*ep_3;

syms t 
A_1=a11+a12*t+a13*t^2;
B_1=b11+b12*t+b13*t^2;
C_1=c11+c12*t+c13*t^2;
A_3=a31+a32*t+a33*t^2;
B_3=b31+b32*t+b33*t^2;
C_3=c31+c32*t+c33*t^2;
eqn=(A_3*C_1-A_1*C_3)^2+(B_3*C_1-B_1*C_3)^2-(A_3*B_1-A_1*B_3)^2;
T=double(vpa(solve(eqn,t)));
T=T(imag(T)==0);
alpha=2*atan(T)';
Phi=zeros(3,length(alpha));

A1=2*Xx_1*pp_1+2*pp_1*ep_1*cos(alpha);
B1=2*Yy_1*pp_1+2*pp_1*ep_1*sin(alpha);
C1=Xx_1^2+Yy_1^2+pp_1^2+ep_1^2+2*ep_1*(Xx_1*cos(alpha)+Yy_1*sin(alpha))-1;
A3=2*Xx_3*pp_3+2*pp_3*ep_3*cos(2*pi/3+alpha);
B3=2*Yy_3*pp_3+2*pp_3*ep_3*sin(2*pi/3+alpha);
C3=Xx_3^2+Yy_3^2+pp_3^2+ep_3^2+2*ep_3*(Xx_3*cos(2*pi/3+alpha)+Yy_3*sin(2*pi/3+alpha))-1;

CP_2=(B3.*C1-B1.*C3)./(A3.*B1-A1.*B3);
SP_2=-(A3.*C1-A1.*C3)./(A3.*B1-A1.*B3);
Phi(2,:)=atan2(SP_2,CP_2);

CP_1=Xx_1+pp_1*cos(Phi(2,:))+ep_1*cos(alpha);
SP_1=Yy_1+pp_1*sin(Phi(2,:))+ep_1*sin(alpha);
Phi(1,:)=atan2(SP_1,CP_1);

CP_3=Xx_3+pp_3*cos(Phi(2,:))+ep_3*cos(2*pi/3+alpha);
SP_3=Yy_3+pp_3*sin(Phi(2,:))+ep_3*sin(2*pi/3+alpha);
Phi(3,:)=atan2(SP_3,CP_3);

Ex(1,:)=Fx(1)+a(1)*cos(thetas(1))+p(1)*cos(Phi(1,:));
Ex(2,:)=Fx(2)+a(2)*cos(thetas(2))+p(2)*cos(Phi(2,:));
Ex(3,:)=Fx(3)+a(3)*cos(thetas(3))+p(3)*cos(Phi(3,:));
Ex=Ex';
Ey(1,:)=Fy(1)+a(1)*sin(thetas(1))+p(1)*sin(Phi(1,:));
Ey(2,:)=Fy(2)+a(2)*sin(thetas(2))+p(2)*sin(Phi(2,:));
Ey(3,:)=Fy(3)+a(3)*sin(thetas(3))+p(3)*sin(Phi(3,:));
Ey=Ey';
O(:,1)=mean(Ex,2);
O(:,2)=mean(Ey,2);