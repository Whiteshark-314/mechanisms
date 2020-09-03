function [Phi,aplha,Mnx,Mny,O]=Three_RRR(active_link,passive_link,end_effector_side,base_side,theta)
a=active_link;
p=passive_link;
h=end_effector_side;
b=base_side;
Phi=zeros(3,8);


c0=
c1=
c2=
c3=
c4=
c5=
c6=

alpha=roots([c0,c1,c2,c3,c4,c5,c6,]);
alpha=alpha';

A1=-2*b*p^3+2*h*p^3*cos(alpha)+2*a*p^3*cos(theta(1))-2*a*p^3*cos(theta(2));
B1=2*h*p^3*sin(alpha)+2*a*p^3*sin(theta(1))-2*a*p^3*sin(theta(2));
C1=2*a^2*p^2*cos(theta(1)-theta(2))+2*h*a*p^2*cos(alpha-theta(2))...
    -2*h*a*p^2*cos(alpha-theta(1))+2*h*b*p^2*cos(alpha)...
    +2*b*a*p^2*cos(theta(1))-2*b*a*p^2*cos(theta(2))-(h^2+b^2+2*a^2)*p^2;
A3=-b*p^3+h*p^3*cos(alpha)-2*a*p^3*cos(theta(2))+2*a*p^3*cos(theta(3))...
    +sqrt(3)*h*p^3*sin(alpha);
B3=sqrt(3)*b*p^3-sqrt(3)*h*p^3*cos(alpha)-2*a*p^3*sin(theta(2))...
    +2*a*p^3*sin(theta(3))+h*p^3*sin(alpha);
C3=2*a^2*p^2*cos(theta(1)-theta(2))+h*a*p^2*cos(alpha-theta(2))...
    -h*a*p^2*cos(alpha-theta(3))+2*h*b*p^2*cos(alpha)...
    +b*a*p^2*cos(theta(3))-b*a*p^2*cos(theta(2))-(h^2+b^2+2*a^2)*p^2 ...
    +sqrt(3)*(h*a*p^2*sin(alpha-theta(2))+h*a*p^2*sin(alpha+theta(3))...
    +b*a*p^2*sin(theta(2))-b*a^2*sin(theta(3)));
x=(C1.*B3-C3.*B1)./(A1.*B3-A3.*B1);
y=(C1.*A3-C3.*A1)./(A3.*B1-A1.*B3);
Phi(2,:)=atan2(y,x);

S_phi1=(-a*sin(theta(1))-h*sin(alpha)+p*sin(Phi(2))+a*sin(theta(2)))/p;
C_phi1=(-a*cos(theta(1))-h*cos(alpha)+p*cos(Phi(2))+a*cos(theta(2))+b)/p;
Phi(1,:)=atan2(S_phi1,C_phi1);

S_phi3=(-a*sin(theta(3))+h*sin(alpha+2*pi/3)+p*sin(Phi(2))+a*sin(theta(2))+sqrt(3)*b/2)/p;
C_phi3=(-a*cos(theta(3))+h*cos(alpha+2*pi/3)+p*cos(Phi(2))+a*cos(theta(2))+b-b/2)/p;
Phi(3,:)=atan2(S_phi3,C_phi3);

Mnx=[a*cos(theta(1))+p*cos(Phi(1,:));...
    a*cos(theta(2))+p*cos(Phi(2,:))+b;...
    a*cos(theta(3))+p*cos(Phi(3,:))+b/2];
Mny=[a*sin(theta(1))+p*sin(Phi(1,:));...
    a*sin(theta(2))+p*sin(Phi(2,:));...
    a*sin(theta(3))+p*sin(Phi(3,:))+sqrt(3)*b/2];
O=[sum(Mnx);sum(Mny)]/3;
%x=C1B3-C3B1/A1B3-A3B1
%y=C1A3-C3A1/A3B1-A1B3