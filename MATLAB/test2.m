a=[100 100 100];
p=[75 75 75];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];

X1=((F(2,1)-F(1,1))^2-e(1)^2+(a(1)+p(1))^2-(a(2)+p(2))^2)/(2*(F(2,1)-F(1,1))*(a(1)+p(1)));
Y1=(e(1)*(a(2)+p(2)))/((F(2,1)-F(1,1))*(a(1)+p(1)));

X2=((F(2,1)-F(1,1))^2-e(1)^2/3+(a(1)+p(1)+e(1)/sqrt(3))^2-(a(2)+p(2))^2)/(2*(F(2,1)-F(1,1))*(a(1)+p(1)+e(1)/sqrt(3)));
Y2=(e(1)*(a(2)+p(2))/sqrt(3))/((F(2,1)-F(1,1))*(a(1)+p(1)+e(1)/sqrt(3)));

X3=((F(2,1)-F(1,1))^2-e(1)^2+(a(1)+p(1)+e(1))^2-(a(2)+p(2))^2)/(2*(F(2,1)-F(1,1))*(a(1)+p(1)+e(1)));
Y3=(e(1)*(a(2)+p(2)))/((F(2,1)-F(1,1))*(a(1)+p(1)+e(1)));

A_1=Y1-0.5*sqrt(3)*Y2;
B_1=0.5*Y2;
C_1=X2-X1;

A_2=Y1-0.5*Y3;
B_2=0.5*sqrt(3)*Y3;
C_2=X3-X1;

phi=[atan2(B_1,A_1)+atan2(sqrt(A_1^2+B_1^2-C_1^2),C_1), atan2(B_2,A_2)+atan2(sqrt(A_2^2+B_2^2-C_2^2),C_2)]

theta=acos(X1+Y1*cos(phi))

xi=pi-(2*pi-5*pi/6-theta-phi)

xxx=[175*cos(theta(1))+100/sqrt(3)*cos(theta(1)), 175*cos(theta(2))+100/sqrt(3)*cos(theta(2)-pi/6)]
yyy=[175*sin(theta(1))+100/sqrt(3)*sin(theta(1)), 175*sin(theta(2))+100/sqrt(3)*sin(theta(2)-pi/6)]

x=[-(xxx-100),0,xxx(end:-1:1)-100]
y=[yyy,175*sin(acos(2/7))+50/sqrt(3),yyy(end:-1:1)]

plot(x,y,'.k')