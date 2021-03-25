function [thetas,Ex,Ey] = Three_RRR_ik(active_links,passive_links,end_effector_sideLengths,fixed_coordinates,O,alpha)
a=active_links;
p=passive_links;
e=end_effector_sideLengths;
F=fixed_coordinates;
Fx=F(:,1);
Fy=F(:,2);

iota=acos((e(1)^2-e(2)^2+e(3)^2)/(2*e(1)*e(3)));
A_1=[1,1,1;-1,1,0;-1,0,1];
A_2=[1,1,1;-1,1,0;-1,0,1];
b_1=[3*O(1);e(1)*cos(alpha);e(3)*cos(iota+alpha)];
b_2=[3*O(2);e(1)*sin(alpha);e(3)*sin(iota+alpha)];
ex=A_1\b_1;
ey=A_2\b_2;
Ex=ex';
Ey=ey';
eFx=(ex-Fx)';
eFy=(ey-Fy)';

A=eFx.^2+eFy.^2+a.^2+2*eFx.*a-p.^2;
B=-4*eFy.*a;
C=eFx.^2+eFy.^2+a.^2-2*eFx.*a-p.^2;

syms t1 t2 t3
At1=A(1)*t1^2;
Bt1=B(1)*t1;
Ct1=C(1);
At2=A(2)*t2^2;
Bt2=B(2)*t2;
Ct2=C(2);
At3=A(3)*t3^2;
Bt3=B(3)*t3;
Ct3=C(3);
eqn_t1=At1+Bt1+Ct1;
eqn_t2=At2+Bt2+Ct2;
eqn_t3=At3+Bt3+Ct3;
T1=double(vpa(solve(eqn_t1,t1)));
T2=double(vpa(solve(eqn_t2,t2)));
T3=double(vpa(solve(eqn_t3,t3)));
T1=T1(imag(T1)==0);
T2=T2(imag(T2)==0);
T3=T3(imag(T3)==0);
if isempty(T1)||isempty(T2)||isempty(T3)
    thetas=[];
else
    thetas=[2*atan(T1)';2*atan(T2)';2*atan(T3)'];
end
end

