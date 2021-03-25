function [max_deltaO, I, deltaO]=Three_RRR_tol_max(tol_val,n)
a=[100 100 100];
p=[100 100 100];
h=[100 100 100];
b=[0,0;200,0;100,100*sqrt(3)];
t=[65 130 115];
thetas=t*pi/180;
all_possible=getcondvects3(9,tol_val);
all_possible_a=a+all_possible(:,7:9);
all_possible_p=p+all_possible(:,4:6);
all_possible_h=h+all_possible(:,1:3);
E=zeros(19683,2);
parfor i=1:19683
    A=all_possible_a(i,:);
    P=all_possible_p(i,:);
    H=all_possible_h(i,:);
    [~,~,~,~,O]=Three_RRR(A,P,H,b,thetas);
    E(i,:)=O(n,:);
end
x_sq=(E(:,1)-E(1,1)).^2;
y_sq=(E(:,2)-E(1,2)).^2;
deltaO=sqrt(x_sq+y_sq);
[max_deltaO, I]=max(deltaO);
deltaO=sort(deltaO,'descend');

