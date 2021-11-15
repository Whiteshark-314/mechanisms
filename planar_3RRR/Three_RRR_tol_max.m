function [max_deltaO, I, deltaO]=Three_RRR_tol_max(tol_val,actuation_angles,n)
a=[100 100 100];
p=[100 100 100];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
thetas=actuation_angles*pi/180;
all_possible=getcondvects_n_k(3,9,tol_val);
all_possible_a=a+all_possible(:,7:9);
all_possible_p=p+all_possible(:,4:6);
all_possible_e=e+all_possible(:,1:3);
E=zeros(19683,2);
test=Three_RRR(a,p,e,F);
test=fk(test,thetas);
if ~isempty(test.O)
    parfor i=1:19683
        a=all_possible_a(i,:);
        p=all_possible_p(i,:);
        e=all_possible_e(i,:);
        mech=Three_RRR(a,p,e,F);
        mech=fk(mech,thetas);
        E(i,:)=mech.O(n,:);
    end
    x_sq=(E(:,1)-E(1,1)).^2;
    y_sq=(E(:,2)-E(1,2)).^2;
    deltaO=sqrt(x_sq+y_sq);
    [max_deltaO, I]=max(deltaO);
    deltaO=sort(deltaO,'descend');
else
    max_deltaO=NaN;
    I=NaN;
    deltaO=NaN;
end
