A=[100 100 100];
P=[100 100 100];
E=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
all_possible=getcondvects_n_k(3,9,[0,0.1,-0.1]);
a=A+all_possible(:,7:9);
p=P+all_possible(:,4:6);
e=E+all_possible(:,1:3);
m=size(a,1);
inc=5;
v=0:inc:360;
n=length(v);
D=getcondvects_n_k(3,9,[0,0.1,-0.1],'cell');
ii=1;
pos_errors=zeros(m,n);

parfor i=1:n
    theta=deg2rad(ithRowCondvect_n_k(i,n,3,v));
    mech_ref=Three_RRR(A,P,E,F);
    mech_ref=fk(mech_ref,theta);
    if ~(sum(isnan(mech_ref.alpha))==8)
        mech=Three_RRR(a,p,e,F);
        mech=fk(mech,theta);
        pos_errors(:,i)=min(sqrt((mech.O(:,:,1)-reshape(mech_ref.O(1,:),1,1,8)).^2+...
            (mech.O(:,:,2)-reshape(mech_ref.O(2,:),1,1,8)).^2),[],[2, 3]);
    else
        pos_errors(:,i)=NaN;
    end
end
[Max,I]=max(pos_errors);
%save('fullpart','Max','I')
% ex1=cell2mat(points_ex1);
% ex2=cell2mat(points_ex2);
% ex3=cell2mat(points_ex3);
% ey1=cell2mat(points_ey1);
% ey2=cell2mat(points_ey2);
% ey3=cell2mat(points_ey3);
% Ox=cell2mat(points_Ox);
% Oy=cell2mat(points_Oy);
% 
% xcell=cell(length(Ox),6);
% count=0;
% for j=1:length(points_Ox)
%     sols=length(points_Ox{j});
%     for k=1:sols
%         xcell(count+k,:)={D(j),th{j},[points_ex1{j}(k),...
%             points_ey1{j}(k),points_ex2{j}(k),points_ey2{j}(k),...
%             points_ex3{j}(k),points_ey3{j}(k),points_Ox{j}(k),...
%             points_Oy{j}(k),Max_pos_error{j}(k)],ph{j},alpha{j},index(j)};
%     end
%     count=count+sols;
% end
% writecell(xcell,'TOL2_0.1_max.xls')
