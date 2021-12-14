a=[100 100 100];
p=[100 100 100];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
all_possible=getcondvects_n_k(3,9,[0,0.1,-0.1]);
a=a+all_possible(:,7:9);
p=p+all_possible(:,4:6);
e=e+all_possible(:,1:3);
m=size(all_possible,1);
inc=15;
v=0:inc:360;
n=length(v);
points_ex1=cell(m,1);
points_ex2=cell(m,1);
points_ex3=cell(m,1);
points_ey1=cell(m,1);
points_ey2=cell(m,1);
points_ey3=cell(m,1);
points_Ox=cell(m,1);
points_Oy=cell(m,1);
th=cell(m,1);
ph=cell(m,1);
alpha=cell(m,1);
Max_pos_error=cell(m,1);

index=zeros(n^3,1);
thetas=deg2rad(getcondvects_n_k(n,3,v));
M=size(thetas,1);
D=getcondvects_n_k(3,9,[0,0.1,-0.1],'cell');
ii=1;    
mech_ref=Three_RRR(a(1,:),p(1,:),e(1,:),F);
mech_ref=fk(mech_ref,thetas);
for i=2:size(thetas,1)
    mech=Three_RRR(a(i,:),p(i,:),e(i,:),F);
%    if true%~isempty(mech_ref.O)
        mech=fk(mech,thetas);
        pos_errors=sqrt((mech.O(:,ii,1)-mech_ref.O(:,ii,1)).^2+(mech.O(:,ii,2)-mech_ref.O(:,ii,2)).^2);
        [Max_pos_error{i},index(i)]=max(pos_errors,[],1,'linear');
        
        points_Ox{i}= mech.O(index(i),ii,1);
        points_Oy{i}= mech.O(index(i),ii,2);
        points_ex1{i}=mech.Ex(index(i),ii,1);
        points_ex2{i}=mech.Ex(index(i),ii,2);
        points_ex3{i}=mech.Ex(index(i),ii,3);
        points_ey1{i}=mech.Ey(index(i),ii,1);
        points_ey2{i}=mech.Ey(index(i),ii,2);
        points_ey3{i}=mech.Ey(index(i),ii,3);
        th{i}=rad2deg(mech.thetas(index(i),:));
        alpha{i}=rad2deg(mech.alpha(index(i),ii));
        ph{i}=squeeze(rad2deg(mech.phis(index(i),ii,:)));
%     else
%         points_Ox{i}= NaN;
%         points_Oy{i}= NaN;
%         points_ex1{i}=NaN;
%         points_ex2{i}=NaN;
%         points_ex3{i}=NaN;
%         points_ey1{i}=NaN;
%         points_ey2{i}=NaN;
%         points_ey3{i}=NaN;
%         ph{i}=NaN;
%    end
end
ex1=cell2mat(points_ex1);
ex2=cell2mat(points_ex2);
ex3=cell2mat(points_ex3);
ey1=cell2mat(points_ey1);
ey2=cell2mat(points_ey2);
ey3=cell2mat(points_ey3);
Ox=cell2mat(points_Ox);
Oy=cell2mat(points_Oy);

xcell=cell(length(Ox),6);
count=0;
for j=1:length(points_Ox)
    sols=length(points_Ox{j});
    for k=1:sols
        xcell(count+k,:)={D{j},th{j},[points_ex1{j}(k),...
            points_ey1{j}(k),points_ex2{j}(k),points_ey2{j}(k),...
            points_ex3{j}(k),points_ey3{j}(k),points_Ox{j}(k),...
            points_Oy{j}(k),Max_pos_error{j}(k)],ph{j},alpha{j},index(j)};
    end
    count=count+sols;
end
writecell(xcell,'TOL_0.1_max.xls')

% figure(1)
% plot(ex1,ey1,'or',ex2,ey2,'og',ex3,ey3,'ob')
% axis equal
% hold on %figure(2)
% plot(Ox,Oy,'.k')
% axis equal