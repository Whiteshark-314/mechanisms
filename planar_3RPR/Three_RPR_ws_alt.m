e=[1 1 1];
F=[0,0;2,0;1,1*sqrt(3)];
mech=Three_RPR(e,F);
min_rho=sqrt(1/3);
delta=0.1;
inc=0.1;
v=1:inc:1.6;
n=length(v);
points_ex1=cell(n^3,1);
points_ex2=cell(n^3,1);
points_ex3=cell(n^3,1);
points_ey1=cell(n^3,1);
points_ey2=cell(n^3,1);
points_ey3=cell(n^3,1);
points_Ox=cell(n^3,1);
points_Oy=cell(n^3,1);
Max_pos_error=cell(n^3,1);
th=cell(n^3,1);
condvects=getcondvects_n_k(n,3,v);
M=size(condvects,1);
lambda=getcondvects_n_k(5,6,72:72:360);
N=size(lambda,1);
D=getcondvects_n_k(n,3,v,'cell');
ii=2;
%combi=cell(length(lambda),length(condvects));
parfor i=1:size(condvects,1)
    mech=Three_RPR(e,F);
    mech_ref=Three_RPR(e,F);
    mech_base=fk(mech_ref,condvects(i,:));
    if ~isempty(mech_base.O)
        mech=fk(mech,condvects(i,:),delta,lambda);
        pos_errors=sqrt((mech.O(1,ii,:)-mech_base.O(1,ii)).^2+(mech.O(2,ii,:)-mech_base.O(2,ii)).^2);
        [max_pos_err,Index]=max(pos_errors,[],'all','linear');

        points_Ox{i}= sort(mech.O(1,ii,Index));
        points_Oy{i}= sort(mech.O(2,ii,Index));
        points_ex1{i}=sort(mech.Ex(1,ii,Index));
        points_ex2{i}=sort(mech.Ex(2,ii,Index));
        points_ex3{i}=sort(mech.Ex(3,ii,Index));
        points_ey1{i}=sort(mech.Ey(1,ii,Index));
        points_ey2{i}=sort(mech.Ey(2,ii,Index));
        points_ey3{i}=sort(mech.Ey(3,ii,Index));
        th{i}=rad2deg(mech.thetas(:,ii,Index)');
        Max_pos_error{i}=max_pos_err;
    else
        points_Ox{i}= NaN;
        points_Oy{i}= NaN;
        points_ex1{i}=NaN;
        points_ex2{i}=NaN;
        points_ex3{i}=NaN;
        points_ey1{i}=NaN;
        points_ey2{i}=NaN;
        points_ey3{i}=NaN;
        Max_pos_error{i}=NaN;
        th{i}=NaN;
    end
end
ex1=cell2mat(points_ex1');
ex2=cell2mat(points_ex2');
ex3=cell2mat(points_ex3');
ey1=cell2mat(points_ey1');
ey2=cell2mat(points_ey2');
ey3=cell2mat(points_ey3');
Ox=cell2mat(points_Ox');
Oy=cell2mat(points_Oy');
MPE=cell2mat(Max_pos_error');
thetas=cell2mat(th);

xcell=cell(length(Ox),2);
count=0;
for j=1:length(points_Ox)
    sols=length(points_Ox{j});
    for k=1:sols
        xcell(count+k,:)={D{j},[points_ex1{j}(k),...
            points_ey1{j}(k),points_ex2{j}(k),points_ey2{j}(k),...
            points_ex3{j}(k),points_ey3{j}(k),points_Ox{j}(k),...
            points_Oy{j}(k),Max_pos_error{j}(k),th{j}]};
    end
    count=count+sols;
end
writecell(xcell,'C0.1_clear_max.xls')

% figure(1)
% plot(ex1,ey1,'or',ex2,ey2,'og',ex3,ey3,'ob')
% axis equal
% figure(2)
% plot(Ox,Oy,'*k')
% axis equal
