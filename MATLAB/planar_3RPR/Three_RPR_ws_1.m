e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
mech=Three_RPR_1(e,F);
min_rho=sqrt(1/3);
inc=10;
v=100:inc:160;
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
for i=1:size(condvects,1)
    thetas=condvects(i,:)*pi/180;
    mech=Three_RPR_1(e,F);
    mech=fk(mech,condvects(i,:));
    if ~isempty(mech.O)
        points_Ox{i}= mech.O(:,1);
        points_Oy{i}= mech.O(:,2);
        points_ex1{i}=mech.Ex(:,1);
        points_ex2{i}=mech.Ex(:,2);
        points_ex3{i}=mech.Ex(:,3);
        points_ey1{i}=mech.Ey(:,1);
        points_ey2{i}=mech.Ey(:,2);
        points_ey3{i}=mech.Ey(:,3);
    else
        points_Ox{i}= NaN;
        points_Oy{i}= NaN;
        points_ex1{i}=NaN;
        points_ex2{i}=NaN;
        points_ex3{i}=NaN;
        points_ey1{i}=NaN;
        points_ey2{i}=NaN;
        points_ey3{i}=NaN;
    end
end
ex1=cell2mat(points_ex1);
ex2=cell2mat(points_ex2);
ex3=cell2mat(points_ex3);
ey1=cell2mat(points_ey1);
ey2=cell2mat(points_ey2);
ey3=cell2mat(points_ey3);
Ox=cell2mat(points_Ox);
Oy=cell2mat(points_Oy);

xcell=cell(length(Ox),2);
count=0;
for j=1:length(points_Ox)
    sols=length(points_Ox{j});
    for k=1:sols
        xcell(count+k,:)={D{j},[points_ex1{j}(k),...
            points_ey1{j}(k),points_ex2{j}(k),points_ey2{j}(k),...
            points_ex3{j}(k),points_ey3{j}(k),points_Ox{j}(k),points_Oy{j}(k)]};
    end
    count=count+sols;
end
writecell(xcell,'C.xls')

figure(1)
plot(ex1,ey1,'or',ex2,ey2,'og',ex3,ey3,'ob')
axis equal
figure(2)
plot(Ox,Oy,'*k')
axis equal
