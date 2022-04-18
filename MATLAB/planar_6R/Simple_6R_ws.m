a=[150 150 150];
p=[80 80 80];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
inc=2;
n=length(-180:inc:180);
points_ex1=cell(n^3,1);
points_ex2=cell(n^3,1);
points_ex3=cell(n^3,1);
points_ey1=cell(n^3,1);
points_ey2=cell(n^3,1);
points_ey3=cell(n^3,1);
points_Ox=cell(n^3,1);
points_Oy=cell(n^3,1);
condvects=getcondvects_n_k(n,3,-180:inc:180);
parfor i=1:size(condvects,1)
    thetas=condvects(i,:)*pi/180;
    s6r=Three_RRR(a,p,e,F);
    s6r=fk(s6r,thetas);
    if ~isempty(s6r.O)
%         abc=a(1)+p(1)+sqrt(sum((e_xy(1,:)-O_xy).^2));
%         fgh=sqrt(mech.O(:,1).^2+mech.O(:,2).^2);
%         if sum(fgh>abc)>0
%             fprintf("Check FK for %f %f %f\n",thetas)
%         end
        points_Ox{i}= s6r.O(:,1);
        points_Oy{i}= s6r.O(:,2);
        points_ex1{i}=s6r.Ex(:,1);
        points_ex2{i}=s6r.Ex(:,2);
        points_ex3{i}=s6r.Ex(:,3);
        points_ey1{i}=s6r.Ey(:,1);
        points_ey2{i}=s6r.Ey(:,2);
        points_ey3{i}=s6r.Ey(:,3);
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
O_x=cell2mat(points_Ox);
O_y=cell2mat(points_Oy);
figure(3)
% Three_RRR_ws_t
% %plot(ex1,ey1,'ob',ex2,ey2,'ob',ex3,ey3,'ob')
% hold on
scatter(O_x,O_y,'.b')
% F_prime=[F;F(1,:)];
% plot(F_prime(:,1),F_prime(:,2),'k')
axis equal
% hold off