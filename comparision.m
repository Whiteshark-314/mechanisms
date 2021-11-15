a=[100 100 100];
p=[100 100 100];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
Three_RRR_ws_t
Fx=F(:,1);
Fy=F(:,2);
Cx=sum(Fx)/3;
Cy=sum(Fy)/3;
t=[0:1:359]';
x_values=50*cosd(t)+Cx;
y_values=50*sind(t)+Cy;
alpha_values=zeros(length(t),1);
tol_val=0.05;
all_possible=getcondvects_n_k(3,9,[0,tol_val,-tol_val]);
all_possible_a=a+all_possible(:,7:9);
all_possible_p=p+all_possible(:,4:6);
all_possible_e=e+all_possible(:,1:3);
%load baseline.mat
deltaO=cell(1,19683);
meandeltaO=zeros(1,19683);
start_index=1;
stop_index=10;
parfor i=start_index:stop_index
    deviation=Three_RRR(all_possible_a(i,:),all_possible_p(i,:),all_possible_e(i,:),F);
    deviation=traj_repeat(deviation,baseline.traj_thetas_1(2:end,1),baseline.traj_thetas_2(2:end,1),baseline.traj_thetas_3(2:end,1),baseline.traj_alpha);
    dO=sqrt((x_values-deviation.traj_x).^2+(y_values-deviation.traj_y).^2);
    meandeltaO(i)=mean(deltaO);
end
filename=['meandeltaO',num2str(start_index),'_',num2str(stop_index),'.mat'];
save(filename,'meandeltaO')
fprintf(['DONE computing combinations ',num2str(start_index),' through ',num2str(stop_index),'!!! \n']);