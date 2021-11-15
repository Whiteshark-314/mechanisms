a=[100 100 100];
p=[100 100 100];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
Fx=F(:,1);
Fy=F(:,2);
tol_val=[-0.1,0,0.1];
all_possible=getcondvects_n_k(3,9,tol_val);
all_possible_a=a+all_possible(:,7:9);
all_possible_p=p+all_possible(:,4:6);
all_possible_e=e+all_possible(:,1:3);
baseline=Three_RRR(a,p,e,F);
baseline=workspace(baseline);
R=baseline.r;
Area=polyarea(baseline.ws_x,baseline.ws_y);
baseline=trajectory(baseline,baseline.ws_x*0.5,baseline.ws_y*0.5,zeros(1,length(baseline.ws_x)));
load baseline.mat
listdeltaR_1=zeros(1,19683);
listdeltaR_2=zeros(1,19683);
listdeltaR_3=zeros(1,19683);
dArea=zeros(1,19683);
start_index=1;
stop_index=19683;
parfor i=start_index:stop_index
    deviation=Three_RRR(all_possible_a(i,:),all_possible_p(i,:),all_possible_e(i,:),F);
    deviation=workspace(deviation);
    dR=abs(R-deviation.r);
    listdeltaR_1(i)=dR(1);
    listdeltaR_2(i)=dR(2);
    listdeltaR_3(i)=dR(3);
    dArea(i)=polyarea(deviation.ws_x,deviation.ws_y)-Area
end
listdeltaR=[listdeltaR_1;listdeltaR_2;listdeltaR_3];
[sorted_diff_area,idx]=sort(dArea,'descend');
max_deviation=Three_RRR(all_possible_a(idx(1),:),all_possible_p(idx(1),:),all_possible_e(idx(1),:),F);
max_deviation=workspace(max_deviation);
min_deviation=Three_RRR(all_possible_a(idx(end),:),all_possible_p(idx(end),:),all_possible_e(idx(end),:),F);
min_deviation=workspace(min_deviation);
plot(baseline.ws_x,baseline.ws_y,max_deviation.ws_x,max_deviation.ws_y)%,min_deviation.ws_x,min_deviation.ws_y)
axis equal