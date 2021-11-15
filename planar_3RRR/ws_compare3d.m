a=[150 150 150];
p=[80 80 80];
e=[100 100 100];
F=[0,0;200,0;100,100*sqrt(3)];
Fx=F(:,1);
Fy=F(:,2);
baseline=Three_RRR(a,p,e,F);
baseline=workspace(baseline);
linktols=ones(1,9);
figure(3)
hold on
fill3(baseline.ws_x,baseline.ws_y,0*ones(1,length(baseline.ws_x)),b)
for i=0.1:0.1:1
    deviation=Three_RRR(a+i*linktols(1,7:9),p+i*linktols(1,4:6),e+i*linktols(1,1:3),F);
    deviation=workspace(deviation);
    fill3(deviation.ws_x,deviation.ws_y,-10*i*ones(1,length(baseline.ws_x)),[rand, rand, rand])
end
axis equal
