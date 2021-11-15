tol=0.001;
angular_interval=10;
tolerence_list=[0,tol,-tol];
[max_pos_err,I,Fx,Fy]=jansen_tolerence(tolerence_list, angular_interval);
L=getcondvects_n_k(length(tolerence_list),11,tolerence_list);
[r,c]=size(L);
Lprime=[zeros(r,1),L(:,1:10),zeros(r,1),L(:,11)];
thetas=pi*[0:angular_interval:360]/180;
xcell=cell(length(thetas)+1,6);
xcell(1,:)={'Angle','Combination','Combination_values','Fx','Fy','Max_pos_error'};
for i=1:37
    xcell(i+1,:)={thetas(i)*180/pi,I(i),Lprime(I(i),:),Fx(i,i),Fy(i,i),max_pos_err(i)};
end
writecell(xcell,['jansen_tol_',num2str(tol),'_',num2str(angular_interval),'.xls'])