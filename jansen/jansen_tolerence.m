function [max_pos_err,I,Fx,Fy]=jansen_tolerence(tolerence_list, angle_intervals, link_lengths)
    arguments
        tolerence_list (1,:) double
        angle_intervals (1,1) double
        link_lengths double =[38,41.5,39.3,40.1,55.8,39.4,36.7,65.7,49,50,61.9,7.8,15];
    end
L=getcondvects_n_k(length(tolerence_list),11,tolerence_list);
[r,c]=size(L);
Lprime=[zeros(r,1),L(:,1:10),zeros(r,1),L(:,11)];
thetas(1,1,:)=pi*[0:angle_intervals:360]/180;
E=zeros(2,ceil(360/angle_intervals)+1,r);
plotting=false;
if plotting==true
    figure(1)
    clf
    hold on
    axis equal
end
for i=1:r
    F=jansen_FK(thetas,Lprime(i,:)+link_lengths);
    E(:,:,i)=squeeze(F(:,9,:));
    if plotting==true && i<10
        plot(E(1,:,i),E(2,:,i))
    end
end
base=E(:,:,1);
[max_pos_err,I]=max(squeeze(sqrt((base(1,:)-E(1,:,:)).^2+(base(2,:)-E(2,:,:)).^2)),[],2);
Fx=squeeze(E(1,:,I));
Fy=squeeze(E(2,:,I));