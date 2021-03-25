a=[100 100 100];
p=[100 100 100];
h=[100 100 100];
b=[0,0;200,0;100,100*sqrt(3)];
t=[65 130 115]; 
tt=[0 2*pi/3 4*pi/3];
%thetas=[0 2*pi/3 4*pi/3];
thetas=t*pi/180;
[Phi,alpha,Ex,Ey,O]=Three_RRR_fk(a,p,h,b,thetas);
if ~isempty(alpha)
    figure(2)
    l=length(alpha);
    tri=zeros(4,2);
    if l==2
        m=1;n=2;
    elseif l==4
        m=2;n=2;
    elseif l==6
        m=2;n=3;
    elseif l==8
        m=3;n=3;
    end
    for i=1:length(alpha)
        subplot(m,n,i)
        plot(b(:,1),b(:,2),'o')
        hold on
        axis equal
        for j=1:3
            a_l=[b(j,:);b(j,1)+a(j)*cos(thetas(j)) b(j,2)+a(j)*sin(thetas(j))];
            plot(a_l(:,1),a_l(:,2))
            p_l=[a_l(2,:);a_l(2,1)+p(j)*cos(Phi(j,i)) a_l(2,2)+p(j)*sin(Phi(j,i))];
            plot(p_l(:,1),p_l(:,2))
            %e_l=[p_l(2,:);p_l(2,1)+h(j)*cos(alpha(i)+tt(j)) p_l(2,2)+h(j)*sin(alpha(i)+tt(j))];
            %plot(e_l(:,1),e_l(:,2))
            tri(j,:)=p_l(2,:);
        end
        tri(4,:)=tri(1,:);
        plot(tri(:,1),tri(:,2))
        plot(O(i,1),O(i,2),'o')
        hold off
    end
end
hold off