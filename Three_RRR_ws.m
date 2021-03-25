a=[100 100 100];
p=[100 100 100];
h=[100 100 100];
b=[0,0;200,0;100,100*sqrt(3)];
q=1;
points=[];
for i=0:30:359
    for j=0:30:359
        for k=0:30:359
            thetas=[i,j,k]*pi/180;
            [~,~,~,~,O]=Three_RRR_fk(a,p,h,b,thetas);
            if ~isempty(O)
                points=[points;O];
            end
        end
    end
end