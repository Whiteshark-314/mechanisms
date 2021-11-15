function xldata= Three_RRR_tol_anla()
i=1;
xldata=cell(10,3);
while i<=10
    actuation_angles=randi(359,1,3);
    clear max_deltaO I
    fprintf('%f\n',i);
    [max_deltaO, I, ~]=Three_RRR_tol_max([0,-0.1,0.1],actuation_angles,1);
    if isnan(max_deltaO)
        continue;
    else
        xldata{i,1}=actuation_angles;
        xldata{i,2}=max_deltaO;
        xldata{i,3}=I;
        i=i+1;
    end
end
