count=0;
B=[];
C=[];
for i=1:30-2
    B=[B sum(A(i:i+2))];
    if i>1
        if B(i)>B(i-1)
            count=count+1;
        end
    end
end
