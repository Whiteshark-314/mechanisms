cats = categories(stormData.State);
L=length(cats);
costStatewise=table;
costStatewise.State = categorical(zeros(L,1));
costStatewise.PropertyCost = zeros(L,1);
for i=1:L
costStatewise.State(i) = cats{i};
costStatewise.PropertyCost(i) = sum(stormData.Property_Cost(stormData.State==cats{i}), 'omitnan');
end
costStatewise = sortrows(costStatewise, 'PropertyCost', 'descend');
head(costStatewise, 5)