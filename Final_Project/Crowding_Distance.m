function crowding_distance = Crowding_Distance(pop_obj,front_rank)

[n,m]    = size(pop_obj);
crowding_distance = zeros(n,1);
fronts   = setdiff(unique(front_rank),inf);
for f = 1 : length(fronts)
    front = find(front_rank==fronts(f));
    fmax  = max(pop_obj(front,:),[],1);
    fmin  = min(pop_obj(front,:),[],1);
    for i = 1 : m
        [~,rank] = sortrows(pop_obj(front,i));
        crowding_distance(front(rank(1)))   = inf;
        crowding_distance(front(rank(end))) = inf;
        for j = 2 : length(front)-1
            crowding_distance(front(rank(j))) = crowding_distance(front(rank(j)))+(pop_obj(front(rank(j+1)),i)-pop_obj(front(rank(j-1)),i))/(fmax(i)-fmin(i));
        end
    end
end
end

