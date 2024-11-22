function [front_rank,rank_num] = NonDominated_Rank(pop_obj,num_sort)
% calculate the Pareto front ranks of each individual

[pop_obj,~,Loc] = unique(pop_obj,'rows');
Table          = hist(Loc,1:max(Loc));
[N,M]          = size(pop_obj);
[pop_obj,rank]  = sortrows(pop_obj);
front_rank        = inf(1,N);
rank_num         = 0;
while sum(Table(front_rank<inf)) < min(num_sort,length(Loc))
    rank_num = rank_num + 1;
    for i = 1 : N
        if front_rank(i) == inf
            Dominated = false;
            for j = i-1 : -1 : 1
                if front_rank(j) == rank_num
                    m = 2;
                    while m <= M && pop_obj(i,m) >= pop_obj(j,m)
                        m = m + 1;
                    end
                    Dominated = m > M;
                    if Dominated || M == 2
                        break;
                    end
                end
            end
            if ~Dominated
                front_rank(i) = rank_num;
            end
        end
    end
end
front_rank(rank) = front_rank;
front_rank       = front_rank(Loc');
front_rank = front_rank';
end

