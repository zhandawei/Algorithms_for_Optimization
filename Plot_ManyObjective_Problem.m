clearvars;clc;close all;
n = 1000;
num_nondominated_point = zeros(1,20);
for m = 1:20
    rand_point = rand(n,m);
    non_dominated_point = Pareto_Set(rand_point);
    num_nondominated_point(m) = size(non_dominated_point,1);
end
plot(1:20,num_nondominated_point,'r-o')
