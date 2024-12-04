clearvars;close all;
load('3D_data.mat');
figure;
scatter3(pareto_front(:,1),pareto_front(:,2),pareto_front(:,3),'b.'); hold on;
scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
writematrix(pareto_front,'pareto_front.dat','delimiter','\t');
writematrix(non_dominated_front,'non_dominated_front.dat','delimiter','\t');