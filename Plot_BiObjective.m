clearvars;close all;
lower_bound = [0,0];
upper_bound = [3,3];
grid_num = 101;
x1_line = linspace(lower_bound(1),upper_bound(1),grid_num);
x2_line = linspace(lower_bound(2),upper_bound(2),grid_num);
x_point = zeros(grid_num*grid_num,2);
f_point = zeros(grid_num*grid_num,2);
for ii = 1:grid_num
    for jj = 1:grid_num
        x_point((ii-1)*grid_num+jj,:) = [x1_line(ii),x2_line(jj)];
        f_point((ii-1)*grid_num+jj,:) = Bi_Objective(x_point((ii-1)*grid_num+jj,:));
    end
end
figure;
subplot(1,2,1);
scatter(x_point(:,1),x_point(:,2),'.');
subplot(1,2,2);
scatter(f_point(:,1),f_point(:,2),'.');
