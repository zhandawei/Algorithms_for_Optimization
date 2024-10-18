clearvars;close all;
lower_bound = [-10,-10];
upper_bound = [10,10];
grid_num = 101;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
loss_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        point = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        loss_mesh(ii,jj) = Regression_Loss(point);
    end
end
loss2_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        point = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        loss2_mesh(ii,jj) = Regression_Loss(point) + 50*sum(point.^2);
    end
end
figure;
subplot(1,2,1);
contour(x1_mesh,x2_mesh,loss_mesh,50);
axis equal;
subplot(1,2,2);
contour(x1_mesh,x2_mesh,loss2_mesh,50);
axis equal;


