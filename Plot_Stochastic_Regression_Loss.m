clearvars;close all;
lower_bound = [-10,-10];
upper_bound = [10,10];
grid_num = 101;
batch_size = 5;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
loss_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        loss_mesh(ii,jj) = Stochastic_Regression_Loss(x,batch_size);
    end
end
figure;
subplot(1,2,1);
mesh(x1_mesh,x2_mesh,loss_mesh);
subplot(1,2,2);
contour(x1_mesh,x2_mesh,loss_mesh,50);
axis equal;
hold on;
[min1,index1] = min(loss_mesh);
[min2,index2] = min(min1);
xmin = [x1_mesh(index1(index2),index2),x2_mesh(index1(index2),index2)];
scatter(xmin(1),xmin(2),'r','filled');



