clearvars;close all;
lower_bound = [-10,-10];
upper_bound = [10,10];
grid_num = 101;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
loss_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        loss_mesh(ii,jj) = 55*x(1)^2+5*x(2)^2+30*x(1)*x(2)-315*x(1)-91.4*x(2)+461.57 + 1*sum(x.^2);
    end
end
figure;
subplot(1,2,1);
mesh(x1_mesh,x2_mesh,loss_mesh);
subplot(1,2,2);
contour(x1_mesh,x2_mesh,loss_mesh,50);
axis equal;
hold on;
scatter(0,0,'b','filled');
scatter(2,3.2,'b','filled');
[min1,index1] = min(loss_mesh);
[min2,index2] = min(min1);
xmin = [x1_mesh(index1(index2),index2),x2_mesh(index1(index2),index2)];
scatter(xmin(1),xmin(2),'r','filled');



