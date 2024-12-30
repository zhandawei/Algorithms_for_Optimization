clearvars;
lower_bound = [-3,-3];
upper_bound = [3,3];
grid_num = 101;
sigma = 100;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
f_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        f_mesh(ii,jj) = norm(x) + sin(4*atan(x(2)/x(1))) - 1/(sigma*(-x(1)^2-x(2)^2+2));
    end
end
figure;
subplot(1,2,1);
mesh(x1_mesh,x2_mesh,f_mesh);
title(sprintf('sigma = %g',sigma));
subplot(1,2,2);
contour(x1_mesh,x2_mesh,f_mesh,100);
axis equal;
hold on;
