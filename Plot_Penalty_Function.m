clearvars;
lower_bound = [-2,-2];
upper_bound = [2,2];
grid_num = 101;
sigma = 0.1;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
f_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        f_mesh(ii,jj) = x(1) + x(2) + sigma*(x(1)^2 + x(2)^2-2)^2;
    end
end
figure;
subplot(1,2,1);
mesh(x1_mesh,x2_mesh,f_mesh);
title(sprintf('sigma = %g',sigma));
subplot(1,2,2);
contour(x1_mesh,x2_mesh,f_mesh,200);
axis equal;
hold on; 
scatter(-1,-1,'ro','filled');
[temp,index1] = min(f_mesh,[],2);
[min_f,index2] = min(temp);
xmin = [x1_mesh(index2,index1(index2)),x2_mesh(index2,index1(index2))];
scatter(xmin(1),xmin(2),'bo','filled');
title(sprintf('sigma = %g',sigma));