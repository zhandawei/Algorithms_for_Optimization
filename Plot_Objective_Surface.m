clearvars;close all;
fun_name = 'Ackley';
num_vari = 2;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
grid_num = 101;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
f_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        point = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        f_mesh(ii,jj) = feval(fun_name,point);
    end
end
figure;
subplot(1,2,1);
mesh(x1_mesh,x2_mesh,f_mesh);
subplot(1,2,2);
contour(x1_mesh,x2_mesh,f_mesh,50);