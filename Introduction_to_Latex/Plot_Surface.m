clearvars; close all;
lower_bound = [0,0];
upper_bound = [5,5];
grid_num = 51;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
f_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        f_mesh(ii,jj) = Sasena(x);
    end
end
figure;
mesh(x1_mesh,x2_mesh,f_mesh);
xlabel('x_1');
ylabel('x_2');
zlabel('f');
% output data
output_x1_mesh = x1_mesh';
output_x2_mesh = x2_mesh';
output_f_mesh = f_mesh';
data_mesh = [output_x1_mesh(:) output_x2_mesh(:) output_f_mesh(:)];
writematrix(data_mesh,'f_mesh.dat','delimiter','\t');
