clearvars;close all;clc;
fun_name = 'Stochastic_Regression_Loss';
num_vari = 2;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
grid_num = 101;
[x1_mesh,x2_mesh] = meshgrid(linspace(lower_bound(1),upper_bound(1),grid_num),linspace(lower_bound(2),upper_bound(2),grid_num));
loss_mesh = zeros(grid_num,grid_num);
for ii = 1:grid_num
    for jj = 1:grid_num
        x = [x1_mesh(ii,jj),x2_mesh(ii,jj)];
        loss_mesh(ii,jj) = feval(fun_name,x);
    end
end
figure;
contour(x1_mesh,x2_mesh,loss_mesh,50);
axis equal;
hold on;
% gradient descent
max_evaluation = 100;
eita = 0.01;
x_old = rand(1,2).*(upper_bound-lower_bound) + lower_bound;
[f,df] = Stochastic_Regression_Loss(x_old);
xall = x_old;
fmin_record = f;
fprintf('evaluation: %d, fmin: %0.2f\n',1,f);
scatter(x_old(:,1),x_old(:,2),'ro','filled');
% pause(1);
for ii = 2:max_evaluation
    x_new = x_old - eita*df;
    [f,df] = Stochastic_Regression_Loss(x_new);
    xall = [xall;x_new];
    x_old = x_new;
    fmin_record = [fmin_record;f];
    scatter(xall(:,1),xall(:,2),'go','filled');
    scatter(x_new(:,1),x_new(:,2),'ro','filled');
    drawnow;
    % pause(1);
    fprintf('evaluation: %d, fmin: %0.2f\n',ii,f);
end
figure;
plot(fmin_record,'r-s');


