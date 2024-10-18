clearvars;close all;clc;
fun_name = 'Ackley';
lower_bound = [-10,-10];
upper_bound = [10,10];
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
% random search
max_evaluation = 100;
r = 1;
% xmin = rand(1,2).*(upper_bound-lower_bound) + lower_bound;
xmin = [-8,-8];
fmin = feval(fun_name,xmin);
xall = xmin;
fmin_record = fmin;
f_record = fmin;
fprintf('evaluation: %d, fmin: %0.2f\n',1,fmin);
scatter(xall(:,1),xall(:,2),'go','filled');
scatter(xmin(:,1),xmin(:,2),'ro','filled');
for ii = 2:max_evaluation
    rand_vector = randn(1,2);
    y = rand_vector/sqrt(sum(rand_vector.^2));
    u = rand(1,1);
    xnew = xmin + r*y*u^(1/2);
    while any(xnew<lower_bound) || any(xnew>upper_bound)
        rand_vector = randn(1,2);
        y = rand_vector/sqrt(sum(rand_vector.^2));
        u = rand(1,1);
        xnew = xmin + r*y*u^(1/2);
    end
    fnew = feval(fun_name,xnew);
    xall = [xall;xnew];
    if fnew < fmin
        xmin = xnew;
        fmin = fnew;
    end
    fmin_record = [fmin_record;fmin];
    f_record = [f_record;fnew];
    scatter(xall(:,1),xall(:,2),'go','filled');
    scatter(xmin(:,1),xmin(:,2),'ro','filled');
    drawnow;
    %pause;
    fprintf('evaluation: %d, fmin: %0.2f\n',ii,fmin);
end
figure;
plot(f_record,'b-*'); hold on;
plot(fmin_record,'r-s');


