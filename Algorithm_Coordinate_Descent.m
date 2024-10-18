clearvars;close all;clc;
fun_name = 'Sixhump';
lower_bound = [-2,-2];
upper_bound = [2,2];
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
% coordinate descent
max_evaluation = 1000;
xmin = rand(1,2).*(upper_bound-lower_bound) + lower_bound;
fmin = feval(fun_name,xmin);
evaluation = 1;
xall = xmin;
fmin_record = fmin;
f_record = fmin;
fprintf('evaluation: %d, fmin: %0.2f\n',evaluation,fmin);
scatter(xall(:,1),xall(:,2),'go','filled');
scatter(xmin(:,1),xmin(:,2),'ro','filled');
while evaluation < max_evaluation
    % line search
    grid_n = 100;
    % search along x1
    x_search = [linspace(lower_bound(1),upper_bound(1),grid_n)',xmin(2)*ones(grid_n,1)];
    f_search = feval(fun_name,x_search);
    [fmin,index] = min(f_search);
    xmin = x_search(index,:);
    xall = [xall;x_search(index,:)];
    fmin_record = [fmin_record;fmin];
    scatter(xall(:,1),xall(:,2),'go','filled');
    scatter(xmin(:,1),xmin(:,2),'ro','filled');
    drawnow;
    pause;
    evaluation = evaluation+grid_n;
    fprintf('evaluation: %d, fmin: %0.2f\n',evaluation,fmin);
    % search along x2
    x_search = [xmin(1)*ones(grid_n,1),linspace(lower_bound(2),upper_bound(2),grid_n)'];
    f_search = feval(fun_name,x_search);
    [fmin,index] = min(f_search);
    xmin = x_search(index,:);
     xall = [xall;x_search(index,:)];
    fmin_record = [fmin_record;fmin];
    scatter(xall(:,1),xall(:,2),'go','filled');
    scatter(xmin(:,1),xmin(:,2),'ro','filled');
    drawnow;
    pause;
    evaluation = evaluation+grid_n;
    fprintf('evaluation: %d, fmin: %0.2f\n',evaluation,fmin);    
end
figure;
plot(f_record,'b-*'); hold on;
plot(fmin_record,'r-s');


