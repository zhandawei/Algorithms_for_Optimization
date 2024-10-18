clearvars;close all;clc;
num_vari = 2;
fun_name = 'Regression_Loss';
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
% pattern search
max_evaluation = 100;
eita = 1;
xmin = rand(1,2).*(upper_bound-lower_bound) + lower_bound;
fmin = feval(fun_name,xmin);
fmin_record = fmin;        
evaluation = 1;
fprintf('evaluation: %d, fmin: %0.2f\n',evaluation,fmin);
for ii = 1:max_evaluation/4
    pattern_point = [xmin(1)+eita,xmin(2);xmin(1)-eita,xmin(2);xmin(1),xmin(2)+eita;xmin(1),xmin(2)-eita];
    pattern_y = feval(fun_name,pattern_point);
    evaluation = evaluation + 4;
    scatter(xmin(:,1),xmin(:,2),'ro','filled'); hold on;
    scatter(pattern_point(:,1),pattern_point(:,2),'go','filled');
    drawnow;
    pause;
    if min(pattern_y) < fmin
        [~,index] = min(pattern_y);
        xmin = pattern_point(index,:);
        fmin = min(pattern_y);
    else
        eita = 0.5*eita;
    end
    fmin_record = [fmin_record;fmin];
    fprintf('evaluation: %d, fmin: %0.2f\n',evaluation,fmin);
end
figure;
plot(fmin_record,'r-s');


