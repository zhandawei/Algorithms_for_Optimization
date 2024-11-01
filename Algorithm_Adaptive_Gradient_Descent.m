clearvars;close all;clc;
fun_name = 'Valley_Function';
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
max_evaluation = 20;
eita = 0.1;
x_old = [-8,-8];
[f,df] = feval(fun_name,x_old);
G = df;
xall = x_old;
fmin_record = f;
fprintf('evaluation: %d, fmin: %0.2f\n',1,f);
scatter(x_old(:,1),x_old(:,2),'ro','filled');
pause(0.5);
for ii = 2:max_evaluation
    x_new = x_old - eita*df./(sqrt(sum(G.^2,1))/sum(sqrt(sum(G.^2))));
    [f,df] = feval(fun_name,x_new);
    xall = [xall;x_new];
    G = [G;df];
    x_old = x_new;
    fmin_record = [fmin_record;f];
    plot(xall(:,1),xall(:,2),'b-');
    scatter(xall(:,1),xall(:,2),'bo','filled');
    scatter(x_new(:,1),x_new(:,2),'ro','filled');
    drawnow;
    pause(0.5);
    fprintf('evaluation: %d, fmin: %0.2f\n',ii,f);
end
figure;
plot(fmin_record,'r-s');


