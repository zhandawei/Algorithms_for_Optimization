clearvars;close all;clc;
num_vari = 100;
max_evaluation = 100000;
fun_name = 'Rosenbrock';
max_run = 1;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
% random search
fmin_record = zeros(max_evaluation,max_run);
for run  = 1:max_run
    r = 1;
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    fmin_record(1,run) = fmin;
    fprintf('random search on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,1,fmin);
    for ii = 2:max_evaluation
        rand_vector = randn(1,num_vari);
        y = rand_vector/sqrt(sum(rand_vector.^2));
        u = rand(1,1);
        xnew = xmin + r*y*u^(1/num_vari);
        while any(xnew<lower_bound) || any(xnew>upper_bound)
            rand_vector = randn(1,num_vari);
            y = rand_vector/sqrt(sum(rand_vector.^2));
            u = rand(1,1);
            xnew = xmin + r*y*u^(1/num_vari);
        end
        fnew = feval(fun_name,xnew);
        if fnew < fmin
            xmin = xnew;
            fmin = fnew;
        end
        fmin_record(ii,run) = fmin;
        fprintf('random search on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,ii,fmin);
    end
end


% coordinate descent
grid_n = 100;
fmin_record2 = zeros(max_evaluation+1,max_run);
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    evaluation = 1;
    fmin_record2(evaluation,run) = fmin;
    fprintf('coordinate descent on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,evaluation,fmin);
    while evaluation < max_evaluation
        % line search
        for ii = 1:num_vari
            x_search = repmat(xmin,grid_n,1);
            x_search(:,ii) = linspace(lower_bound(ii),upper_bound(ii),grid_n)';
            f_search = feval(fun_name,x_search);
            [fmin,index] = min(f_search);
            xmin = x_search(index,:);            
            evaluation = evaluation+grid_n;
            fmin_record2(evaluation,run) = fmin;
            fprintf('coordinate descent on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,evaluation,fmin);
        end
    end
end




figure;
plot(1:1:max_evaluation,mean(fmin_record,2),'b-'); hold on;
plot(1:grid_n:max_evaluation+1,mean(fmin_record2(1:grid_n:max_evaluation+1,:),2),'r-');
legend('random search','coordinate descent')