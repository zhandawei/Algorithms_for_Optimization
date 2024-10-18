clearvars;close all;clc;
fun_name = 'Ackley';
num_vari = 100;
max_evaluation = 10000;
max_run = 10;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
% random search
r = 0.01*norm(upper_bound - lower_bound);
fmin_record_RS= zeros(max_evaluation,max_run);
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    iteration = 1;
    evaluation = 1;
    fmin_record_RS(iteration,run) = fmin;
    while evaluation < max_evaluation
        rand_vector = randn(1,num_vari);
        y = rand_vector/sqrt(sum(rand_vector.^2));
        xnew = xmin + r*y*rand^(1/num_vari);
        while any(xnew<lower_bound) || any(xnew>upper_bound)
            rand_vector = randn(1,num_vari);
            y = rand_vector/sqrt(sum(rand_vector.^2));
            xnew = xmin + r*y*rand^(1/num_vari);
        end
        fnew = feval(fun_name,xnew);
        iteration = iteration+1;
        evaluation = evaluation+1;
        if fnew < fmin
            xmin = xnew;
            fmin = fnew;
        end
        fmin_record_RS(iteration,run) = fmin;
        fprintf('Random Search on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,run,evaluation,fmin);
    end
end
% simulated annealing
fmin_record_SA = zeros(max_evaluation,max_run);
T0 = 1000;
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    ymin = fmin;
    iteration = 1;
    evaluation = 1;
    fmin_record_SA(iteration,run) = ymin;
    while evaluation < max_evaluation
        rand_vector = randn(1,num_vari);
        y = rand_vector/sqrt(sum(rand_vector.^2));
        xnew = xmin + r*y*rand^(1/num_vari);
        while any(xnew<lower_bound) || any(xnew>upper_bound)
            rand_vector = randn(1,num_vari);
            y = rand_vector/sqrt(sum(rand_vector.^2));
            xnew = xmin + r*y*rand^(1/num_vari);
        end
        fnew = feval(fun_name,xnew);
        ymin = min(ymin,fnew);
        iteration = iteration+1;
        evaluation = evaluation+1;
        if fnew < fmin
            xmin = xnew;
            fmin = fnew;
        else
            p = 1/(1+exp((fnew-fmin)/(T0*0.999^iteration)));
            if rand < p
                xmin = xnew;
                fmin = fnew;
            end
        end
        fmin_record_SA(iteration,run) = ymin;
        fprintf('Simulated Annealing on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,run,evaluation,ymin);
    end
end
% plot convergence figure
figure;
plot(1:1:max_evaluation,prctile(fmin_record_RS,50,2),'k');hold on;
plot(1:1:max_evaluation,prctile(fmin_record_SA,50,2),'r');
title(strcat(num2str(num_vari),'D-',fun_name));
lower_bound_PS = prctile(fmin_record_RS,25,2);
upper_bound_PS = prctile(fmin_record_RS,75,2);
lower_bound_SA = prctile(fmin_record_SA,25,2);
upper_bound_SA = prctile(fmin_record_SA,75,2);
h = fill([(1:1:max_evaluation)';(max_evaluation:-1:1)'],[lower_bound_PS;upper_bound_PS(end:-1:1)],'k');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(1:1:max_evaluation)';(max_evaluation:-1:1)'],[lower_bound_SA;upper_bound_SA(end:-1:1)],'r');
set(h,'edgealpha',0,'facealpha',0.3);
legend('Random Search','Simulated Annealing');
xlabel('number of evaluations');
ylabel('current minimum objective function');





