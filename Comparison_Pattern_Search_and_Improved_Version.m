clearvars;close all;clc;
fun_name = 'Ackley';
num_vari = 100;
max_evaluation = 10000;
max_run = 100;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
% pattern search
fmin_record_PS = zeros(max_evaluation/(2*num_vari)+1,max_run);
eita = 1;
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    iteration = 1;
    evaluation = 0;
    fmin_record_PS(iteration,run) = fmin;
    while evaluation < max_evaluation
        pattern_point = ones(2*num_vari,1)*xmin;
        for ii = 1:num_vari
            pattern_point((ii-1)*2+1,ii) = pattern_point((ii-1)*2+1,ii) + eita;
            pattern_point(ii*2,ii) = pattern_point(ii*2,ii) - eita;
        end
        pattern_y = feval(fun_name,pattern_point);
        iteration = iteration + 1;
        evaluation = evaluation + 2*num_vari;
        if min(pattern_y) < fmin
            [~,index] = min(pattern_y);
            xmin = pattern_point(index,:);
            fmin = min(pattern_y);
        else
            eita = 0.5*eita;
        end
        fmin_record_PS(iteration,run) = fmin;
        fprintf('Pattern Search on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,run,evaluation,fmin);
    end
end
% improved version
fmin_record_improved = zeros(max_evaluation/(2*num_vari)+1,max_run);
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    iteration = 1;
    evaluation = 0;
    fmin_record_improved(iteration,run) = fmin;
    while evaluation < max_evaluation
        temp = randn(2*num_vari,num_vari);
        direction = temp./sqrt(sum(temp.^2,2));
        eita_max = min(max((upper_bound-xmin)./direction,(lower_bound-xmin)./direction),[],2);
        pattern_point = xmin + rand(2*num_vari,1).*eita_max.*direction;
        pattern_y = feval(fun_name,pattern_point);
        iteration = iteration + 1;
        evaluation = evaluation + 2*num_vari;
        if min(pattern_y) < fmin
            [~,index] = min(pattern_y);
            xmin = pattern_point(index,:);
            fmin = min(pattern_y);
        end
        fmin_record_improved(iteration,run) = fmin;
        fprintf('Improved Pattern Search on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,run,evaluation,fmin);
    end
end
% plot convergence figure
figure;
plot(0:2*num_vari:max_evaluation,prctile(fmin_record_PS,50,2),'k');hold on;
plot(0:2*num_vari:max_evaluation,prctile(fmin_record_improved,50,2),'r');
title(strcat(num2str(num_vari),'D-',fun_name));
lower_bound_PS = prctile(fmin_record_PS,25,2);
upper_bound_PS = prctile(fmin_record_PS,75,2);
lower_bound_improved = prctile(fmin_record_improved,25,2);
upper_bound_improved = prctile(fmin_record_improved,75,2);
h = fill([(0:2*num_vari:max_evaluation)';(max_evaluation:-2*num_vari:0)'],[lower_bound_PS;upper_bound_PS(end:-1:1)],'k');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(0:2*num_vari:max_evaluation)';(max_evaluation:-2*num_vari:0)'],[lower_bound_improved;upper_bound_improved(end:-1:1)],'r');
set(h,'edgealpha',0,'facealpha',0.3);
legend('Pattern Search','Improved Pattern Search');
xlabel('number of evaluations');
ylabel('current minimum objective function');





