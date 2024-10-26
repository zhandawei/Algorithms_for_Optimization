clearvars;close all;clc;
fun_name = 'Rastrigin';
num_vari = 10;
max_evaluation = 10000;
max_run = 10;
[lower_bound,upper_bound] = Test_Function(fun_name,num_vari);
% random search
fmin_record_RS = zeros(max_evaluation,max_run);
r = 0.01*norm(upper_bound - lower_bound);
for run  = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    fmin_record_RS(1,run) = fmin;
    fprintf('Random Search on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,1,fmin);
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
        fmin_record_RS(ii,run) = fmin;
        fprintf('random search on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,ii,fmin);
    end
end


% coordinate descent
grid_n = 100;
fmin_record_CD = zeros(max_evaluation/grid_n+1,max_run);
for run = 1:max_run
    xmin = rand(1,num_vari).*(upper_bound-lower_bound) + lower_bound;
    fmin = feval(fun_name,xmin);
    evaluation = 1;
    iteration = 1;
    fmin_record_CD(iteration,run) = fmin;
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
            iteration = iteration +1;
            fmin_record_CD(iteration,run) = fmin;
            fprintf('coordinate descent on %s problem, run: %d, evaluation: %d, fmin: %0.2f\n',fun_name,run,evaluation,fmin);
        end
    end
end

% pattern search
fmin_record_PS = zeros(max_evaluation/(2*num_vari)+1,max_run);
eita = 0.1*norm(upper_bound - lower_bound);
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



% simulated annealing
r = 0.01*norm(upper_bound - lower_bound);
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

% genetic algorithm
pop_size = 200;
max_gen = max_evaluation/pop_size;
fmin_record_GA = zeros(max_evaluation/pop_size,max_run);
for run = 1:max_run
    generation = 1;
    % the initial generation
    pop_size = pop_size - rem(pop_size,2);
    pop_vari = lhsdesign(pop_size, num_vari).*(upper_bound - lower_bound) + lower_bound;
    % calculate the objective values
    pop_fitness = feval(fun_name,pop_vari);
    fmin_record_GA(generation,run) = min(pop_fitness);
    evaluation = pop_size;
    % the evoluation of the generation
    while generation < max_gen
        % parent selection using k-tournament (default k=2) selection
        k = 2;
        temp = randi(pop_size,pop_size,k);
        [~,index] = min(pop_fitness(temp),[],2);
        pop_parent = pop_vari(sum(temp.*(index == 1:k),2),:);
        % crossover (simulated binary crossover)
        % dic_c is the distribution index of crossover
        % crossover rate is 0.9
        dis_c = 20;
        pro_c = 0.9;
        mu  = rand(pop_size/2,num_vari);
        parent1 = pop_parent(1:2:pop_size,:);
        parent2 = pop_parent(2:2:pop_size,:);
        beta = 1 + 2*min(min(parent1,parent2)-lower_bound,upper_bound-max(parent1,parent2))./max(abs(parent2-parent1),1E-6);
        alpha = 2 - beta.^(-dis_c-1);
        betaq = (alpha.*mu).^(1/(dis_c+1)).*(mu <= 1./alpha) + (1./(2-alpha.*mu)).^(1/(dis_c+1)).*(mu > 1./alpha);
        % the crossover is performed randomly on each variable
        betaq = betaq.*(-1).^randi([0,1],pop_size/2,num_vari);
        betaq(rand(pop_size/2,num_vari)>pro_c) = 1;
        offspring1 = 0.5*((1+betaq).*parent1 + (1-betaq).*parent2);
        offspring2 = 0.5*((1-betaq).*parent1 + (1+betaq).*parent2);
        pop_crossover = [offspring1;offspring2];
        % mutation (polynomial mutation)
        % dis_m is the distribution index of polynomial mutation
        % mutation rate is 1/d
        dis_m = 20;
        pro_m = 1/num_vari;
        rand_var = rand(pop_size,num_vari);
        mu  = rand(pop_size,num_vari);
        deta = min(pop_crossover-lower_bound, upper_bound-pop_crossover)./(upper_bound-lower_bound);
        detaq = zeros(pop_size,num_vari);
        position1 = rand_var<=pro_m & mu<=0.5;
        position2 = rand_var<=pro_m & mu>0.5;
        detaq(position1) = ((2*mu(position1) + (1-2*mu(position1)).*(1-deta(position1)).^(dis_m+1)).^(1/(dis_m+1))-1);
        detaq(position2) = (1 - (2*(1-mu(position2))+2*(mu(position2)-0.5).*(1-deta(position2)).^(dis_m+1)).^(1/(dis_m+1)));
        pop_mutation = pop_crossover + detaq.*(upper_bound-lower_bound);
        pop_mutation  = max(min(pop_mutation,upper_bound),lower_bound);
        % fitness calculation
        pop_mutation_fitness = feval(fun_name,pop_mutation);
        evaluation = evaluation + pop_size;
        % environment selection
        pop_vari_iter = [pop_vari;pop_mutation];
        pop_fitness_iter = [pop_fitness;pop_mutation_fitness];
        [~,win_num] = sort(pop_fitness_iter);
        pop_vari = pop_vari_iter(win_num(1:pop_size),:);
        pop_fitness = pop_fitness_iter(win_num(1:pop_size),:);
        % update the evaluation number of generation number
        generation = generation + 1;
        fmin_record_GA(generation,run) = min(pop_fitness);
        fprintf('Genetic Algorithm on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,run,evaluation,min(pop_fitness));
    end
end

% plot convergence figure
figure;
plot(1:1:max_evaluation,prctile(fmin_record_RS,50,2),'k');hold on;
plot(0:grid_n:max_evaluation,prctile(fmin_record_CD,50,2),'g');
plot(0:2*num_vari:max_evaluation,prctile(fmin_record_PS,50,2),'b');
plot(0:2*num_vari:max_evaluation,prctile(fmin_record_improved,50,2),'m');
plot(1:1:max_evaluation,prctile(fmin_record_SA,50,2),'c');
plot(pop_size:pop_size:max_evaluation,prctile(fmin_record_GA,50,2),'r');
title(strcat(num2str(num_vari),'D-',fun_name));

lower_bound_RS = prctile(fmin_record_RS,25,2);
upper_bound_RS = prctile(fmin_record_RS,75,2);
lower_bound_CD = prctile(fmin_record_CD,25,2);
upper_bound_CD = prctile(fmin_record_CD,75,2);
lower_bound_PS = prctile(fmin_record_PS,25,2);
upper_bound_PS = prctile(fmin_record_PS,75,2);
lower_bound_improved = prctile(fmin_record_improved,25,2); 
upper_bound_improved = prctile(fmin_record_improved,75,2);
lower_bound_SA = prctile(fmin_record_SA,25,2);
upper_bound_SA = prctile(fmin_record_SA,75,2);
lower_bound_GA = prctile(fmin_record_GA,25,2);
upper_bound_GA = prctile(fmin_record_GA,75,2);
h = fill([(1:1:max_evaluation)';(max_evaluation:-1:1)'],[lower_bound_RS;upper_bound_RS(end:-1:1)],'k');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(0:grid_n:max_evaluation)';(max_evaluation:-grid_n:0)'],[lower_bound_CD;upper_bound_CD(end:-1:1)],'g');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(0:2*num_vari:max_evaluation)';(max_evaluation:-2*num_vari:0)'],[lower_bound_PS;upper_bound_PS(end:-1:1)],'b');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(0:2*num_vari:max_evaluation)';(max_evaluation:-2*num_vari:0)'],[lower_bound_improved;upper_bound_improved(end:-1:1)],'m');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(1:1:max_evaluation)';(max_evaluation:-1:1)'],[lower_bound_SA;upper_bound_SA(end:-1:1)],'c');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(pop_size:pop_size:max_evaluation)';(max_evaluation:-pop_size:pop_size)'],[lower_bound_GA;upper_bound_GA(end:-1:1)],'r');
set(h,'edgealpha',0,'facealpha',0.3);
legend('Random Search','Coordinate Descent','Pattern Search','Improved PS','Simulated Annealing','Genetic Algorithm');
xlabel('number of evaluations');
ylabel('current minimum objective function');





