clearvars;clc;close all;
fun_name = 'DTLZ2';
num_obj = 3;
num_vari = 10;
pop_size = 20;
max_evaluation = 200;
lower_bound = zeros(1,num_vari);
upper_bound = ones(1,num_vari);
pareto_front = Calculate_Pareto_Front(fun_name,10000,num_obj);
pareto_front_plot = Calculate_Pareto_Front(fun_name,1000,num_obj);
evaluation = 0;
generation = 1;
pop_vari = rand(pop_size,num_vari).*(upper_bound-lower_bound) + lower_bound;
pop_obj = feval(fun_name, pop_vari, num_obj);
evaluation = evaluation + size(pop_obj,1);
non_dominated_front = Pareto_Set(pop_obj);
IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
if num_obj == 2
    scatter(pop_obj(:,1),pop_obj(:,2),'ro');hold on;
    scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
    scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
    title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    drawnow;hold off;
elseif num_obj == 3
    scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
    scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
    scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
    title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    view(135,30);drawnow;hold off;
end
while evaluation < max_evaluation
    % tournament selection
    front_rank = NonDominated_Rank(pop_obj,pop_size);
    crowd_distance = Crowding_Distance(pop_obj,front_rank);
    k = 2;
    [~,rank] = sortrows([front_rank,-crowd_distance]);
    [~,rank] = sort(rank);
    rand_index  = randi(size(pop_vari,1),k,pop_size);
    [~,best] = min(rank(rand_index),[],1);
    index    = rand_index(best+(0:pop_size-1)*k)';
    parent = pop_vari(index,:);
    % simulated binary crossover
    pro_c = 1;
    dis_c = 20;
    parent    = parent([1:size(parent,1),1:ceil(size(parent,1)/2)*2-size(parent,1)]',:);
    parent1 = parent(1:pop_size/2,:);
    parent2 = parent(pop_size/2+1:end,:);
    beta = zeros(pop_size/2,num_vari);
    mu   = rand(pop_size/2,num_vari);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(dis_c+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(dis_c+1));
    beta = beta.*(-1).^randi([0,1],pop_size/2,num_vari);
    beta(rand(pop_size/2,num_vari)<0.5) = 1;
    beta(repmat(rand(pop_size/2,1)>pro_c,1,num_vari)) = 1;
    pop_vari_crossover = [(parent1+parent2)/2+beta.*(parent1-parent2)/2;(parent1+parent2)/2-beta.*(parent1-parent2)/2];
    % polynomial mutation
    pro_m = 1;
    dis_m = 20;
    pop_vari_mutation = pop_vari_crossover;
    lower = repmat(lower_bound,pop_size,1);
    upper = repmat(upper_bound,pop_size,1);
    site  = rand(pop_size,num_vari) < pro_m/num_vari;
    mu    = rand(pop_size,num_vari);
    temp  = site & mu<=0.5;
    pop_vari_mutation(temp) = pop_vari_mutation(temp)+(upper(temp)-lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(pop_vari_mutation(temp)-lower(temp))./(upper(temp)-lower(temp))).^(dis_m+1)).^(1/(dis_m+1))-1);
    temp = site & mu>0.5;
    pop_vari_mutation(temp) = pop_vari_mutation(temp)+(upper(temp)-lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(upper(temp)-pop_vari_mutation(temp))./(upper(temp)-lower(temp))).^(dis_m+1)).^(1/(dis_m+1)));
    pop_vari_mutation  = max(min(pop_vari_mutation,upper),lower);
    % calculate the objective values of offsprings
    pop_obj_mutation = feval(fun_name, pop_vari_mutation, num_obj);
    % environment selection
    pop_vari_inter = [pop_vari;pop_vari_mutation];
    pop_obj_inter = [pop_obj;pop_obj_mutation];
    [front_rank,rank_num] = NonDominated_Rank(pop_obj_inter,pop_size);
    crowd_distance = Crowding_Distance(pop_obj_inter,front_rank);
    [B,index] = sortrows([front_rank,-crowd_distance]);
    pop_vari = pop_vari_inter(index(1:pop_size),:);
    pop_obj = pop_obj_inter(index(1:pop_size),:);
    evaluation = evaluation + pop_size;
    generation = generation + 1;
    non_dominated_front = Pareto_Set(pop_obj);
    IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
    if num_obj == 2
        scatter(pop_obj(:,1),pop_obj(:,2),'ro');hold on;
        scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
        scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
        title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        drawnow;hold off;
    elseif num_obj == 3
        scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
        scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
        scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
        title(sprintf('NSGA-II on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        view(135,30);drawnow;hold off;
    end
end


