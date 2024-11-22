clearvars; close all;
fun_name = 'DTLZ2';
num_vari = 10;
num_obj = 3;
num_weight = 20;
max_gen = 5;
num_neighbor = 10;
lower_bound = zeros(1,num_vari);
upper_bound = ones(1,num_vari);
[weight,num_weight]= UniformPoint(num_weight,num_obj);
neighbor = pdist2(weight,weight);
[~,neighbor] = sort(neighbor,2);
neighbor = neighbor(:,1:num_neighbor);
% the first population
generation = 1;
evaluation = num_weight;
pop_vari = lhsdesign(num_weight, num_vari).*(upper_bound - lower_bound) + lower_bound;
pop_obj = feval(fun_name,pop_vari,num_obj);
ideal_point = min(pop_obj,[],1);
non_dominated_front = Pareto_Set(pop_obj);
pareto_front = Calculate_Pareto_Front(fun_name, 10000, num_obj);
pareto_front_plot = Calculate_Pareto_Front(fun_name, 1000, num_obj);
IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
    figure;hold on;
if num_obj == 2
    scatter(pop_obj(:,1),pop_obj(:,2),'ro');
    scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
    scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
    title(sprintf('MOEA/D on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    drawnow;hold off;
elseif num_obj == 3
    scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro');
    scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
    scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
    title(sprintf('MODA/D  on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    view(135,30);drawnow;hold off;
end

while generation < max_gen
    for ii = 1:num_weight
        this_neighbor = neighbor(ii,:);
        parent = pop_vari(this_neighbor(randperm(num_neighbor,2)),:);
        % simulated binary crossover
        dis_c = 20;
        mu  = rand(1,num_vari);
        parent1 = parent(1,:);
        parent2 = parent(2,:);
        beta = 1 + 2*min(min(parent1,parent2)-lower_bound,upper_bound-max(parent1,parent2))./max(abs(parent2-parent1),1E-6);
        alpha = 2 - beta.^(-dis_c-1);
        betaq = (alpha.*mu).^(1/(dis_c+1)).*(mu <= 1./alpha) + (1./(2-alpha.*mu)).^(1/(dis_c+1)).*(mu > 1./alpha);
        % the crossover is performed randomly on each variable
        betaq = betaq.*(-1).^randi([0,1],1,num_vari);
        offspring1 = 0.5*((1+betaq).*parent1 + (1-betaq).*parent2);
        offspring2 = 0.5*((1-betaq).*parent1 + (1+betaq).*parent2);
        crossover = [offspring1;offspring2];
        % mutation (ploynomial mutation)
        % dis_m is the distribution index of polynomial mutation
        dis_m = 20;
        pro_m = 1/num_vari;
        rand_var = rand(2,num_vari);
        mu  = rand(2,num_vari);
        deta = min(crossover-lower_bound, upper_bound-crossover)./(upper_bound-lower_bound);
        detaq = zeros(2,num_vari);
        position1 = rand_var<=pro_m & mu<=0.5;
        position2 = rand_var<=pro_m & mu>0.5;
        detaq(position1) = ((2*mu(position1) + (1-2*mu(position1)).*(1-deta(position1)).^(dis_m+1)).^(1/(dis_m+1))-1);
        detaq(position2) = (1 - (2*(1-mu(position2))+2*(mu(position2)-0.5).*(1-deta(position2)).^(dis_m+1)).^(1/(dis_m+1)));
        mutation = crossover + detaq.*(upper_bound-lower_bound);
        mutation  = max(min(mutation,upper_bound),lower_bound);
        offspring = mutation(1,:);
        offspring_obj = feval(fun_name,offspring,num_obj);
        % update the ideal point
        ideal_point =  min([pop_obj;offspring_obj],[],1);
        % update the populatuion
        g_old = max(abs(pop_obj(this_neighbor,:)-repmat(ideal_point,num_neighbor,1)).*weight(this_neighbor,:),[],2);
        g_new = max(repmat(abs(offspring_obj-ideal_point),num_neighbor,1).*weight(this_neighbor,:),[],2);
        pop_vari(this_neighbor(g_new<=g_old),:) = repmat(offspring,sum(g_new<=g_old),1);
        pop_obj(this_neighbor(g_new<=g_old),:) = repmat(offspring_obj,sum(g_new<=g_old),1);
    end
    generation = generation+1;
    evaluation = evaluation + num_weight;
    non_dominated_front = Pareto_Set(pop_obj);
    IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
    if num_obj == 2
        scatter(pop_obj(:,1),pop_obj(:,2),'ro');
        scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
        scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
        title(sprintf('MODA/D on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        drawnow;hold off;
    elseif num_obj == 3
        scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
        scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
        scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
        title(sprintf('MODA/D  on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        view(135,30);drawnow;hold off;
    end
end