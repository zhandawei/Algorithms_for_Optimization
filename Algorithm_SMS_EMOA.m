clearvars;clc;close all;
fun_name = 'DTLZ5';
num_obj = 3;
num_vari = 10;
pop_size = 50;
max_evalution = 1000;
pareto_front = Calculate_Pareto_Front(fun_name,10000,num_obj);
pareto_front_plot = Calculate_Pareto_Front(fun_name,1000,num_obj);
lower_bound = zeros(1,num_vari);
upper_bound = ones(1,num_vari);
pop_vari = rand(pop_size,num_vari).*(upper_bound-lower_bound) + lower_bound;
pop_obj = feval(fun_name,pop_vari,num_obj);
evaluation = pop_size;
generation = 1;
non_dominated_front = Pareto_Set(pop_obj);
IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
if num_obj == 2
    scatter(pop_obj(:,1),pop_obj(:,2),'ro');hold on;
    scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
    scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
    title(sprintf('SMS-EMOA on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    drawnow;hold off;
elseif num_obj == 3
    scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
    scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
    scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
    title(sprintf('SMS-EMOA on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
    view(135,30);drawnow;hold off;
end
while evaluation < max_evalution
    % parent selection (tournament selection)
    front_rank = NonDominated_Rank(pop_obj,num_obj);
    rank_1 = find(front_rank==1);
    select_index = randi(size(rank_1,1),2);
    parent1 = pop_vari(rank_1(select_index(1)),:);
    parent2 = pop_vari(rank_1(select_index(2)),:);
    % crossover (simulated binary crossover)
    dis_c = 20;
    mu  = rand(1,num_vari);
    beta = 1 + 2*min(min(parent1,parent2)-lower_bound,upper_bound-max(parent1,parent2))./max(abs(parent2-parent1),1E-6);
    alpha = 2 - beta.^(-dis_c-1);
    betaq = (alpha.*mu).^(1/(dis_c+1)).*(mu <= 1./alpha) + (1./(2-alpha.*mu)).^(1/(dis_c+1)).*(mu > 1./alpha);
    betaq = betaq.*(-1).^randi([0,1],1,num_vari);
    betaq(rand(1,num_vari)>0.5) = 1;
    offspring1 = 0.5*((1+betaq).*parent1 + (1-betaq).*parent2);
    offspring2 = 0.5*((1-betaq).*parent1 + (1+betaq).*parent2);
    pop_crossover = [offspring1;offspring2];
    % mutation (ploynomial mutation)
    dis_m = 20;
    pro_m = 1/num_vari;
    rand_var = rand(2,num_vari);
    mu  = rand(2,num_vari);
    deta = min(pop_crossover-lower_bound, upper_bound-pop_crossover)./(upper_bound-lower_bound);
    detaq = zeros(2,num_vari);
    position1 = rand_var<=pro_m & mu<=0.5;
    position2 = rand_var<=pro_m & mu>0.5;
    detaq(position1) = ((2*mu(position1) + (1-2*mu(position1)).*(1-deta(position1)).^(dis_m+1)).^(1/(dis_m+1))-1);
    detaq(position2) = (1 - (2*(1-mu(position2))+2*(mu(position2)-0.5).*(1-deta(position2)).^(dis_m+1)).^(1/(dis_m+1)));
    pop_mutation = pop_crossover + detaq.*(upper_bound-lower_bound);
    pop_mutation  = max(min(pop_mutation,upper_bound),lower_bound);
    temp = randi(2);
    x = pop_mutation(temp,:);
    % calculate the objective values of offsprings
    y = feval(fun_name,x,num_obj);
    pop_vari = [pop_vari;x];
    pop_obj = [pop_obj;y];
    % discard the worst individual
    front_rank = NonDominated_Rank(pop_obj,num_obj);
    rank_v = max(front_rank);
    if rank_v >1
        select_index = find(front_rank == rank_v);
        R_v = pop_obj(select_index,:);
        number_of_dominated_points = zeros(size(R_v,1),2);
        number_of_dominated_points(:,1) = select_index;
        for i = 1:size(R_v,1)
            obj = repmat(R_v(i,:),size(pop_obj,1),1);
            C = (obj >= pop_obj);
            A = find(sum(C,2)==num_obj);
            number_of_dominated_points(i,2) = size(A,1);
        end
        [~,remove_index] = max(number_of_dominated_points(:,2));
        remove = number_of_dominated_points(remove_index,1);
    else
        R_1 = zeros(size(pop_obj,1),1);
        ref_p = max(pop_obj,[],1) + 1;
        HVsum = Hypervolume(pop_obj,ref_p);
        for i = 1: size(pop_obj,1)
            temp = pop_obj;
            temp(i,:) = [];
            R_1(i) = HVsum - Hypervolume(temp,ref_p);
        end
        [~,remove] = min(R_1);
    end
    pop_vari(remove,:) = [];
    pop_obj(remove,:) = [];
    evaluation = evaluation + 1;
    generation = generation + 1;
    non_dominated_front = Pareto_Set(pop_obj);
    IGD = mean(min(pdist2(pareto_front,non_dominated_front),[],2));
    if num_obj == 2
        scatter(pop_obj(:,1),pop_obj(:,2),'ro');hold on;
        scatter(non_dominated_front(:,1),non_dominated_front(:,2),'ro','filled');
        scatter(pareto_front_plot(:,1),pareto_front_plot(:,2),'b.');
        title(sprintf('SMS-EMOA on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        drawnow;hold off;
    elseif num_obj == 3
        scatter3(pop_obj(:,1),pop_obj(:,2),pop_obj(:,3),'ro'); hold on;
        scatter3(non_dominated_front(:,1),non_dominated_front(:,2),non_dominated_front(:,3),'ro','filled');
        scatter3(pareto_front_plot(:,1),pareto_front_plot(:,2),pareto_front_plot(:,3),'b.');
        title(sprintf('SMS-EMOA on %d-objective %s \n generation: %d, evaluations: %d, IGD: %0.4g',num_obj,fun_name,generation,evaluation,IGD));
        view(135,30);drawnow;hold off;
    end
end



