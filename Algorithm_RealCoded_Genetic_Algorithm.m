% genetic algorithm
num_vari = 2;
fun_name = 'Ellipsoid';
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
pop_size = 20;
max_gen = 20;
fmin_record_GA = zeros(max_gen,1);
generation = 1;
% the initial generation
pop_size = pop_size - rem(pop_size,2);
pop_vari = lhsdesign(pop_size, num_vari).*(upper_bound - lower_bound) + lower_bound;
% calculate the objective values
figure(1);
contour(x1_mesh,x2_mesh,loss_mesh,50);
axis equal;
hold on;
scatter(pop_vari(:,1),pop_vari(:,2),'ro','filled');
title(sprintf('generation:%d',generation));
hold off;
drawnow;
pause(1);
pop_fitness = feval(fun_name,pop_vari);
fmin_record_GA(generation,1) = min(pop_fitness);
evaluation = pop_size;
% the evoluation of the generation
while generation <= max_gen
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
    figure(1);
    contour(x1_mesh,x2_mesh,loss_mesh,50);
    axis equal;
    hold on;
    scatter(pop_vari(:,1),pop_vari(:,2),'ro','filled');
    title(sprintf('generation:%d',generation));
    hold off;
    drawnow;
    pause(1);
    % update the evaluation number and generation number
    generation = generation + 1;
    fmin_record_GA(generation,1) = min(pop_fitness);
    fprintf('Genetic Algorithm on %d-D %s function, run: %d, evaluation: %d, fmin: %0.2f\n',num_vari,fun_name,1,evaluation,min(pop_fitness));
end

