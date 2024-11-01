clearvars;clc; format shortG;
num_vari = 16;
pop_size = 1000;
max_gen = 200;
% initialization
fmin_record = zeros(max_gen,1);
generation = 1;
chrom_length = 4*num_vari;
pop_chrom = randi([0,1],pop_size,chrom_length);
% decode
pop_x = zeros(pop_size,num_vari);
for ii = 1:num_vari
    pop_x(:,ii) = pop_chrom(:,4*(ii-1)+1:4*ii)*[8,4,2,1]' + 1;
end
pop_y = test_fun(pop_x);
fmin_record(1) = min(pop_y);
fprintf('generation: %d, fmin: %0.2f\n',generation,min(pop_y));
for generation  = 2:max_gen
    % parent selection
    parent_chrom = zeros(pop_size,chrom_length);
    for ii = 1:pop_size
        compete_index = randi([1,pop_size],1,2);
        [min_value,min_index] = min(pop_y(compete_index));
        select_index = compete_index(min_index);
        parent_chrom(ii,:) = pop_chrom(select_index,:);
    end
    % crossover
    cross_chrom = zeros(pop_size,chrom_length);
    for ii = 1:pop_size/2
        chrom_1 = parent_chrom(2*(ii-1)+1,:);
        chrom_2 = parent_chrom(2*ii,:);
        cross_point = randi([1,chrom_length-1]);
        cross_1 = [chrom_1(1:cross_point),chrom_2(cross_point+1:end)];
        cross_2 = [chrom_2(1:cross_point),chrom_1(cross_point+1:end)];
        cross_chrom(2*(ii-1)+1,:) = cross_1;
        cross_chrom(2*ii,:) = cross_2;
    end
    % mutation
    mutation_chrom = cross_chrom;
    for ii = 1:pop_size
        mutation_point = randi([1,chrom_length]);
        mutation_chrom(ii,mutation_point) = mod(mutation_chrom(ii,mutation_point)+1,2);
    end
    % decode
    mutation_x = zeros(pop_size,num_vari);
    for ii = 1:num_vari
        mutation_x(:,ii) = mutation_chrom(:,4*(ii-1)+1:4*ii)*[8,4,2,1]' + 1;
    end
    mutation_y = test_fun(mutation_x);
    % environmental selection
    combine_y = [pop_y;mutation_y];
    combine_x = [pop_x;mutation_x];
    combine_chrom = [pop_chrom;mutation_chrom];
    [sort_y,sort_index] = sort(combine_y);
    pop_chrom = combine_chrom(sort_index(1:pop_size),:);
    pop_x = combine_x(sort_index(1:pop_size),:);
    pop_y = combine_y(sort_index(1:pop_size),:);
    min_y  = pop_y(1);
    min_x = pop_x(1,:);
    fmin_record(generation) = min_y;
    fprintf('generation: %d, fmin: %0.2f\n',generation,min(pop_y));
end

figure;
plot(fmin_record);



function y = test_fun(x)
y = sum((x - (1:16)).^(2:2:32),2);
end








