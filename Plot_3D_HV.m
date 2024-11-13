clearvars;clc;close all;
m = 3;
n = 100;
v1 = normrnd(0,1,[n,m]);
front = 10*abs(v1)./vecnorm(v1,2,2);
ref_point = 10*ones(1,m);
figure;
scatter3(ref_point(:,1),ref_point(:,2),ref_point(:,3),'bo','filled');hold on;
scatter3(front(:,1),front(:,2),front(:,3),'ro','filled');
for ii = 1:n
    plot_cube(front(ii,:),ref_point - front(ii,:),ii);
end

% Walking Fishing Group method
hv = Hypervolume(front,ref_point);
% Monte Carlo method
N = 1000;
lower_bound = min(front);
A = prod(ref_point - lower_bound);
rand_sample = rand(N,m).*(ref_point-lower_bound)+lower_bound;
is_dominated = zeros(N,n);
for ii = 1:n
    is_dominated(:,ii) = sum(rand_sample >=  front(ii,:),2) == m;
end
k = sum(sum(is_dominated,2)~=0);
a = A*k/N;
disp([a,hv]);
















function plot_cube(origin,length,color)
% Define the vertexes of the unit cubic
ver = [1 1 0;
    0 1 0;
    0 1 1;
    1 1 1;
    0 0 1;
    1 0 1;
    1 0 0;
    0 0 0];

%  Define the faces of the unit cubic
fac = [1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8];
cube = [ver(:,1)*length(1),ver(:,2)*length(2),ver(:,3)*length(3)];
cube = [cube(:,1)+origin(1),cube(:,2)+origin(2),cube(:,3)+origin(3)];
patch('Faces',fac,'Vertices',cube,'FaceVertexCData',color,'FaceColor','flat');
end