function pareto_front = Calculate_Pareto_Front(obj_fun, num_n, num_obj)

switch obj_fun
    case 'ZDT1'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1 - f(:,1).^0.5;
        pareto_front = f;
    case 'ZDT2'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^2;
        pareto_front = f;
    case 'ZDT3'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^0.5-f(:,1).*sin(10*pi*f(:,1));
        pareto_front = f;
    case 'ZDT4'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^0.5;
        pareto_front = f;
    case 'ZDT6'
        minf1     = 0.280775;
        f(:,1)    = (minf1:(1-minf1)/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^2;
        pareto_front = f;
    case 'DTLZ1'
        f = UniformPoint(num_n,num_obj)/2;
        pareto_front = f;
    case 'DTLZ2'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ3'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ4'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ5'
        f = [0:1/(num_n-1):1;1:-1/(num_n-1):0]';
        f = f./repmat(sqrt(sum(f.^2,2)),1,size(f,2));
        f = [f(:,ones(1,num_obj-2)),f];
        pareto_front = f./sqrt(2).^repmat([num_obj-2,num_obj-2:-1:0],size(f,1),1);
    case 'DTLZ6'
        f = [0:1/(num_n-1):1;1:-1/(num_n-1):0]';
        f = f./repmat(sqrt(sum(f.^2,2)),1,size(f,2));
        f = [f(:,ones(1,num_obj-2)),f];
        pareto_front = f./sqrt(2).^repmat([num_obj-2,num_obj-2:-1:0],size(f,1),1);
    case 'DTLZ7'
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
%         X            = ReplicatePoint(10000,num_obj-1);
        if num_obj  > 2
            sample_num = (ceil(num_n^(1/(num_obj - 1))))^(num_obj - 1);
            Gap       = 0:1/(sample_num^(1/(num_obj - 1))-1):1;
            eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:(num_obj - 1))))
            eval(sprintf('X=[%s];',sprintf('c%d(:),',1:(num_obj - 1))))
        else
            X = (0:1/(num_n-1):1)';
        end
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        pareto_front            = [X,2*(num_obj-sum(X/2.*(1+sin(3*pi.*X)),2))];
end
end
function [W,N] = UniformPoint(N,M)
%UniformPoint - Generate a set of uniformly distributed points on the unit
%hyperplane
%
%   [W,N] = UniformPoint(N,M) returns approximate N uniformly distributed
%   points with M objectives.
%
%   Example:
%       [W,N] = UniformPoint(275,10)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    H1 = 1;
    while nchoosek(H1+M,M-1) <= N
        H1 = H1 + 1;
    end
    W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
    W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
    if H1 < M
        H2 = 0;
        while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
            H2 = H2 + 1;
        end
        if H2 > 0
            W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
            W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
            W  = [W;W2/2+1/(2*M)];
        end
    end
    W = max(W,1e-6);
    N = size(W,1);
end