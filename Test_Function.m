function [lower_bound,upper_bound] = Test_Function(fun_name,num_vari)

switch fun_name
    case 'Ellipsoid'
        lower_bound = -5.12*ones(1,num_vari);
        upper_bound = 5.12*ones(1,num_vari);
    case 'Rosenbrock'
        lower_bound = -2.048*ones(1,num_vari);
        upper_bound = 2.048*ones(1,num_vari);
    case 'Ackley'
        lower_bound = -32.768*ones(1,num_vari);
        upper_bound = 32.768*ones(1,num_vari);
    case 'Griewank'
        lower_bound = -600*ones(1,num_vari);
        upper_bound = 600*ones(1,num_vari);
    case 'Rastrigin'
        lower_bound = -5.12*ones(1,num_vari);
        upper_bound = 5.12*ones(1,num_vari);
    case 'Sixhump'
        lower_bound = -2*ones(1,num_vari);
        upper_bound = 2*ones(1,num_vari);
    case 'Regression_Loss'
        lower_bound = -10*ones(1,num_vari);
        upper_bound = 10*ones(1,num_vari);
    otherwise
        lower_bound = -10*ones(1,num_vari);
        upper_bound = 10*ones(1,num_vari);

end

