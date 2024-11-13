function ref_point =  Calculate_Reference_Point(fun_name,num_obj)
switch fun_name
    case 'DTLZ1'
        ref_point = 400*ones(1,num_obj);
    case 'DTLZ2'
        ref_point = 2.5*ones(1,num_obj);
    case 'DTLZ3'
        ref_point = 1500*ones(1,num_obj);
    case 'DTLZ4'
        ref_point = 2.5*ones(1,num_obj);
    case 'DTLZ5'
        ref_point = 11*ones(1,num_obj);
    case 'DTLZ6'
        ref_point = 11*ones(1,num_obj);
    case 'DTLZ7'
        ref_point = [ones(1,num_obj-1),70];
end