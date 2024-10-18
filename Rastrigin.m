function f = Rastrigin(x)
% the Rastrigin function
% xi = [-5.12,5.12]
d = size(x,2);
f = 10*d + sum(x.^2 - 10*cos(2*pi*x),2);

end