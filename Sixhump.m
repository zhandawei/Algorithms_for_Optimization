function value = Sixhump(x)
%--------------------------------------------------------
% 6-Hump Camel Test Function for Nonlinear Optimization
%
% Taken from "Towards Global Optimisation 2",edited by L.C.W. Dixon and G.P.
% Szego, North-Holland Publishing Company, 1978. ISBN 0 444 85171 2
%
% -2 <= x1 <= 2                  
% -2 <= x2 <= 2                  
% fmin = -1.0316284535
% xmin = 0.08984201  -0.08984201
%       -0.71265640   0.71265640
%--------------------------------------------------------

%---------------------------------------------------------%
% http://www4.ncsu.edu/~definkel/research/index.html  
%---------------------------------------------------------%
x1 = x(:,1);
x2 = x(:,2);
value =(4-2.1.*x1.^2+x1.^4./3).*x1.^2+x1.*x2+(-4+4.*x2.^2).*x2.^2;        

