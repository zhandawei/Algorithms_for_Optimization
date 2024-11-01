function [f,df,Hf] = Sixhump(x)
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
f =(4-2.1.*x1.^2+x1.^4./3).*x1.^2+x1.*x2+(-4+4.*x2.^2).*x2.^2;        
df = [2*x1.^5-8.2*x1.^3+8*x1+x2,16*x2.^3-8*x2+x1];
Hf = [10*x1.^4 - 24.6*x1.^2+8,1;1,48*x2.^2-8];

