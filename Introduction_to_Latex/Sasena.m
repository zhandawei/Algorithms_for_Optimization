function y=Sasena(x)
%----------------------------------------------------------
% Sasena function (called mystery function by Sasena
%
% Sasena, M. J., ¡°Flexibility and Efficiency Enhancements for
% Constrained Global Design Optimization with Kriging Approximations,¡±
% Ph.D. Thesis, Univ. of Michigan, Ann Arbor, MI, 2002.
%
% 0 <= x1 <= 5
% 0 <= x2 <= 5
% fmin = -1.4565
%----------------------------------------------------------

x1=x(:,1);x2=x(:,2);
y=2+0.01*(x2-x1.^2).^2+(1-x1).^2+2*(2-x2).^2+7*sin(0.5*x1).*sin(0.7*x1.*x2);


end