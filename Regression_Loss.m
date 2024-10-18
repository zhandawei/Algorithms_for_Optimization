function y = Regression_Loss(x)

y = 55*x(:,1).^2 + 5*x(:,2).^2 + 30*x(:,1).*x(:,2) -315*x(:,1) -91.4*x(:,2) +461.57;

end

