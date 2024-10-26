function [f,df] = Regression_Loss(w)
w = w';
X = [1,1;2,1;3,1;4,1;5,1];
y = [4.5;7.1;9.9;11.9;12.3];
f = sum((X*w - y).^2)/5;
df = 2 *X'*(X*w-y)/5;
df = df';

end

