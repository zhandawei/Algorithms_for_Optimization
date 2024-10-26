function [f,df] = Batch_Regression_Loss(w)
w = w';
X = [1,1;2,1;3,1;4,1;5,1];
y = [4.5;7.1;9.9;11.9;12.3];
rand_num = randperm(5,5);
X_rand = X(rand_num,:);
y_rand = y(rand_num,:);
f = sum((X_rand*w - y_rand).^2);
df = 2 *X_rand'*(X_rand*w-y_rand)/5;
df = df';
end
