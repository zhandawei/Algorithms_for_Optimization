clearvars; close all;
N = 10000;
mu1 = 5; 
mu2 = 10;
sigma1 = 10; 
sigma2 = 10;
X1 = normrnd(mu1,sigma1,[N,1]);
X2 = normrnd(mu2,sigma2,[N,1]);
X = max(X1,X2);
figure(1);
histogram(X1,100); hold on;
histogram(X2,100); 
histogram(X,100); 
legend('X1','X2','max(X1,X2)')
disp(mean(X))