clearvars;clc;close all;
load('optimization_result_f2.mat');

median_0 = prctile(fmin0,50,2);
lower_0 = prctile(fmin0,25,2);
upper_0 = prctile(fmin0,75,2);

median_1 = prctile(fmin1,50,2);
lower_1 = prctile(fmin1,25,2);
upper_1 = prctile(fmin1,75,2);

figure;
plot(200:1000,median_0,'b'); hold on;
plot(200:1000,median_1,'r');
h = fill([(200:1000)';(1000:-1:200)'],[lower_0;upper_0(end:-1:1)],'b');
set(h,'edgealpha',0,'facealpha',0.3);
h = fill([(200:1000)';(1000:-1:200)'],[lower_1;upper_1(end:-1:1)],'r');
set(h,'edgealpha',0,'facealpha',0.3);
xlabel('number of function evaluations');
ylabel('current minimum objective value');
legend('X', 'Y');

writematrix([(200:1000)',median_0,median_1],'result_f1.dat','delimiter','\t');
writematrix([[(200:1000)';(1000:-1:200)'],[lower_0;upper_0(end:-1:1)],...
    [lower_1;upper_1(end:-1:1)]],'result_f1_fill.dat','delimiter','\t');