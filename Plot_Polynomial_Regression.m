clearvars;close all;
num_sample = 9;
sample_x = linspace(0,1,num_sample)';
sample_y = 5*sample_x.^2+4*sample_x-5 + (rand(num_sample,1)*2-1);

X1 = [sample_x.^2,sample_x,ones(num_sample,1)];
y = sample_y;
w1 = X1\y;
loss1 = sum((X1*w1-y).^2);

X2 = [sample_x.^8,sample_x.^7,sample_x.^6,sample_x.^5,sample_x.^4,sample_x.^3,sample_x.^2,sample_x,ones(num_sample,1)];
y = sample_y;
w2 = X2\y;
loss2 = sum((X2*w2-y).^2);


num_plot = 100;
plot_x = linspace(0,1,num_plot)';
plot_y1 = [plot_x.^2,plot_x,ones(num_plot,1)]*w1;
plot_y2 = [plot_x.^8,plot_x.^7,plot_x.^6,plot_x.^5,plot_x.^4,plot_x.^3,plot_x.^2,plot_x,ones(num_plot,1)]*w2;
figure;
scatter(sample_x,sample_y,'filled'); hold on;
plot(plot_x,plot_y1,'b-','Linewidth',1);
plot(plot_x,plot_y2,'r-','Linewidth',1);


