function [f,df,Hf] = Valley_Function(x)
f = 0.1*(x(:,1)-1).^2 + 2*(x(:,2)-2).^2;
df = [0.2*(x(:,1)-1),4*(x(:,2)-2)];
Hf = [0.2,0;0,4];


end

