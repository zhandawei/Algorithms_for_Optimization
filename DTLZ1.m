function y=DTLZ1(x,m)
% 0<xi<1   i=1,2...6
% x is a Vector
n=size(x,2);     %number of design variable
k=n-m+1;
g=100*(k+sum((x(:,m:end)-0.5).^2-cos(20*pi*(x(:,m:end)-0.5)),2));

y(:,1)=0.5*prod(x(:,1:m-1),2).*(1+g);
if m>2
    for ii=2:m-1
        y(:,ii)=0.5*prod(x(:,1:m-ii),2).*(1-x(:,m-ii+1)).*(1+g);
    end
end
y(:,m)=0.5*(1-x(:,1)).*(1+g);


end