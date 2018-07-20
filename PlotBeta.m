par.Td=20;
par.tau=8*24;
par.lambda=32;

x=zeros(1,481);
y=zeros(1,481);
for t=1:481
    x(t)=t-1;
    y(t)=beta(t-1,par);
end
figure(10); hold on;
plot(x,y)   