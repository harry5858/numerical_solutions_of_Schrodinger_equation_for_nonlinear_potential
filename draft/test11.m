w_i=3;
w_f=6;
T=tmax/2;
a=(w_f+w_i)/2;
b=(w_f-w_i)/2;

dT=0.1:0.1:1.2;
n_i=numel(dT);
for j=1:n_i
    w=a+b.*(tanh((t-T)./dT(j)));
    figure(1)
    plot(t,w,'b')
    hold on
    axis([0 tmax 0 11])
    drawnow;
end






