clear;

tmin=0;
tmax=12;
nt=120000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

c_i=1;
c_f=5;
a=(c_f+c_i)/2;
b=(c_f-c_i)/2;
dT=1;
c=a+b*tanh((t-6)/dT);

k=0:0.5:10;
n_i=numel(k);
beta2=zeros(1,n_i);

for n=1:n_i
    
w=sqrt(k(n)^2+c.^2);
w_in=sqrt(k(n)^2+c_i^2);
w_out=sqrt(k(n)^2+c_f^2);
% figure(1)
% plot(t,w)

% setup for solutions
f=zeros(1,nt);
g=zeros(1,nt);

%initial condition
f(1)=sqrt(1/(2*w_in));
g(1)=-1i*w_in*sqrt(1/(2*w_in));

%initial run
f(2)=f(1)-dt*((w(1)^2)*g(1));
g(2)=g(1)+dt*f(1);

for i=2:nt-1
    
    f(i+1)=f(i-1)-dt*((w(i)^2)*g(i));
    g(i+1)=g(i-1)+dt*f(i);

    g(nt)=g(nt-1);
    f(nt)=f(nt-1);
    
end
% figure(2)
% plot(t,real(g))

beta=zeros(1,nt);

for j=1:nt-1
    f_out=sqrt(1/(2*w_out))*exp(1i*w_out*t(j));
    df_out=1i*w_out*f_out;
    
    dg=(g(j+1)-g(j))/dt;
    
    beta(i)=1i*(g(i)*df_out-dg*f_out);
    
    beta(nt)=beta(nt-1);
end
% figure(3)
%     plot(t,abs(beta))

beta2(n)=(abs(beta(nt-10000)))^2;

end

plot(k,beta2,'*','Markersize',15)



