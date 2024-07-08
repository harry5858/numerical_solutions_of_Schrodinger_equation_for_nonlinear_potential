clear;

%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = -1*w^2(x)*y----> h(x,y,z)
for nt=[240000 120000 60000]
tmin=0;
tmax=12;
dt=2*((tmax-tmin)/nt);
% t=zeros(1,nt);
% for i=1:nt
%     t(i)=i*dt;
% end
t=linspace(tmin,tmax,nt);
%  w initial and w final

w_i=3;
w_f=6;

a=(w_f+w_i)/2;
b=(w_f-w_i)/2;

dT=0.5;
% w=2*ones(1,nt);
w=a+(b*tanh((t-(tmax/2))/dT));
% figure(1)
% plot(t,w)

% RK 4
f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1)=1/sqrt(2*w_i);
g(1)=-1i*w_i*1/sqrt(2*w_i);
for n=1:2:nt-2
    
    k1=dt*(g(n));
    l1=dt*(-1*f(n)*w(n)^2);
    k2=dt*(g(n)+l1/2);
    l2=dt*(-1*(f(n)+k1/2)*w(n+1)^2);
    k3=dt*(g(n)+l2/2);
    l3=dt*(-1*(f(n)+k2/2)*w(n+1)^2);
    k4=dt*(g(n)+l3/2);
    l4=dt*(-1*(f(n)+k3/2)*w(n+2)^2);
    f(n+2)=f(n)+(1/6)*(k1+2*k2+2*k3+k4);
    g(n+2)=g(n)+(1/6)*(l1+2*l2+2*l3+l4);    
end
t1=linspace(tmin,tmax,nt/2);
f_new=zeros(1,nt/2);
f_new(1)=f(1);
for i=2:(nt-2)/2
f_new(i)=f(2*i+1);
end
f_new(nt/2)=f_new(nt/2-1);
f_ex=(1/2)*exp(-1i*2*t);

% plot(t,f_ex);
% hold on
% plot(t,f)

if nt==240000
    f1=f_new;
elseif nt==120000
    f2=f_new;
elseif nt==60000
    f4=f_new;
end
end
p_num=zeros(1,30000);
p_den=zeros(1,30000);
p=zeros(1,30000);
for j=1:30000
p_num(j) =abs( f4(j)- f2(2*j));
p_den(j) =abs( f2(2*j)- f1(4*j));
p(j)=p_num/p_den;
end
figure(1)
plot(t1,log2(p))

















