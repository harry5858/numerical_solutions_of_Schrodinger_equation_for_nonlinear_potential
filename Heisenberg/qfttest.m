clear all;
%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = w^2(x)*y----> h(x,y,z)

tmin=0;
tmax=12;
nt=600000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end
%  set up for C(eta)

c_i=1;
c_f=4;
C(1:nt/2)=c_i;
C(nt/2:nt)=c_f;

% range of k

k=0:0.25:10;
n_i=numel(k);

beta2=zeros(1,n_i);
alpha2=zeros(1,n_i);
N=zeros(1,n_i);
for j=1:n_i

    w_i=sqrt(k(j)^2+c_i^2);
    w_f=sqrt(k(j)^2+c_f^2);
    
    w=sqrt(k(j)^2+C.^2);

% RK 4

f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1) = sqrt(1/2*w_i);
g(1) = -1i*w_i*sqrt(1/2*w_i);

for n=1:nt-2
    
    t(n+1) = t(n)+dt;
    k1=dt*(g(n));
    l1=dt*(-1*f(n)*(w(n))^2);
    k2=dt*(g(n)+l1/2);
    l2=dt*(-1*(f(n)+k1/2)*(w(n+1))^2);
    k3=dt*(g(n)+l2/2);
    l3=dt*(-1*(f(n)+k2/2)*(w(n+1))^2);
    k4=dt*(g(n)+l3/2); 
    l4=dt*(-1*(f(n)+k3/2)*(w(n+2))^2);
    f(n+1)=f(n)+(1/6)*(k1+2*k2+2*k3+k4);
    g(n+1)=g(n)+(1/6)*(l1+2*l2+2*l3+l4);
    f(nt)=f(nt-1);
    g(nt)=g(nt-1);
end

alpha=zeros(1,nt);
beta=zeros(1,nt);

for i=2:nt-1

f_out_a=sqrt(1/(2*w_f))*exp(-1i*w_f*t(i));
f_out_b=conj(f_out_a);
df_out_a=-1i*w_f*f_out_a;
df_out_b=1i*w_f*f_out_b;

df_in= (conj(f(i+1))-conj(f(i-1)))/(2*dt);

alpha(i)=1i*(conj(f(i))*df_out_a-df_in*f_out_a);
beta(i)=1i*((conj(f(i))*df_out_b-df_in*f_out_b));

alpha(nt-1)=alpha(nt-2);
alpha(nt)=alpha(nt-1);
alpha(1)=alpha(2);
beta(nt-1)=beta(nt-2);
beta(nt)=beta(nt-1);
beta(1)=beta(2);
end

beta2(j)=(abs(beta(nt-10000)))^2;
alpha2(j)=(abs(alpha(nt-10000)))^2;
N(j)=alpha2(j)-beta2(j);
end

figure(1)
plot(k,beta2,'*','Markersize',15)
xlabel('k')
ylabel('|beta(k)|^2')
title('Average number of particles in different modes')
set(gca,'fontsize',20)
figure(2)
plot(k,alpha2,'*','Markersize',15)
figure(3)
plot(k,N,'*','Markersize',15)

