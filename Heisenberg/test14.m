clear;

%solve d2f/dt2+w(t)^2*f=0
%solve a system of equations
%rename f to y and t to x( not the space coordinate )
%system of equations
% dydx = z -------> g(x,y,z)
% dzdx = w^2(x)*y----> h(x,y,z)

tmin=0;
tmax=12;
nt=1200000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

%  w initial and w final  % sudden transition

w_i=3;
w_f=6;
w(1:nt/2)=w_i;
w(nt/2:nt)=w_f;

% RK 4
f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1) = sqrt(1/6);
g(1) = -1i*sqrt(3/2);

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

figure(1)
plot(t,real(f),'-b','Linewidth',1.5)
hold 
plot(t,imag(f),'-r','Linewidth',1.5)
xlabel('t')
title('real and imaginary part of f(t) vs t')
set(gca,'fontsize',20)
legend('real part of f(t)','imag part of f(t)')

alpha=zeros(1,nt);
beta=zeros(1,nt);
for i=2:nt-1

f_in_a=sqrt(1/(2*w_i))*exp(-1i*w_i*t(i));
f_in_b=conj(f_in_a);
df_in_a=-1i*w_i*f_in_a;
df_in_b=conj(df_in_a);

df_t= ((f(i+1))-(f(i-1)))/(2*dt);

alpha(i)=1i*(f_in_b*df_t-df_in_b*f(i));
beta(i)=1i*(f_in_a*df_t-df_in_a*f(i));

alpha(nt-1)=alpha(nt-2);
alpha(nt)=alpha(nt-1);
alpha(1)=alpha(2);
beta(nt-1)=beta(nt-2);
beta(nt)=beta(nt-1);
beta(1)=beta(2);
end


figure(2)

alpha_1=abs(alpha);

plot(t,alpha_1,'Linewidth',1.5)
title('absolute value of alpha vs t')
xlabel('t')
ylabel('abs(alpha)')
set(gca,'fontsize',20)

figure(3)

alpha_real=real(alpha);
alpha_imag=imag(alpha);

plot(t,alpha_real,'Linewidth',1.5)
hold on
plot(t,alpha_imag,'Linewidth',1.5)
title('real and imaginary part of alpha vs t')
xlabel('t')
set(gca,'fontsize',20)

figure(4)

beta_1=abs(beta);

plot(t,beta_1,'Linewidth',1.5)
title('absolute vlaue of beta vs t')
xlabel('t')
ylabel('abs(beta)')
set(gca,'fontsize',20)

figure(5)

beta_real=real(beta);
beta_imag=imag(beta);

plot(t,beta_real,'Linewidth',1.5)
hold on
plot(t,beta_imag,'Linewidth',1.5)
title('real and imaginary part of beta vs t')
xlabel('t')
set(gca,'fontsize',20)

figure(6)

beta_2=beta_1.^2;

plot(t,beta_2,'Linewidth',1.5)
title('out Number operator (abs(beta))^2 vs t')
xlabel('t')
ylabel('out Number operator')
set(gca,'fontsize',20)

g=zeros(1,nt);

for j=1:nt

   g(j)=abs(alpha(j))^2-abs(beta(j))^2;
   
end

figure(7)

plot(t,g,'Linewidth',1.5)
title('abs(alpha)^2 - abs(beta)^2 vs t')
xlabel('t')
ylabel('abs(alpha)^2 - abs(beta)^2')
set(gca,'fontsize',20)

% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=alpha(nt);
beta_p=beta(nt);

A=-1*conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0=(abs(N))^2;
C2=(abs(N*sqrt(2)*(A)))^2;
C4=(abs(N*sqrt(6)*(A^2)))^2;

n=0:1:4;

figure(6)
plot(n(1),C0,'b*','Markersize', 15)
hold on
plot(n(3),C2,'b*','Markersize', 15)
plot(n(5),C4,'b*','Markersize', 15)
axis([-0.1 5 0 1])
xlabel('n','Fontsize', 24)
ylabel('probability |Cn|^2','Fontsize', 24)
title('probability |<2n out | 0 in>|^2 (sudden transition) ','Fontsize', 20)

