%Numerov method
clear all;

tmin=0;
tmax=24;
nt=120000;
dt=(tmax-tmin)/nt;
t=zeros(1,nt);
for i=1:nt
    t(i)=i*dt;
end

%  w initial and w final

w_i=3;
w_f=6;

dT=0.1:0.1:2;
n_i=numel(dT);

C0=zeros(1,n_i);
C2=zeros(1,n_i);
C4=zeros(1,n_i);
Norm=zeros(1,n_i);
for j=1:n_i
    
a=(w_f+w_i)/2;
b=(w_f-w_i)/2;

w=a+(b*tanh((t-(tmax/2))/dT(j)));
% RK 4

f=zeros(1,nt);
g=zeros(1,nt);

% initial condition
f(1)=1/sqrt(2*w_i);
g(1)=-1i*w_i*1/sqrt(2*w_i);

%RK 2 for another initial condition

t(2) = t(1)+dt;
    k1=dt*(g(1));
    l1=dt*(-1*f(1)*(w(1)^2));
    k2=dt*(g(1)+l1/2);
    l2=dt*(-1*(f(1)+k1/2)*(w(2)^2));
    f(2)=f(1)+0.5*(k1+k2);

% new set of initial conditions is f(1) and f(2)

for n=2:nt-1
    
f(n+1) = (2*(1-(5/12)*(dt^2)*(w(n)^2))*f(n)-(1+(1/12)*(dt^2)*(w(n-1)^2))*f(n-1))/(1+(1/12)*(dt^2)*(w(n+1)^2));

end

f(nt)=f(nt-1);

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

Norm(j)=abs(alpha(nt))^2-abs(beta(nt))^2;
   
% Proability
% Only consider the first 3 even terms for comparison with the Schrodinger
% numerical solution

alpha_p=alpha(nt);
beta_p=beta(nt);

A=-1*conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0(j)=(abs(N))^2;
C2(j)=(abs(N*sqrt(2)*(A)))^2;
C4(j)=(abs(N*sqrt(6)*(A^2)))^2;

end
figure(1)
semilogy(dT,C0,'b*','Markersize', 15)
hold on
semilogy(dT,C2,'r')
semilogy(dT,C4,'k')
% axis([-5 20 0 1])
legend('C0','C2','C4')
% legend(' n=2 ',' n=4 ')
xlabel('\Delta T','Fontsize', 24)
title('Transition probability |<2n out | 0 in>|^2 , n=0,1,2 ','Fontsize', 20)
set(gca,'fontsize',20)
figure(2)
plot(dT,Norm)