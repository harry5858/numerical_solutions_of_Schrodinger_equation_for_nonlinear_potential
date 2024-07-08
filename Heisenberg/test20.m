trange=linspace(0,12,60000);
dt=max(trange)/60000;
nt=numel(trange);
w_f=6;
w_i=3;
dT=1;
[t, state_values] = ode45(@(t,x) myfunc(t,x,w_f,w_i,dT),trange,[sqrt(1/(2*w_i)),-1i*(w_i)*sqrt(1/(2*w_i))]);
f = state_values(:,1);
fdot = state_values(:,2);
figure(1)
plot(t,real(f))
hold on
plot(t,imag(f))
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

figure(2)

alpha_1=abs(alpha);
beta_1=abs(beta);

plot(t,alpha_1,'-b','Linewidth',1.5)
hold on
plot(t,beta_1,'-r','Linewidth',1.5)
title('absolute value of alpha and beta vs t')
xlabel('t')
set(gca,'fontsize',20)
legend('abs(Alpha)','abs(Beta)')

figure(3)

alpha_real=real(alpha);
alpha_imag=imag(alpha);

plot(t,alpha_real,'Linewidth',1.5)
hold on
plot(t,alpha_imag,'Linewidth',1.5)
title('real and imaginary part of alpha vs t')
xlabel('t')
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
% 
alpha_p=alpha_1(nt-6000);
beta_p=beta_1(nt-6000);

A=conj(beta_p)/(2*alpha_p);
N=sqrt(1/abs(alpha_p));

C0=(abs(N))^2;
C2=(abs(N*sqrt(2)*(A)))^2;
C4=(abs(N*sqrt(6)*(A^2)))^2;

n=0:1:4;

figure(8)
plot(n(1),C0,'b*','Markersize', 15)
hold on
plot(n(3),C2,'b*','Markersize', 15)
plot(n(5),C4,'b*','Markersize', 15)
axis([-0.1 5 0 1])
